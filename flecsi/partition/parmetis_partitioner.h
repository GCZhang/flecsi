/*~-------------------------------------------------------------------------~~*
 *  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
 * /@@/////  /@@          @@////@@ @@////// /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
 * //       ///  //////   //////  ////////  //
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/

#ifndef FLECSI_PARMETIS_PARTITIONER_H
#define FLECSI_PARMETIS_PARTITIONER_H

#include <algorithm>
#include <numeric>

#include <cinchtest.h>
#include <parmetis.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/gather.h>
#include <thrust/transform.h>
#include <thrust/functional.h>

#include "flecsi/topology/mesh_topology.h"
#include "flecsi/partition/index_partition.h"

namespace flecsi {

using namespace flecsi::topology;
using namespace flecsi::dmp;

class parmetis_partitioner {
public:
  // TODO: incompatibility between ParMetis' idx_t and FleCSI's id_t or size_t
  using mesh_graph_partitions = std::vector<mesh_graph_partition<idx_t>>;
  using index_partition_t = index_partition__<idx_t>;

  parmetis_partitioner(MPI_Comm _mpi_comm = MPI_COMM_WORLD) :
    mpi_comm(_mpi_comm) {
    MPI_Comm_size(mpi_comm, &comm_size);
    MPI_Comm_rank(mpi_comm, &rank);
  }

  // Input parameter graph_partitions can not be made const because ParMetis
  // is not const correct.
  void partition(mesh_graph_partition<idx_t> &graph_partition,
                 index_partition_t &index_partition) {
    // input parameters for ParMetis
    idx_t wgtflag = 0;  // no weight
    idx_t numflag = 0;  // C/C++ style index from zero
    idx_t ncon = 1;
    idx_t nparts = comm_size;   // partition into # of ranks pieces

    // the fraction of "vertices" in each piece.
    std::vector<real_t> tpwgts(comm_size, 1.0/comm_size);

    // No, these arrays can not be NULL
    real_t ubvec[1] = {1.05};
    idx_t options[1] = {0};

    idx_t edgecut;      // number of edge cut returned from ParMetis

    // Reserve memory of number of cells on this rank for output from ParMetis
    int number_of_cells =
      graph_partition.partition[rank+1] -
      graph_partition.partition[rank];

    std::vector<idx_t> part(number_of_cells);

    // FIXME: How does ParMetis find out how many vertices are there on each node?
    auto ret = ParMETIS_V3_PartKway(
      graph_partition.partition.data(),  // vtxdist
      graph_partition.offset.data(),     // xadj
      graph_partition.index.data(),      // adjncy
      nullptr, /* vwgt */
      nullptr, /* adjwgt */
      &wgtflag,
      &numflag,
      &ncon,
      &nparts,
      tpwgts.data(),
      ubvec,
      options,
      &edgecut,
      part.data(),
      &mpi_comm
    );

    redistribute(part, index_partition);
  }

  void partition(mesh_graph_partitions &graph_partitions,
                 index_partition_t &index_partition) {
    partition(graph_partitions[rank], index_partition);
  }

  // TODO: should we make partition a reference? Do we want to reuse partition
  // in the caller after redistribute() mangle with it?
  void redistribute(std::vector<idx_t> partition,
                    index_partition_t &index_partition) {
    // number of cells on this rank before redistribution.
    int init_num_cells = partition.size();

    // we have to initialized cell_id_start to 0 since MPI_Excan does not update
    // it on Rank 0.
    // TODO: is this just the vtxdist in ParMetis?
    int cell_id_start = 0;
    // Exclusive scan to find out the start cell entity id on this rank.
    // This assumes:
    //     1. cell entity id starts from 0
    //     2. cells are numbered consecutively in the initial partitioning DCSR
    //        as input to ParMetis.
    MPI_Exscan(&init_num_cells, &cell_id_start, 1, MPI_INT, MPI_SUM, mpi_comm);

    // actually enumerate "global" cell entity ids
    std::vector<idx_t> cell_ids(init_num_cells);
    std::iota(cell_ids.begin(), cell_ids.end(), cell_id_start);

    // sort cell ids according to their destination partitions
    thrust::stable_sort_by_key(partition.begin(), partition.end(), cell_ids.begin());

    // calculate pairs of (destination, number of cells).
    std::vector<int> destinations(comm_size);
    std::iota(destinations.begin(), destinations.end(), 0);

    std::vector<int> send_counts(comm_size);

    std::vector<int> lbs(comm_size);
    thrust::lower_bound(partition.begin(), partition.end(),
                        thrust::counting_iterator<int>(0),
                        thrust::counting_iterator<int>(0)+ comm_size,
                        lbs.begin());
    std::vector<int> ubs(comm_size);
    thrust::upper_bound(partition.begin(), partition.end(),
                        thrust::counting_iterator<int>(0),
                        thrust::counting_iterator<int>(0)+ comm_size,
                        ubs.begin());

    thrust::transform(ubs.begin(), ubs.end(),
                      lbs.begin(), send_counts.begin(),
                      thrust::minus<int>());

    // All to all communication to tell others how many cells are coming.
    std::vector<int> recv_counts(comm_size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                 recv_counts.data(), 1, MPI_INT,
                 MPI_COMM_WORLD);

    // calculate send and receive displacements from send and receive counts
    std::vector<int> send_disp(comm_size);
    thrust::exclusive_scan(send_counts.begin(), send_counts.end(),
                           send_disp.begin());

    std::vector<int> recv_disp(comm_size);
    thrust::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                           recv_disp.begin());

    // total number of cell received. This is also the total number of cells
    // on this rank after redistribution.
    int total_recv = thrust::reduce(recv_counts.begin(), recv_counts.end());

    // All to all communication to exchange cell ids.
    std::vector<int> my_cell_ids(total_recv);
    MPI_Alltoallv(cell_ids.data(),    send_counts.data(), send_disp.data(), MPI_INT,
                  my_cell_ids.data(), recv_counts.data(), recv_disp.data(), MPI_INT,
                  MPI_COMM_WORLD);

    index_partition.exclusive = std::move(my_cell_ids);
  }
private:
  MPI_Comm mpi_comm;
  int comm_size;
  int rank;
};
}
#endif //FLECSI_PARMETIS_PARTITIONER_H
