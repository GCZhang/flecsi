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

#include <algorithm>
#include <numeric>

#include <cinchtest.h>
#include <parmetis.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CPP
#include <thrust/iterator/constant_iterator.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/gather.h>

#include "../dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace testing;

class mpi_parmetis_2way : public Test {
protected:

  virtual void SetUp() {
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto cells_per_rank = dolfin.num_cells()/comm_size;
    for (auto i = 0; i < comm_size - 1 ; i++) {
      cell_sizes.push_back(cells_per_rank);
    }
    cell_sizes.push_back(dolfin.num_cells() - cells_per_rank*(comm_size-1));

    // create Distributed CSR from dolfin's cell to cell connectivity.
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);

    // input parameters for ParMetis
    idx_t wgtflag = 0;  // no weight
    idx_t numflag = 0;  // no flag
    idx_t ncon = 1;
    idx_t nparts = 2;   // partition into two pieces

    // No, these three arrays can not be NULL
    real_t tpwgts[2] = {0.5, 0.5};
    real_t ubvec[1] = {1.05};
    idx_t options[1] = {0};

    idx_t edgecut;      // number of edge cut returned from ParMetis
    MPI_Comm comm = MPI_COMM_WORLD;

    // Use resize() to reserve memory for ParMetis output and also move end().
    part.resize(5);

    // FIXME: How does ParMetis find out how many vertices are there on each node?
    auto ret = ParMETIS_V3_PartKway(
      cell_partitions[rank].partition.data(),  // vtxdist
      cell_partitions[rank].offset.data(),     // xadj
      cell_partitions[rank].index.data(),      // adjncy
      nullptr, /* vwgt */
      nullptr, /* adjwgt */
      &wgtflag,
      &numflag,
      &ncon,
      &nparts,
      tpwgts,
      ubvec,
      options,
      &edgecut,
      part.data(),
      &comm
    );

    ASSERT_EQ(ret, METIS_OK);
  }

  // WARNING: we deliberately use Metis' idx_t (which is either 32 or 64-bits
  // signed integer) for cell_sizes and mesh_graph_partition. This may truncate
  // the high order bits of FleCSI's id_t (which is essentially a 64-bit
  // unsigned integer).
  // TODO: how many of these instance variables are truely needed?
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;
  std::vector<idx_t> cell_sizes;
  std::vector<flecsi::topology::mesh_graph_partition<idx_t>> cell_partitions;

  int comm_size;
  int rank;

  std::vector<idx_t> part;
};

TEST_F(mpi_parmetis_2way, comm_size_should_be_2) {
  ASSERT_THAT(comm_size, Eq(2));
}

TEST_F(mpi_parmetis_2way, vertices_are_in_partition_0_or_1) {
  ASSERT_THAT(part, Each(AnyOf(Eq(0), Eq(1))));
}

TEST_F(mpi_parmetis_2way, redistribute_cell_ids) {
  // number of cells on this rank.
  int num_cells = cell_sizes[rank];

  // we have to initialized cell_id_start to 0 since MPI_Excan does not update
  // it on Rank 0.
  int cell_id_start = 0;
  // Exclusive scan to find out the start cell entity id on this rank.
  // This assumes:
  //     1. cell entity id starts from 0
  //     2. cells are numbered consecutively in the initial partitioning DCSR
  //        as input to ParMetis.
  MPI_Exscan(&num_cells, &cell_id_start, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    ASSERT_THAT(cell_id_start, Eq(0));
  } else {
    ASSERT_THAT(cell_id_start, Eq(5));
  }

  // actually enumerate cell entity ids
  std::vector<idx_t> cell_ids(num_cells);
  std::iota(cell_ids.begin(), cell_ids.end(), cell_id_start);

  // sort cell ids according to their destination partitions
  thrust::stable_sort_by_key(part.begin(), part.end(), cell_ids.begin());
  if (rank == 0) {
    ASSERT_THAT(cell_ids, ElementsAreArray({0, 1, 2, 3, 4}));
  } else {
    ASSERT_THAT(cell_ids, ElementsAreArray({8, 9, 5, 6, 7}));
  }

  // calculate pairs of (destination, number of cells).
  std::vector<int> destinations(comm_size);
  std::vector<int> send_counts(comm_size);

  thrust::reduce_by_key(part.begin(), part.end(),
                        thrust::constant_iterator<int>(1),
                        destinations.begin(),
                        send_counts.begin());

  ASSERT_THAT(destinations, ElementsAre(0, 1));
  if (rank == 0) {
    ASSERT_THAT(send_counts, ElementsAre(3, 2));
  } else {
    ASSERT_THAT(send_counts, ElementsAre(2, 3));
  }

  // All to all communication to tell others how many cells are coming.
  std::vector<int> recv_counts(comm_size);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  if (rank == 0) {
    ASSERT_THAT(recv_counts, ElementsAreArray({3, 2}));
  } else {
    ASSERT_THAT(recv_counts, ElementsAreArray({2, 3}));
  }

  // calculate send and receive displacements from send and receive counts
  std::vector<int> send_disp(comm_size);
  thrust::exclusive_scan(send_counts.begin(), send_counts.end(),
                         send_disp.begin());
  if (rank == 0) {
    ASSERT_THAT(send_disp, ElementsAreArray({0, 3}));
  } else {
    ASSERT_THAT(send_disp, ElementsAreArray({0, 2}));
  }

  std::vector<int> recv_disp(comm_size);
  thrust::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                         recv_disp.begin());
  if (rank == 0) {
    ASSERT_THAT(recv_disp, ElementsAreArray({0, 3}));
  } else {
    ASSERT_THAT(recv_disp, ElementsAreArray({0, 2}));
  }

  // All to all communication to exchange cell ids.
  int total_recv = thrust::reduce(recv_counts.begin(), recv_counts.end());
  std::vector<int> my_cell_ids(total_recv);
  MPI_Alltoallv(cell_ids.data(),    send_counts.data(), send_disp.data(), MPI_INT,
                my_cell_ids.data(), recv_counts.data(), recv_disp.data(), MPI_INT,
                MPI_COMM_WORLD);
  if (rank == 0) {
    ASSERT_THAT(my_cell_ids, ElementsAreArray({0, 1, 2, 8, 9}));
  } else {
    ASSERT_THAT(my_cell_ids, ElementsAreArray({3, 4, 5, 6, 7}));
  }
}


TEST_F(mpi_parmetis_2way, redistribute_cell_2_cell_connectivity) {
  // give each element in adjncy list its partition number
  auto num_adjncy = cell_partitions[rank].index.size();
  std::vector<int> gather_list(num_adjncy);
  thrust::upper_bound(cell_partitions[rank].offset.begin()+1,
                      cell_partitions[rank].offset.end(),
                      thrust::counting_iterator<int>(0),
                      thrust::counting_iterator<int>(num_adjncy),
                      gather_list.begin());
  ASSERT_THAT(gather_list, ElementsAreArray({0, 0, 1, 1, 2, 2, 3, 3, 3, 4, 4}));

  std::vector<int> part_for_adjncy(num_adjncy);
  thrust::gather(gather_list.begin(), gather_list.end(),
                 part.begin(),
                 part_for_adjncy.begin());
  if (rank == 0)
    ASSERT_THAT(part_for_adjncy, ElementsAreArray({0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1}));
  else
    ASSERT_THAT(part_for_adjncy, ElementsAreArray({1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0}));


  // sort ajdncy according to their destination partition
  thrust::stable_sort_by_key(part_for_adjncy.begin(), part_for_adjncy.end(),
                             cell_partitions[rank].index.begin());
  if (rank == 0)
    ASSERT_THAT(cell_partitions[rank].index,
                ElementsAreArray({1, 9, 2, 0, 3, 1, 4, 8, 2, 5, 3}));
  else
    ASSERT_THAT(cell_partitions[rank].index,
                ElementsAreArray({9, 3, 7, 0, 8, 6, 4, 7, 5, 8, 6,}));

  // calculate how many elements of adjncy to send to each destination
  // in terms of pairs of (destination, number of cells).
  std::vector<int> destinations(comm_size);
  std::vector<int> send_counts(comm_size);
  thrust::reduce_by_key(part_for_adjncy.begin(), part_for_adjncy.end(),
                        thrust::constant_iterator<int>(1),
                        destinations.begin(),
                        send_counts.begin());

  ASSERT_THAT(destinations, ElementsAreArray({0, 1}));
  if (rank == 0)
    ASSERT_THAT(send_counts, ElementsAreArray({6, 5}));
  else
    ASSERT_THAT(send_counts, ElementsAreArray({5, 6}));

  // All to all communication to tell others how many elements of adjncy are
  // comming.
  std::vector<int> recv_counts(comm_size);
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  if (rank == 0)
    ASSERT_THAT(recv_counts, ElementsAreArray({6, 5}));
  else
    ASSERT_THAT(recv_counts, ElementsAreArray({5, 6}));

  // calculate send and recv displacements from send and recv counts
  std::vector<int> send_disp(comm_size);
  thrust::exclusive_scan(send_counts.begin(), send_counts.end(),
                         send_disp.begin());
  if (rank == 0)
    ASSERT_THAT(send_disp, ElementsAreArray({0, 6}));
  else
    ASSERT_THAT(send_disp, ElementsAreArray({0, 5}));
  // TODO: figure out if send_disp is always the same as recv_disp

  std::vector<int> recv_disp(comm_size);
  thrust::exclusive_scan(recv_counts.begin(), recv_counts.end(),
                         recv_disp.begin());
  if (rank == 0)
    ASSERT_THAT(recv_disp, ElementsAreArray({0, 6}));
  else
    ASSERT_THAT(recv_disp, ElementsAreArray({0, 5}));

  // All to all communication to exchange adjncy elements
  int total_recv = thrust::reduce(recv_counts.begin(), recv_counts.end());
  std::vector<int> my_adjncy(total_recv);
  MPI_Alltoallv(cell_partitions[rank].index.data(),
                send_counts.data(), send_disp.data(), MPI_INT,
                my_adjncy.data(),
                recv_counts.data(), recv_disp.data(), MPI_INT,
                MPI_COMM_WORLD);
  if (rank == 0)
    ASSERT_THAT(my_adjncy, ElementsAreArray({1, 9, 2, 0, 3, 1, 9, 3, 7, 0, 8}));
  else
    ASSERT_THAT(my_adjncy, ElementsAreArray({4, 8, 2, 5, 3, 6, 4, 7, 5, 8, 6,}));

  // number of cells on this rank.
  // TODO: how to make this more general? is this before redistribution or
  // after redistribution of cells?
  int num_cells = cell_sizes[rank];

  // calculate the "degrees" of each cell as the "length"
  std::vector<int> cell_degrees(num_cells);
  // Note: we assume xadj starts from 0 and shift the calculation by one element.
  //std::adjacent_difference(xadj[rank]+1, xadj[rank]+6, cell_degrees.begin());
  std::adjacent_difference(cell_partitions[rank].offset.begin()+1,
                           cell_partitions[rank].offset.end(),
                           cell_degrees.begin());
  ASSERT_THAT(cell_degrees, ElementsAreArray({2, 2, 2, 3, 2}));

  // sort cell degrees according to their destination partition
  thrust::stable_sort_by_key(part.begin(), part.end(), cell_degrees.begin());
  if (rank == 0)
    ASSERT_THAT(cell_degrees, ElementsAreArray({2, 2, 2, 3, 2}));
  else
    ASSERT_THAT(cell_degrees, ElementsAreArray({3, 2, 2, 2, 2}));

  // shuffle again
  thrust::reduce_by_key(part.begin(), part.end(),
                        thrust::constant_iterator<int>(1),
                        destinations.begin(),
                        send_counts.begin());
  if (rank == 0)
    ASSERT_THAT(send_counts, ElementsAreArray({3, 2}));
  else
    ASSERT_THAT(send_counts, ElementsAreArray({2, 3}));
  MPI_Alltoall(send_counts.data(), 1, MPI_INT,
               recv_counts.data(), 1, MPI_INT,
               MPI_COMM_WORLD);
  if (rank == 0)
    ASSERT_THAT(recv_counts, ElementsAreArray({3, 2}));
  else
    ASSERT_THAT(recv_counts, ElementsAreArray({2, 3}));

  thrust::exclusive_scan(send_counts.begin(), send_counts.end(), send_disp.begin());
  thrust::exclusive_scan(recv_counts.begin(), recv_counts.end(), recv_disp.begin());

  // exchange cell degrees.
  total_recv = thrust::reduce(recv_counts.begin(), recv_counts.end());
  std::vector<int> my_cell_degrees(total_recv);
  MPI_Alltoallv(cell_degrees.data(),    send_counts.data(), send_disp.data(), MPI_INT,
                my_cell_degrees.data(), recv_counts.data(), recv_disp.data(), MPI_INT,
                MPI_COMM_WORLD);
  if (rank == 0)
    ASSERT_THAT(my_cell_degrees, ElementsAreArray({2, 2, 2, 3, 2}));
  else
    ASSERT_THAT(my_cell_degrees, ElementsAreArray({3, 2, 2, 2, 2}));

  // reconstruct xadj from cell degrees
  std::vector<int> my_xadj(total_recv+1, 0);
  thrust::inclusive_scan(my_cell_degrees.begin(), my_cell_degrees.end(),
                         my_xadj.begin()+1);
  if (rank == 0)
    ASSERT_THAT(my_xadj, ElementsAreArray({0, 2, 4, 6, 9, 11}));
  else
    ASSERT_THAT(my_xadj, ElementsAreArray({0, 3, 5, 7, 9, 11}));

}