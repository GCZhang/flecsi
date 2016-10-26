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

#
#include "flecsi/partition/parmetis_partitioner.h"
#include "flecsi/specializations/dolfin_triangle/dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace flecsi::topology;
using namespace testing;

class parmetis_partitioner_ugrid2d : public Test {
protected:
  static constexpr int N = 4;

  virtual  void SetUp() {
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // FIXME: assuming N*N is evenly divisible by comm_size;
    auto cells_per_rank = N*N/comm_size;
    auto cell_id_start = cells_per_rank * rank;

    std::vector<idx_t> xadj;
    std::vector<idx_t> adjncy;
    std::vector<idx_t> vtxdist;
    xadj.push_back(0);
    for (auto id = cell_id_start, offset = 0; id < cell_id_start + cells_per_rank; id++) {
      const int i = id % N;
      const int j = id / N;
      // W, E, S, N
      if (i != 0) {
        // not a left boundary cell, add its left neighbor
        adjncy.push_back(id-1);
        offset++;
      }
      if (i != N-1) {
        // not right boundary cell, add its right neighbor
        adjncy.push_back(id+1);
        offset++;
      }
      if (j != 0) {
        // not a bottom boundary cell, add its bottom neighbor
        adjncy.push_back(id - N);
        offset++;
      }
      if (j != N-1) {
        // not a top boundary cell, add its top neighbor
        adjncy.push_back(id + N);
        offset++;
      }

      xadj.push_back(offset);
    }

    for (auto id = 0; id <= N*N; id += cells_per_rank) {
      vtxdist.push_back(id);
    }

    cell_partition.offset = std::move(xadj);
    cell_partition.index = std::move(adjncy);
    cell_partition.partition = std::move(vtxdist);

    auto partitioner = parmetis_partitioner();
    partitioner.partition(cell_partition, index_partition);
  }

  //using mesh_graph_partitions = std::vector<mesh_graph_partition<idx_t>>;
  //mesh_graph_partitions cell_partitions;
  using  mesh_partition = mesh_graph_partition<idx_t>;
  mesh_partition cell_partition;

  using index_partition_t = index_partition__<idx_t>;
  index_partition_t index_partition;

  int comm_size;
  int rank;
};

TEST_F(parmetis_partitioner_ugrid2d, dump) {
  std::cout << "cell id: ";
  for (auto i : index_partition.exclusive) {
    std::cout << i << " ";
  }
  std::cout << std::endl;
}