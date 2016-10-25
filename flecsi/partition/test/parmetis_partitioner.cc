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

#include "flecsi/partition/parmetis_partitioner.h"
#include "flecsi/specializations/dolfin_triangle/dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace testing;

class parmetis_partitioner_fixture : public Test {
protected:
  virtual  void SetUp() {
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto cells_per_rank = dolfin.num_cells()/comm_size;
    for (auto i = 0; i < comm_size - 1 ; i++) {
      cell_sizes.push_back(cells_per_rank);
    }
    cell_sizes.push_back(dolfin.num_cells() - cells_per_rank*(comm_size-1));

    // create Distributed CSR from dolfin's cell to cell connectivity.
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);

    auto partitioner = parmetis_partitioner();
    partitioner.partition(cell_partitions, index_partition);
  }

  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;
  std::vector<idx_t> cell_sizes;

  using mesh_graph_partitions = std::vector<mesh_graph_partition<idx_t>>;
  mesh_graph_partitions cell_partitions;

  using index_partition_t = index_partition__<idx_t>;
  index_partition_t index_partition;

  int comm_size;
  int rank;
};

TEST_F(parmetis_partitioner_fixture, comm_size_should_be_2) {
  ASSERT_THAT(comm_size, Eq(2));
}

TEST_F(parmetis_partitioner_fixture, dump) {
  if (rank == 0) {
    ASSERT_THAT(index_partition.exclusive, ElementsAreArray({0, 1, 2, 8, 9}));
  } else {
    ASSERT_THAT(index_partition.exclusive, ElementsAreArray({3, 4, 5, 6, 7}));
  }
}