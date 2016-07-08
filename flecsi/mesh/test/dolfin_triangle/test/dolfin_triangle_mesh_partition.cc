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

#include <cinchtest.h>

#include "../dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace testing;

class A_Dolfin_Triangle_Partitioned_In_One : public ::testing::Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;

  virtual void SetUp() {
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);
  }

  // WARNING: we deliberately use int (which is either 32 or 64-bits
  // signed integer) for cell_sizes and mesh_graph_partition. This may truncate
  // the high order bits of FleCSI's id_t (which is essentially a 64-bit
  // unsigned integer).
  std::vector<int> cell_sizes = {10};
  std::vector<mesh_graph_partition<int>> cell_partitions;
};

TEST_F(A_Dolfin_Triangle_Partitioned_In_One,
       number_of_cell_partitions_should_be_1) {
  ASSERT_THAT(cell_partitions, SizeIs(1));
}

TEST_F(A_Dolfin_Triangle_Partitioned_In_One,
       dump_cell_2_cell_serial_CSR) {
  CINCH_CAPTURE() << "Serial CSR\n";
  CINCH_CAPTURE() << "xadj\t";
  for (auto x : cell_partitions[0].offset) {
    CINCH_CAPTURE() << x << " ";
  }
  CINCH_CAPTURE() << std::endl;

  CINCH_CAPTURE() << "adjncy\t";
  for (auto x : cell_partitions[0].index) {
    CINCH_CAPTURE() << x << " ";
  }
  CINCH_CAPTURE() << std::endl;
}

class A_Dolfin_Triangle_Partitioned_In_Two : public ::testing::Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;

  virtual void SetUp() {
    // convert and divide vertex to vertex and cell to cell connectivities into
    // two equal partitions in the form of Distributed CSR format as in ParMetis
    // manual.
    dolfin.compute_graph_partition(0, 0, vertex_sizes, vertex_partitions);
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);
  }

  // WARNING: we deliberately use int (which is either 32 or 64-bits
  // signed integer) for cell_sizes and mesh_graph_partition. This may truncate
  // the high order bits of FleCSI's id_t (which is essentially a 64-bit
  // unsigned integer).
  std::vector<int> vertex_sizes = {5, 5};
  std::vector<mesh_graph_partition<int>> vertex_partitions;
  std::vector<int> cell_sizes = {5, 5};
  std::vector<mesh_graph_partition<int>> cell_partitions;
};

TEST_F(A_Dolfin_Triangle_Partitioned_In_Two,
       number_of_vertex_partitions_should_be_2) {
  ASSERT_THAT(vertex_partitions, SizeIs(2));
}

TEST_F(A_Dolfin_Triangle_Partitioned_In_Two,
       number_of_vertex_in_each_partion_should_be_5) {
  ASSERT_THAT(vertex_partitions[0].offset, SizeIs(5));
  ASSERT_THAT(vertex_partitions[1].offset, SizeIs(5));
}

TEST_F(A_Dolfin_Triangle_Partitioned_In_Two,
       number_of_cell_partitions_should_be_2) {
  ASSERT_THAT(cell_partitions, SizeIs(2));
}

TEST_F(A_Dolfin_Triangle_Partitioned_In_Two,
       number_of_cell_in_each_partition_should_be_5) {
  ASSERT_THAT(cell_partitions[0].offset, SizeIs(5));
  ASSERT_THAT(cell_partitions[1].offset, SizeIs(5));
}

TEST_F(A_Dolfin_Triangle_Partitioned_In_Two, dump_cell_2_cell_DCSR) {
  CINCH_CAPTURE() << "Distributed CSR\n";
  for (auto cell_partition : cell_partitions) {
    CINCH_CAPTURE() << "xadj\t";
    for (auto x : cell_partition.offset) {
      CINCH_CAPTURE() << x << " ";
    }
    CINCH_CAPTURE() << std::endl;

    CINCH_CAPTURE() << "adjncy\t";
    for (auto x : cell_partition.index) {
      CINCH_CAPTURE() << x << " ";
    }
    CINCH_CAPTURE() << std::endl;

    CINCH_CAPTURE() << "vtxdist\t";
    for (auto x : cell_partition.partition) {
      CINCH_CAPTURE() << x << " ";
    }
    CINCH_CAPTURE() << std::endl;
  }
}