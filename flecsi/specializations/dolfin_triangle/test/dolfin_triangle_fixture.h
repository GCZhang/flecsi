/*~--------------------------------------------------------------------------~*
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
 *~--------------------------------------------------------------------------~*/

#ifndef FLECSI_DOLFIN_TRIANGLE_FIXTURE_H
#define FLECSI_DOLFIN_TRIANGLE_FIXTURE_H

namespace flecsi {

class A_Dolfin_Triangle : public ::testing::Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;
};

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
  std::vector<topology::mesh_graph_partition<int>> cell_partitions;
};

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

  std::vector<size_t> vertex_sizes = {5, 5};
  std::vector<flecsi::topology::mesh_graph_partition<size_t>>
    vertex_partitions;
  std::vector<size_t> cell_sizes = {5, 5};
  std::vector<flecsi::topology::mesh_graph_partition<size_t>> cell_partitions;
};
}
#endif //FLECSI_DOLFIN_TRIANGLE_FIXTURE_H
