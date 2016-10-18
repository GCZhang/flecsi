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

#include <metis.h>

#include "../dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace testing;

class metis_partition_graph_2way : public Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;

  virtual void  SetUp() override {
    // get cell to cell connectivity, this is the "Graph" in Metis' lingo.
    idx_t num_cells = dolfin.num_cells();
    cell_sizes.push_back(num_cells);
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);

    idx_t ncon = 1;
    idx_t nparts = 2;
    idx_t objval;

    part.resize(num_cells);

    auto ret = METIS_PartGraphKway(
      &num_cells,
      &ncon,
      cell_partitions[0].offset.data(),
      cell_partitions[0].index.data(),
      nullptr, /* vwgt */
      nullptr, /* vsize */
      nullptr, /* adjwt */
      &nparts,
      nullptr, /* tpwgts */
      nullptr, /* ubvec */
      nullptr, /* options */
      &objval,
      part.data()
    );

    ASSERT_EQ(ret, METIS_OK);
  }

  // WARNING: we deliberately use Metis' idx_t (which is either 32 or 64-bits
  // signed integer) for cell_sizes and mesh_graph_partition. This may truncate
  // the high order bits of FleCSI's id_t (which is essentially a 64-bit
  // unsigned integer).
  std::vector<idx_t> cell_sizes;
  std::vector<flecsi::topology::mesh_graph_partition<idx_t>> cell_partitions;
  std::vector<idx_t> part;
};

TEST_F(metis_partition_graph_2way, there_are_two_partitions) {
   ASSERT_THAT(part, Each(AnyOf(Eq(0), Eq(1))));
}

TEST_F(metis_partition_graph_2way, dump) {
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

   CINCH_CAPTURE() << "partition ";
   for (auto i: part) {
     CINCH_CAPTURE() << i << " ";
   }
   CINCH_CAPTURE() << std::endl;
   ASSERT_TRUE(CINCH_EQUAL_BLESSED("metis_partition_graph_kway.blessed"));
 }