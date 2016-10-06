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

#include "flecsi/specializations/dolfin_triangle/dolfin_triangle_mesh.h"

using namespace flecsi;
using namespace testing;

class metis_partition_mesh_dual : public Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;

  virtual void SetUp() override {
<<<<<<< HEAD:flecsi/mesh/test/dolfin_triangle/test/metis_partition_mesh_dual.cc
    // get cell to vertices connectivity
    connectivity_t conn = dolfin.get_connectivity(0, 2, 0);
=======
    flecsi::topology::connectivity_t& conn = dolfin.get_connectivity(0, 2, 0);
>>>>>>> remotes/origin/execution:flecsi/specializations/dolfin_triangle/test/metis_partition_mesh_dual.cc

    // we use vector.resize() there because we want to also change vector.end()
    // when we reserve memory. vector.reserve() DOES NOT change vector.end().
    idx_t num_cells = dolfin.num_cells();
    epart.resize(static_cast<size_t>(num_cells));
    idx_t num_vertices = dolfin.num_vertices();
    npart.resize(static_cast<size_t>(num_vertices));

    // Because the data type of indices (id_t a.k.a size_t) in connectivity_t is
    // not the same as idx_t (a.k.a. int32 or int64) Metis expects, we have
    // to manually copy the elements to convert types.
    auto from_index = conn.get_from_index_vec();
    std::vector<idx_t> eptr(from_index.begin(), from_index.end());
    // The same goes for to_index. However, the elements of to_index are of
    // type id_t which has an implicit type conversion operator to size_t.
    auto to_index = conn.get_entities();
    std::vector<idx_t> eind;

    for(auto to_id : to_index){
      eind.push_back(to_id.entity());
    }

    // METIS_PartMeshDual will partition the mesh base on cell to cell
    // connectivities. It recovers this information from the cell to vertices
    // connectivies using number of common vertices two cells have.
    // num_common = 2 means two cells are connected when they share a common
    // edge.
    idx_t num_common = 2;
    idx_t num_parts = 2;
    idx_t objval;

    auto ret = METIS_PartMeshDual(
      &num_cells,
      &num_vertices,
      eptr.data(),
      eind.data(),
      nullptr,
      nullptr,
      &num_common,
      &num_parts,
      nullptr,
      nullptr,
      &objval,
      epart.data(),
      npart.data()
    );
    ASSERT_EQ(ret, METIS_OK);
  }

  std::vector<idx_t> epart;
  std::vector<idx_t> npart;
};

TEST_F(metis_partition_mesh_dual, cells_are_partitioned_into_partition_0_or_1) {
  ASSERT_THAT(epart, Each(Ge(0)));
  ASSERT_THAT(epart, Each(Le(1)));
}

TEST_F(metis_partition_mesh_dual, there_are_5_cells_in_each_partition) {
  auto count0 = std::count_if(epart.begin(), epart.end(),
                             [](auto part) {return part == 0;});
  ASSERT_THAT(count0, Eq(5));

  auto count1 = std::count_if(epart.begin(), epart.end(),
                              [](auto part) {return part == 1;});
  ASSERT_THAT(count1, Eq(5));
}

TEST_F(metis_partition_mesh_dual, vertices_are_partitioned_into_partition_0_or_1) {
  ASSERT_THAT(npart, Each(Ge(0)));
  ASSERT_THAT(epart, Each(Le(1)));
}

#if IDXTYPEWIDTH == 64
TEST_F(metis_partition_mesh_dual,
       there_are_6_and_4_vertices_in_each_partition_for_64bits_ids) {
  auto count0 = std::count_if(npart.begin(), npart.end(), [](auto part) {return part == 0;});

  ASSERT_THAT(count0, Eq(6));

  auto count1 = std::count_if(npart.begin(), npart.end(), [](auto part) {return part == 1;});
  ASSERT_THAT(count1, Eq(4));
}
#else
TEST_F(metis_partition_mesh_dual,
       there_are_5_vertices_in_each_partition_for_32bits_ids) {
  auto count0 = std::count_if(npart.begin(), npart.end(),
                              [](auto part) {return part == 0;});

  ASSERT_THAT(count0, Eq(5));

  auto count1 = std::count_if(npart.begin(), npart.end(),
                              [](auto part) {return part == 1;});
  ASSERT_THAT(count1, Eq(5));
}
#endif
