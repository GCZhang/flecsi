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

class metis_partition_mesh_nodal : public Test {
protected:
  dolfin_triangle_mesh_t<dolfin_triangle_types_t> dolfin;

  virtual void SetUp() override {
    auto& conn = dolfin.get_connectivity(0, 2, 0);

    // we use vector.resize() there because we want to also change vector.end()
    // when we reserve memory. vector.reserve() DOES NOT change vector.end().
    idx_t num_cells = dolfin.num_cells();
    epart.resize(num_cells);
    idx_t num_vertices = dolfin.num_vertices();
    npart.resize(num_vertices);

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
    
    idx_t num_parts = 2;
    idx_t objval;

    auto ret = METIS_PartMeshNodal(
      &num_cells,
      &num_vertices,
      eptr.data(),
      eind.data(),
      nullptr,
      nullptr,
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

TEST_F(metis_partition_mesh_nodal,
       there_are_more_than_0_cells_in_partition_0) {
  auto count0 = std::count_if(epart.begin(), epart.end(),
                              [](auto part) {return part == 0;});
  ASSERT_THAT(count0, Ne(0));
}

TEST_F(metis_partition_mesh_nodal,
       there_are_more_than_0_cells_in_partition_1) {
  auto count1 = std::count_if(epart.begin(), epart.end(),
                              [](auto part) {return part == 1;});
  ASSERT_THAT(count1, Ne(0));
}

TEST_F(metis_partition_mesh_nodal,
       there_are_the_same_number_of_vertices_in_each_partition) {
  auto count0 = std::count_if(npart.begin(), npart.end(),
                              [](auto part) {return part == 0;});
  auto count1 = std::count_if(npart.begin(), npart.end(),
                              [](auto part) {return part == 1;});
  ASSERT_EQ(count0, count1);
}
