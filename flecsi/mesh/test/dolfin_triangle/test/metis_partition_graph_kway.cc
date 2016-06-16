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
    // TODO: this gives c2c via shared vertices, we want c2c via shared edges.
    //connectivity_t conn = dolfin.get_connectivity(0, 2, 2);

    idx_t num_cells = dolfin.num_cells();

    std::vector<size_t> cell_sizes;
    cell_sizes.push_back(num_cells);
    std::vector<mesh_graph_partition<size_t>> cell_partitions;
    dolfin.compute_graph_partition(0, 2, cell_sizes, cell_partitions);
    std::cout << "num partitions: " << cell_partitions.size() << std::endl;

    auto &from_index = cell_partitions[0].offset;
    std::vector<idx_t> xadj(from_index.begin(), from_index.end());
    for (auto x : from_index) {
      std::cout << x << " ";
    }
    std::cout << std::endl;

    auto &to_index = cell_partitions[0].index;
    std::vector<idx_t> adjncy(to_index.begin(), to_index.end());
    for (auto x : adjncy) {
      std::cout << x << " ";
    }
    std::cout << std::endl;
  }

};

TEST_F(metis_partition_graph_2way, there_are_two_partitions) {

}