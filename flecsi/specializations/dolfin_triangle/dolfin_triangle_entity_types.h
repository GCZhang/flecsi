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

#ifndef FLECSI_DOLFIN_TRIANGLE_ENTITY_TYPES_H
#define FLECSI_DOLFIN_TRIANGLE_ENTITY_TYPES_H

#include "flecsi/topology/mesh_types.h"

namespace flecsi
{
class dolfin_vertex_t : public flecsi::topology::mesh_entity_t<0, 1> {
public:
  dolfin_vertex_t() = default;
  dolfin_vertex_t(flecsi::topology::mesh_topology_base_t &) {}
};

class dolfin_edge_t : public flecsi::topology::mesh_entity_t<1, 1> {
public:
  dolfin_edge_t() = default;
  dolfin_edge_t(flecsi::topology::mesh_topology_base_t &) {}
};

class dolfin_cell_t : public flecsi::topology::mesh_entity_t<2, 1> {
public:
  dolfin_cell_t() = default;
  dolfin_cell_t(flecsi::topology::mesh_topology_base_t &){}
  
  std::vector<size_t>
<<<<<<< HEAD:flecsi/mesh/test/dolfin_triangle/dolfin_triangle_entity_types.h
  create_entities(flecsi::id_t cell_id, size_t dim, domain_connectivity<2> & c,
                  flecsi::id_t * e){
=======
  create_entities(
    flecsi::id_t cell_id,
    size_t dim,
    flecsi::topology::domain_connectivity<2> & c,
    flecsi::id_t * e
  )
  {
>>>>>>> remotes/origin/execution:flecsi/specializations/dolfin_triangle/dolfin_triangle_entity_types.h
    flecsi::id_t* v = c.get_entities(cell_id, 0);
    e[0] = v[0];
    e[1] = v[1];

    e[2] = v[1];
    e[3] = v[2];

    e[4] = v[2];
    e[5] = v[0];
    return {2, 2, 2};
  }
};
}
#endif //FLECSI_DOLFIN_TRIANGLE_ENTITY_TYPES_H
