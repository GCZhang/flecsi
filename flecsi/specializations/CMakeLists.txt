#~----------------------------------------------------------------------------~#
#  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
# /@@/////  /@@          @@////@@ @@////// /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@
# /@@       /@@/@@//// //@@    @@       /@@/@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
# //       ///  //////   //////  ////////  // 
# 
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

set(specializations_HEADERS
  burton/burton_entity_types.h
  burton/burton.h
  burton/burton_io_exodus.h
  burton/burton_mesh.h
  burton/burton_mesh_traits.h
  burton/burton_types.h
  PARENT_SCOPE)

if(ENABLE_IO)
  list(APPEND specializations_HEADERS burton_io_exodus.h)
endif(ENABLE_IO)

set(specializations_SOURCES
  burton/burton_entity_types.cc
  PARENT_SCOPE)

#------------------------------------------------------------------------------#
# Burton specialization.
#------------------------------------------------------------------------------#

cinch_add_unit(burton
  SOURCES burton/test/burton.cc
  INPUTS burton/test/burton.blessed
  LIBRARIES flecsi)

cinch_add_unit(copy
  SOURCES burton/test/copy.cc
  LIBRARIES flecsi)

cinch_add_unit(connectivity
  SOURCES burton/test/connectivity.cc
  LIBRARIES flecsi)

if(LAPACKE_FOUND AND ENABLE_IO)
  cinch_add_unit(burton_gradients
    SOURCES burton/test/burton_gradients.cc
    LIBRARIES flecsi ${LAPACKE_LIBRARIES} ${IO_LIBRARIES})
endif(LAPACKE_FOUND AND ENABLE_IO)

if(ENABLE_IO)
  cinch_add_unit(burton_io
    SOURCES burton/test/burton_io.cc
    LIBRARIES ${IO_LIBRARIES} flecsi)
endif(ENABLE_IO)

#------------------------------------------------------------------------------#
# Basic-tree specialization.
#------------------------------------------------------------------------------#

cinch_add_unit(basic-tree
  SOURCES basic-tree/test/basic-tree.cc
  LIBRARIES flecsi
)

#------------------------------------------------------------------------------#
#  Dolfin specialization.
#------------------------------------------------------------------------------#

cinch_add_unit(dolfin_triangle
  SOURCES dolfin_triangle/test/dolfin_triangle.cc
  INPUTS dolfin_triangle/test/dolfin_triangle.blessed)

cinch_add_unit(dolfin_triangle_mesh_partition
  SOURCES dolfin_triangle/test/dolfin_triangle_mesh_partition.cc)

cinch_add_unit(dolfin_triangle_connectivity
  SOURCES dolfin_triangle/test/dolfin_triangle_connectivity.cc)

if (ENABLE_PARTITION AND METIS_FOUND)
  cinch_add_unit(dolfin_triangle_metis
    SOURCES dolfin_triangle/test/metis_partition_mesh_dual.cc
    LIBRARIES ${PARTITION_LIBRARIES})
  cinch_add_unit(metis_partition_mesh_nodal
    SOURCES dolfin_triangle/test/metis_partition_mesh_nodal.cc
    LIBRARIES ${PARTITION_LIBRARIES})
endif()

#------------------------------------------------------------------------------#
#  Sagittarius specialization.
#------------------------------------------------------------------------------#

cinch_add_unit(sagittarius
  SOURCES sagittarius/test/sagittarius.cc
  INPUTS sagittarius/test/sagittarius.blessed)

cinch_add_unit(sagittarius_mesh_partition
  SOURCES sagittarius/test/sagittarius_mesh_partition.cc)

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
