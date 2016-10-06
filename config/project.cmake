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

project(flecsi)

#------------------------------------------------------------------------------#
# Set application directory
#------------------------------------------------------------------------------#

cinch_add_application_directory("examples")
cinch_add_application_directory("examples/agile")
cinch_add_application_directory("bin")

#------------------------------------------------------------------------------#
# Add library targets
#------------------------------------------------------------------------------#

cinch_add_library_target(flecsi flecsi)

#------------------------------------------------------------------------------#
# Set header suffix regular expression
#------------------------------------------------------------------------------#

set(CINCH_HEADER_SUFFIXES "\\.h")

#------------------------------------------------------------------------------#
# Add options for runtime selection
#------------------------------------------------------------------------------#

set(FLECSI_RUNTIME_MODEL "serial" CACHE STRING
  "Select the runtime model [legion,mpi,mpilegion,serial]")

#------------------------------------------------------------------------------#
# Add option for setting id bits
#------------------------------------------------------------------------------#

set(FLECSI_ID_PBITS "20" CACHE STRING
  "Select the number of bits to use for partition ids. There will be 62-FLECSI_ID_PBITS-FLECSI_ID_FBITS available for entity ids")

set(FLECSI_ID_FBITS "4" CACHE STRING
  "Select the number of bits to use for id flags. There will be 62-FLECSI_ID_PBITS-FLECSI_ID_FBITS available for entity ids")

#------------------------------------------------------------------------------#
# Add option for counter size
#------------------------------------------------------------------------------#

set(FLECSI_COUNTER_TYPE "int32_t" CACHE STRING
  "Select the type that will be used for loop and iterator values")

#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
