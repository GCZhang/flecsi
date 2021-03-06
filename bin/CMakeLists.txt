#-----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#-----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Add source files for example application.
#------------------------------------------------------------------------------#

set(example_src example_app.cc)
set(exec_path ${CMAKE_SOURCE_DIR}/flecsi/execution)

if(FLECSI_RUNTIME_MODEL STREQUAL "serial")
    list(APPEND example_src ${exec_path}/serial/runtime_driver.cc)
elseif(FLECSI_RUNTIME_MODEL STREQUAL "legion")
    list(APPEND example_src ${exec_path}/legion/runtime_driver.cc)
elseif(FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")
    list(APPEND example_src ${exec_path}/mpilegion/runtime_driver.cc
                            ${exec_path}/../partition/init_partitions_task.cc)
else()
    message(FATAL_ERROR "This runtime is not yet supported")
endif()

#------------------------------------------------------------------------------#
# Add a rule to build the executable.
#------------------------------------------------------------------------------#

message(STATUS "example source: ${example_src}")
add_executable(example_app ${example_src})

#------------------------------------------------------------------------------#
# Set compile define for driver.
#------------------------------------------------------------------------------#

target_compile_definitions(example_app PRIVATE
    -DFLECSI_DRIVER=${CMAKE_SOURCE_DIR}/bin/example_driver.h)

#------------------------------------------------------------------------------#
# Add link dependencies.
#------------------------------------------------------------------------------#

target_link_libraries(example_app flecsi ${FLECSI_RUNTIME_LIBRARIES})

#~---------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=4 shiftwidth=4 expandtab :
#~---------------------------------------------------------------------------~-#
