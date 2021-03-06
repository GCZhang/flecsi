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

set(execution_HEADERS
	common/function_handle.h
	common/launch.h
	common/processor.h
	common/task_hash.h
	context.h
	default_driver.h
	execution.h
	function.h
	kernel.h
	legion/context_policy.h
	legion/future.h
	legion/execution_policy.h
	legion/runtime_driver.h
	legion/task_args.h
	legion/task_wrapper.h
	mpilegion/context_policy.h
	mpilegion/execution_policy.h
	mpilegion/legion_handshake.h
	mpilegion/mapper.h
	mpilegion/mpi_legion_interop.h
	mpilegion/runtime_driver.h
	mpilegion/task_ids.h
	mpilegion/task_wrapper.h
	serial/context_policy.h
	serial/execution_policy.h
	serial/runtime_driver.h
	task.h
  PARENT_SCOPE
)


set(execution_SOURCES
  PARENT_SCOPE
)

if(FLECSI_RUNTIME_MODEL STREQUAL "legion")

 cinch_add_unit(task
    SOURCES
      test/task.cc
      test/task_driver.h
      legion/runtime_driver.cc
    DEFINES FLECSI_DRIVER=flecsi/execution/test/task_driver.h
    POLICY LEGION
    LIBRARIES flecsi)

elseif (FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")

#  cinch_add_unit(task
#    SOURCES
#      test/task.cc
#      test/task_driver.h
#      mpilegion/runtime_driver.cc
#      ../partition/init_partitions_task.cc
#    DEFINES FLECSI_DRIVER=flecsi/execution/test/task_driver.h
#    POLICY MPILEGION
#    LIBRARIES flecsi
#  )

#  cinch_add_unit(sprint
#    SOURCES
#      test/task.cc
#      test/sprint.h
#      mpilegion/runtime_driver.cc
#      ../partition/init_partitions_task.cc
#    DEFINES SPRINT_TEST
#    DEFINES FLECSI_DRIVER=flecsi/execution/test/sprint.h
#    POLICY MPILEGION
#    LIBRARIES flecsi
#  )

else() # serial

  cinch_add_unit(task
    SOURCES
      test/task.cc
      test/task_driver.h
      serial/runtime_driver.cc
    DEFINES FLECSI_DRIVER=flecsi/execution/test/task_driver.h
    LIBRARIES flecsi
  )

  cinch_add_unit(simple_function
    SOURCES
      test/function.cc
      test/simple_function.h
      serial/runtime_driver.cc
    DEFINES FLECSI_DRIVER=flecsi/execution/test/simple_function.h
    LIBRARIES flecsi
  )
endif()


if(FLECSI_RUNTIME_MODEL STREQUAL "mpilegion")
  
#  cinch_add_unit(test_legion_handshake
#    SOURCES test/test_legion_handshake.cc
#    POLICY MPILEGION
#  )
#
#  cinch_add_unit(test_mpilegion_interop
#    SOURCES test/test_mpi_legion_interop.cc
#    POLICY MPILEGION
#  )
 
endif()

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
