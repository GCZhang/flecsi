/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_runtime_execution_policy_h
#define flecsi_runtime_execution_policy_h

///
// \file flecsi_runtime_policy.h
// \authors bergen
// \date Initial file creation: Aug 01, 2016
///

///
// This section works with the build system to select the correct runtime
// implemenation for the task model. If you add to the possible runtimes,
// remember to edit config/packages.cmake to include a definition using
// the same convention, e.g., -DFLECSI_RUNTIME_MODEL_new_runtime.
///

// Serial Policy
#if defined(FLECSI_RUNTIME_MODEL_serial)

  #include "flecsi/execution/serial/execution_policy.h"
  #define flecsi_execution_policy_t \
    flecsi::execution::serial_execution_policy_t

// Legion Policy
#elif defined(FLECSI_RUNTIME_MODEL_legion)

  #include "flecsi/execution/legion/execution_policy.h"
  #define flecsi_execution_policy_t \
    flecsi::execution::legion_execution_policy_t

// MPI+Legion Policy
#elif defined(FLECSI_RUNTIME_MODEL_mpilegion)

  #include "flecsi/execution/mpilegion/execution_policy.h"
  #define flecsi_execution_policy_t \
    flecsi::execution::mpilegion_execution_policy_t

// MPI+Legion Policy
#elif defined(FLECSI_RUNTIME_MODEL_mpi)

  #include "flecsi/execution/mpi/execution_policy.h"
  #define flecsi_execution_policy_t \
    flecsi::execution::mpi_execution_policy_t

#endif // FLECSI_RUNTIME_MODEL

#endif // flecsi_runtime_execution_policy_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
