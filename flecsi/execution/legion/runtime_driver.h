/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_execution_legion_runtime_driver_h
#define flecsi_execution_legion_runtime_driver_h

#include <legion.h>

///
// \file legion/runtime_driver.h
// \authors bergen
// \date Initial file creation: Jul 26, 2016
///

namespace flecsi {
namespace execution {

void legion_runtime_driver(const LegionRuntime::HighLevel::Task * task,
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions,
  LegionRuntime::HighLevel::Context ctx,
  LegionRuntime::HighLevel::HighLevelRuntime * runtime);

} // namespace execution 
} // namespace flecsi

#endif // flecsi_execution_legion_runtime_driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
