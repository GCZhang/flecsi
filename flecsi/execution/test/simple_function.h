/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_default_driver_h
#define flecsi_default_driver_h

#include <iostream>

#include "flecsi/utils/common.h"
#include "flecsi/execution/context.h"
#include "flecsi/execution/execution.h"
#include "flecsi/data/data.h"

/*!
 * \file default_driver.h
 * \authors bergen
 * \date Initial file creation: Jul 24, 2016
 */

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
// Function registration.
//----------------------------------------------------------------------------//

double test_function(double r, double e) {
  std::cout << "Executing test_function" << std::endl;
  std::cout << "(r,e): (" << r << "," << e << ")" << std::endl;
  return r*e;
} // function1

register_function(test_function);

//----------------------------------------------------------------------------//
// Driver.
//----------------------------------------------------------------------------//

void driver(int argc, char ** argv) {
  auto handle = function_handle(test_function);

  double result = execute_function(handle, 2.0, 10.0);

  std::cout << "test_function returned: " << result << std::endl;
} // driver

} // namespace execution
} // namespace flecsi

#endif // flecsi_default_driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
