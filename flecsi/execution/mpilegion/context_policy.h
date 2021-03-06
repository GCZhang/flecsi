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

#ifndef flecsi_execution_mpilegion_context_policy_h
#define flecsi_execution_mpilegion_context_policy_h

/*!
 * \file mpilegion/context_policy.h
 * \authors bergen
 * \date Initial file creation: Jul 14, 2016
 */

#include <memory>
#include <functional>
#include <unordered_map>
#include <stack>

#include <legion.h>

#include "flecsi/utils/common.h"
#include "flecsi/utils/const_string.h"
#include "flecsi/utils/tuple_wrapper.h"
#include "flecsi/execution/mpilegion/runtime_driver.h"
#include "flecsi/execution/common/task_hash.h"
#include "flecsi/execution/mpilegion/legion_handshake.h"
#include "flecsi/execution/mpilegion/mpi_legion_interop.h"
#include "flecsi/partition/init_partitions_task.h"

namespace flecsi {
namespace execution {

///
// \struct legion_runtime_runtime_state_t legion/context_policy.h
// \brief legion_runtime_state_t provides storage for Legion runtime
//        information that can be reinitialized as needed to store const
//        data types and references as required by the Legion runtime.
///
struct legion_runtime_state_t {

  legion_runtime_state_t(
    LegionRuntime::HighLevel::Context & context_,
    LegionRuntime::HighLevel::HighLevelRuntime * runtime_,
    const LegionRuntime::HighLevel::Task * task_,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions_
  )
  :
    context(context_),
    runtime(runtime_),
    task(task_),
    regions(regions_)
  {}
    
  LegionRuntime::HighLevel::Context & context;
  LegionRuntime::HighLevel::HighLevelRuntime * runtime;
  const LegionRuntime::HighLevel::Task * task;
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions;

}; // struct legion_runtime_state_t

// Use thread local storage for legion state information. The state_
// is set for each legion task invocation using the task name hash
// as a key. This seems like it should be safe, since multiple concurrent
// invocations of the same task can only occur on seperate threads.
static thread_local std::unordered_map<size_t,
  std::stack<std::shared_ptr<legion_runtime_state_t>>> state_;

///
// \class mpilegion_context_policy_t mpilegion/context_policy.h
// \brief mpilegion_context_policy_t provides...
///
struct mpilegion_context_policy_t
{

  using lr_context_t = LegionRuntime::HighLevel::Context;
  using lr_runtime_t = LegionRuntime::HighLevel::HighLevelRuntime;
  using lr_task_t = LegionRuntime::HighLevel::Task;
  using lr_regions_t =
    std::vector<LegionRuntime::HighLevel::PhysicalRegion>;

  const static LegionRuntime::HighLevel::Processor::Kind lr_loc =
    LegionRuntime::HighLevel::Processor::LOC_PROC;

  const size_t TOP_LEVEL_TASK_ID = 0;

  ext_legion_handshake_t & handshake_ = ext_legion_handshake_t::instance();
  mpi_legion_interop_t interop_helper_;

  /*--------------------------------------------------------------------------*
   * Initialization.
   *--------------------------------------------------------------------------*/

  int
  initialize(
    int argc,
    char ** argv
  )
  {
    handshake_.initialize(ext_legion_handshake_t::IN_EXT, 1,1);

    // Register top-level task
    lr_runtime_t::set_top_level_task_id(TOP_LEVEL_TASK_ID);
    lr_runtime_t::register_legion_task<mpilegion_runtime_driver>(
      TOP_LEVEL_TASK_ID, lr_loc, true, false);

    // FIXME
		// This is Galen's hack to get partitioning working for the sprint
    lr_runtime_t::register_legion_task<flecsi::dmp::parts,
      flecsi::dmp::init_partitions>(
      task_ids_t::instance().init_cell_partitions_task_id,lr_loc, true, false);

   // FIXME
    // This is Galen's hack to get partitioning working for the sprint
    lr_runtime_t::register_legion_task<flecsi::dmp::fill_cells_global_task>(
      task_ids_t::instance().init_cells_global_task_id,lr_loc, true, false);
 
   // FIXME
    // This is Galen's hack to get partitioning working for the sprint
    lr_runtime_t::register_legion_task<flecsi::dmp::find_ghost_task>(
      task_ids_t::instance().find_ghost_task_id,lr_loc, true, false);

    // register connect_to_mpi_task from mpi_legion_interop_t class
    lr_runtime_t::register_legion_task<connect_to_mpi_task>(
      task_ids_t::instance().connect_mpi_task_id, lr_loc, false, true,
      AUTO_GENERATE_ID,
      LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
      "connect_to_mpi_task");

    // register handoff_to_mpi_task from mpi_legion_interop_t class
    lr_runtime_t::register_legion_task<handoff_to_mpi_task>(
      task_ids_t::instance().handoff_to_mpi_task_id, lr_loc,
      false, true, AUTO_GENERATE_ID,
      LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
      "handoff_to_mpi_task");

    // register wait_on_mpi_task from mpi_legion_interop_t class
    lr_runtime_t::register_legion_task<wait_on_mpi_task>(
      task_ids_t::instance().wait_on_mpi_task_id, lr_loc,
      false, true, AUTO_GENERATE_ID,
      LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
      "wait_on_mpi_task");

		// register unset_call_mpi_task from mpi_legion_interop_t class
    lr_runtime_t::register_legion_task<unset_call_mpi_task>(
			task_ids_t::instance().unset_call_mpi_id, lr_loc,
			false, true, AUTO_GENERATE_ID,
			LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
			"unset_call_mpi_task");

    // Register user tasks
    for(auto f: task_registry_) {
      // funky logic: task_registry_ is a map of std::pair
      // f.first is the uintptr_t that holds the user function address
      // f.second is the pair of unique task id and the registration function
      f.second.second(f.second.first);
    } // for
  
    // FIXME: Documentation!
    interop_helper_.initialize();  

    // Start the runtime
    lr_runtime_t::start(argc, argv,true);

    interop_helper_.legion_configure();

    interop_helper_.handoff_to_legion();

    interop_helper_.wait_on_legion();

    //while loop to do some mpi tasks
    while(ext_legion_handshake_t::instance().call_mpi_) {
      #ifdef LEGIONDEBUG
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank);
        std::cout << "inside while loop N " << " rank = "<< rank <<
          "rank_from_handshake = " <<
          ext_legion_handshake_t::instance().rank_ << std::endl;
      #endif

      ext_legion_handshake_t::instance().shared_func_();
      interop_helper_.handoff_to_legion();
      interop_helper_.wait_on_legion();
    } // while

    interop_helper_.wait_on_legion();
   
    return 0;
  } // initialize

  void push_state(
    size_t key,
    LegionRuntime::HighLevel::Context & context,
    LegionRuntime::HighLevel::HighLevelRuntime * runtime,
    const LegionRuntime::HighLevel::Task * task,
    const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions
  )
  {
    #ifndef NDEBUG
      std::cout << "pushing state for " << key << std::endl;
    #endif

    state_[key].push(std::shared_ptr<legion_runtime_state_t>
      (new legion_runtime_state_t(context, runtime, task, regions)));
  } // set_state

  void pop_state(
    size_t key
  )
  {
    #ifndef NDEBUG
      std::cout << "popping state for " << key << std::endl;
    #endif

    state_[key].pop();
  } // set_state

#if 0
  /*!
    Reset the legion runtime state.
   */
  void
  set_state(
    lr_context_t & context,
    lr_runtime_t * runtime,
    const lr_task_t * task,
    const lr_regions_t & regions
  )
  {
    state_.reset(new legion_runtime_state_t(context, runtime, task, regions));
  } // set_state
#endif

  /*--------------------------------------------------------------------------*
   * Task registraiton.
   *--------------------------------------------------------------------------*/

  using task_id_t = LegionRuntime::HighLevel::TaskID;
  using register_function_t = std::function<void(size_t)>;
  using unique_fid_t = unique_id_t<task_id_t>;

  /*!
   */
  bool
  register_task(
    task_hash_key_t key,
    const register_function_t & f
  )
  {
    if(task_registry_.find(key) == task_registry_.end()) {
      task_registry_[key] = { unique_fid_t::instance().next(), f };
      return true;
    } // if

    return false;
  } // register_task

  /*!
   */
  task_id_t
  task_id(
    task_hash_key_t key
  )
  {
    assert(task_registry_.find(key) != task_registry_.end() &&
      "task key does not exist!");

    return task_registry_[key].first;
  } // task_id

  /*--------------------------------------------------------------------------*
   * Function registraiton.
   *--------------------------------------------------------------------------*/

  /*!
   */
  template<typename T>
  bool
  register_function(
    const const_string_t & key,
    T & function
  )
  {
    size_t h = key.hash();
    if(function_registry_.find(h) == function_registry_.end()) {
      function_registry_[h] =
        reinterpret_cast<std::function<void(void)> *>(&function);
      return true;
    } // if

    return false;
  } // register_function
  
  /*!
   */
  std::function<void(void)> *
  function(
    size_t key
  )
  {
    return function_registry_[key];
  } // function

  /*--------------------------------------------------------------------------*
   * Legion runtime accessors.
   *--------------------------------------------------------------------------*/

  LegionRuntime::HighLevel::Context &
  context(
    size_t task_key
  )
  {
    return state_[task_key].top()->context;
  } // context

  LegionRuntime::HighLevel::HighLevelRuntime *
  runtime(
    size_t task_key
  )
  {
    return state_[task_key].top()->runtime;
  } // runtime

  const
  LegionRuntime::HighLevel::Task *
  task(
    size_t task_key
  )
  {
    return state_[task_key].top()->task;
  } // task

  const
  std::vector<LegionRuntime::HighLevel::PhysicalRegion> &
  regions(
    size_t task_key
  )
  {
    return state_[task_key].top()->regions;
  } // regions
  
private:

  /*--------------------------------------------------------------------------*
   * Task registry
   *-------------------------------------------------------------------------*/

  // Define the map type using the task_hash_t hash function.
  std::unordered_map<task_hash_t::key_t,
    std::pair<task_id_t, register_function_t>,
    task_hash_t> task_registry_;

  /*--------------------------------------------------------------------------*
   * Function registry
   *--------------------------------------------------------------------------*/

  std::unordered_map<size_t, std::function<void(void)> *>
    function_registry_;

}; // class mpilegion_context_policy_t

} // namespace execution 
} // namespace flecsi

#endif // flecsi_execution_mpilegion_context_policy_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
