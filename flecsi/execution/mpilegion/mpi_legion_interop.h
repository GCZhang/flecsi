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

#ifndef mpi_legion_interop_h
#define mpi_legion_interop_h

#include <iostream>
#include <string>
#include <cstdio>
#include <mutex>
#include <condition_variable>

#include <mpi.h>
#include <legion.h>
#include <realm.h>

#include "flecsi/execution/mpilegion/legion_handshake.h"
#include "flecsi/execution/mpilegion/mapper.h"
#include "flecsi/execution/mpilegion/task_ids.h"
#include "flecsi/partition/index_partition.h"
#include "flecsi/utils/any.h"

///
// \file mpilegion/mpi_legion_interop.h
// \authors demeshko
// \date Initial file creation: Jul 2016
///

namespace flecsi{
namespace execution{


//----------------------------------------------------------------------------//
///
// a legion task that created User Events/Queques for synchronization 
// between MPI and Legion. This task is an index task that will be executed on
// as many processes as MPI has been called from.
///
//----------------------------------------------------------------------------//

inline
void
connect_to_mpi_task(
  const Legion::Task * legiontask,
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions,
  LegionRuntime::HighLevel::Context ctx,
  LegionRuntime::HighLevel::HighLevelRuntime * runtime
)
{
  #ifdef LEGIONDEBUG
  std::cout <<"inside connect_to_mpi"<<std::endl;
  #endif
  
  // this call creates Legion events/queues for both MPI and Legion runtimes
  // for the later synchronization and then swithces runtime MPI
  ext_legion_handshake_t::instance().legion_init();
} // connect_to_mpi_task

//----------------------------------------------------------------------------//
///
// a legion task that calls legion_handoff_to_ext function from handshake
// class that switches form Legion to MPI runtime. This task is an 
// index task that will be executed on as many processes as MPI 
// has been called from.
///
//----------------------------------------------------------------------------//

inline
void
handoff_to_mpi_task(
  const Legion::Task * legiontask,
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions,
  LegionRuntime::HighLevel::Context ctx,
  LegionRuntime::HighLevel::HighLevelRuntime * runtime)
{
   ext_legion_handshake_t::instance().legion_handoff_to_ext();
} // handoff_to_mpi_task

//----------------------------------------------------------------------------//
///
// a legion task that waits for every MPI thread to finish
///
//----------------------------------------------------------------------------//

inline
void
wait_on_mpi_task(
  const Legion::Task * legiontask,
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> & regions,
  LegionRuntime::HighLevel::Context ctx,
  LegionRuntime::HighLevel::HighLevelRuntime * runtime
)
{
  ext_legion_handshake_t::instance().legion_wait_on_ext();
} // wait_on_mpi_task

//----------------------------------------------------------------------------//
///
// In purpose to tell MPI runtime that there is an MPI task that has to be 
// executed, we have to switch call_mpi boolian variable on every MPI thread.
// his task is an index task that will be executed on as many processes as MPI 
// has been called from.
///
//----------------------------------------------------------------------------//

inline
void
unset_call_mpi_task(
  const Legion::Task *legiontask,
  const std::vector<LegionRuntime::HighLevel::PhysicalRegion> &regions,
  LegionRuntime::HighLevel::Context ctx,
  LegionRuntime::HighLevel::HighLevelRuntime *runtime
)
{
  ext_legion_handshake_t::instance().call_mpi_=false;
} // unset_call_mpi_task


//----------------------------------------------------------------------------//
///
// mpi_legion_interop class is used to wrap legion-mpi handshaking routines 
// in a FLeCSI-friendly interface
///
//----------------------------------------------------------------------------//

struct mpi_legion_interop_t
{
  ///
  // Constructor.
  ///
  mpi_legion_interop_t() {};

  ///
  // Copy constructor.
  ///
  mpi_legion_interop_t(const mpi_legion_interop_t &) = delete;

  ///
  // Assign operator.
  ///
  mpi_legion_interop_t& operator=(const mpi_legion_interop_t &) = delete;

  /// 
  // Destructor.
  //
  ~mpi_legion_interop_t() {};

  ///
  // Initialze() method needs to be called before we start Legion runtime. 
  // The method load a customized mapper that is tuned for MPI-Legion 
  // interoperability.
  ///
  void
  initialize()
  {
    LegionRuntime::HighLevel::HighLevelRuntime::set_registration_callback(
      mapper_registration);
  } // initialize
 
  ///
  // This method is used to communicate to  MPI runtime that there is an
  // MPI task that has to be executed after we switches from Legion to 
  // MPI. It launches legion Index task that switches call_mpi variable
  // to true for every MPI process
  ///
  void
  unset_call_mpi(
   LegionRuntime::HighLevel::Context ctx,
   LegionRuntime::HighLevel::HighLevelRuntime *runtime
  )
  {
  	LegionRuntime::HighLevel::ArgumentMap arg_map;

  	LegionRuntime::HighLevel::IndexLauncher unset_call_mpi_launcher(
    	task_ids_t::instance().unset_call_mpi_id,
    	LegionRuntime::HighLevel::Domain::from_rect<2>(all_processes_),
    	LegionRuntime::HighLevel::TaskArgument(0, 0),
    	arg_map);

  	unset_call_mpi_launcher.tag =  MAPPER_ALL_PROC;

  	LegionRuntime::HighLevel::FutureMap fm5 =
    runtime->execute_index_space(ctx,unset_call_mpi_launcher);

  	fm5.wait_all_results();
	}//unset_call_mpi
 
  ///
  // This method creates pthreads mutex on the MPI side and, in case handshake 
  // is originally created in MPI, waits on when handshake is created on
  // the Legion side.
  ///
  void
  legion_configure()
  {
    ext_legion_handshake_t::instance().ext_init();
  } // legion_configure

  ///
  // Swithches mutex to Legion runtime.
  ///
  void 
  handoff_to_legion()
  {
    ext_legion_handshake_t::instance().ext_handoff_to_legion();
  } // handoff_to_legion

  ///
  // Wait untill mutex switched to MPI runtime and run MPI.
  ///
  void 
  wait_on_legion()
  {
    ext_legion_handshake_t::instance().ext_wait_on_legion();
  } // wait_on_legion

  ///
  // Register all Legion tasks used in the mpi_legion_interop_t class.
  ///
  static
	void 
	register_tasks()
	{
  	LegionRuntime::HighLevel::HighLevelRuntime::register_legion_task
	    <connect_to_mpi_task>(
	    task_ids_t::instance().connect_mpi_task_id,
	    LegionRuntime::HighLevel::Processor::LOC_PROC,
	    false/*single*/, true/*index*/,
	    AUTO_GENERATE_ID,
	    LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
	    "connect_to_mpi_task");

  	LegionRuntime::HighLevel::HighLevelRuntime::register_legion_task
    	<handoff_to_mpi_task>(
	    task_ids_t::instance().handoff_to_mpi_task_id,
    	LegionRuntime::HighLevel::Processor::LOC_PROC,
    	false/*single*/, true/*index*/,
    	AUTO_GENERATE_ID,
    	LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
    	"handoff_to_mpi_task");

  	LegionRuntime::HighLevel::HighLevelRuntime::register_legion_task
    	<wait_on_mpi_task>(
    	task_ids_t::instance().wait_on_mpi_task_id,
    	LegionRuntime::HighLevel::Processor::LOC_PROC,
    	false/*single*/, true/*index*/,
    	AUTO_GENERATE_ID,
    	LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
    	"wait_on_mpi_task");

   	LegionRuntime::HighLevel::HighLevelRuntime::register_legion_task
     	<unset_call_mpi_task>(
     	task_ids_t::instance().unset_call_mpi_id,
     	LegionRuntime::HighLevel::Processor::LOC_PROC,
     	false, true, AUTO_GENERATE_ID,
     	LegionRuntime::HighLevel::TaskConfigOptions(true/*leaf*/),
     	"unset_call_mpi_task");

	}//register_tasks

  ///
  // Helper function that calculates number of availabe processes (the number
  // of processes that MPI has been called from).
  ///
  void 
  calculate_number_of_procs ()
  {
		int num_local_procs=0;
		#ifndef SHARED_LOWLEVEL
		// Only the shared lowlevel runtime needs to iterate over all points
		// on each processor.
		int num_points = 1;
		int num_procs = 0;
		{
			std::set<LegionRuntime::HighLevel::Processor> all_procs;
			Realm::Machine::get_machine().get_all_processors(all_procs);
   		for(std::set<LegionRuntime::HighLevel::Processor>::const_iterator
      	it = all_procs.begin();
      	it != all_procs.end();
     		it++)
			{
          if((*it).kind() == LegionRuntime::HighLevel::Processor::LOC_PROC)
               num_procs++;
      }//end for
  	} // scope
  	num_local_procs=num_procs;
#else
  	int num_procs = LegionRuntime::HighLevel::Machine::get_machine()->
                                            get_all_processors().size();
  	int num_points = rank->proc_grid_size.x[0] * rank->proc_grid_size.x[1]
                                            * rank->proc_grid_size.x[2];
#endif
  	printf("Attempting to connect %d processors with %d points per processor\n",
         num_procs, num_points);
  	LegionRuntime::Arrays::Point<2> all_procs_lo, all_procs_hi;
  	all_procs_lo.x[0] = all_procs_lo.x[1] = 0;
  	all_procs_hi.x[0] = num_procs - 1;
  	all_procs_hi.x[1] = num_points - 1;
  	this->all_processes_ =  LegionRuntime::Arrays::Rect<2>(all_procs_lo,
                                                        all_procs_hi);
  	this->local_procs_ = LegionRuntime::Arrays::Rect<1>(0,num_local_procs);
  	std::cout << all_procs_lo.x[0] << "," << all_procs_lo.x[1] << ","
            << all_procs_hi.x[0] << "," << all_procs_hi.x[1] << std::endl;

	}//calculate_number_of_procs

  ///
  // Calls a legion's index task that created User Events/Queques for 
  // synchronization between MPI and Legion. 
  ///
  void 
  connect_with_mpi(
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *runtime
  )
	{
		calculate_number_of_procs();

  	LegionRuntime::HighLevel::ArgumentMap arg_map;

  	LegionRuntime::HighLevel::IndexLauncher connect_mpi_launcher(
    	task_ids_t::instance().connect_mpi_task_id,
    	LegionRuntime::HighLevel::Domain::from_rect<2>(all_processes_),
    	LegionRuntime::HighLevel::TaskArgument(0, 0),
   	  arg_map);

  	//run legion_init() from each thead
  	LegionRuntime::HighLevel::FutureMap fm1 =
    	runtime->execute_index_space(ctx, connect_mpi_launcher);

  	//run some legion task here
  	fm1.wait_all_results();
	} // connect_with_mpi


  ///
  // This method calls an Index task that switches form Legion to MPI runtime.
  ///
  void 
  handoff_to_mpi(
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *runtime
  )
  {
	 	LegionRuntime::HighLevel::ArgumentMap arg_map;

  	LegionRuntime::HighLevel::IndexLauncher handoff_to_mpi_launcher(
    	task_ids_t::instance().handoff_to_mpi_task_id,
    	LegionRuntime::HighLevel::Domain::from_rect<2>(all_processes_),
    	LegionRuntime::HighLevel::TaskArgument(0, 0),
  	  arg_map);

 	 LegionRuntime::HighLevel::FutureMap fm2 =
    	runtime->execute_index_space(ctx, handoff_to_mpi_launcher);

  	fm2.wait_all_results();
	} // handoff_to_mpi
 
  ///
  // This method waits on all MPI processes to finish all tasks and
  // switch to Legion runtime
  ///
  LegionRuntime::HighLevel::FutureMap 
  wait_on_mpi(
    LegionRuntime::HighLevel::Context ctx,
    LegionRuntime::HighLevel::HighLevelRuntime *runtime
  )
	{
   	LegionRuntime::HighLevel::ArgumentMap arg_map;
  	LegionRuntime::HighLevel::IndexLauncher wait_on_mpi_launcher(
    	task_ids_t::instance().wait_on_mpi_task_id,
    	LegionRuntime::HighLevel::Domain::from_rect<2>(all_processes_),
    	LegionRuntime::HighLevel::TaskArgument(0, 0),
    	arg_map);

  	LegionRuntime::HighLevel::FutureMap fm3 =
    	runtime->execute_index_space(ctx, wait_on_mpi_launcher);

  	fm3.wait_all_results();

  	return fm3;
	}//wait_on_mpi
 
  std::function<void()> shared_func_;

  bool call_mpi_ = false;

  Rect<2> all_processes_;
  Rect<1> local_procs_;

  std::vector<flecsi::utils::any_t> data_storage_;

}; // mpi_legion_interop_t

} // namespace execution
} // namespace flecsi

#endif // mpi_legion_interop_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
