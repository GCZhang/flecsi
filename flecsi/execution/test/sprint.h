/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sprint_h
#define flecsi_sprint_h

#include <iostream>

#include "flecsi/utils/common.h"
#include "flecsi/execution/context.h"
#include "flecsi/execution/execution.h"
#include "flecsi/partition/index_partition.h"
#include "flecsi/partition/init_partitions_task.h"

///
// \file sprint.h
// \authors bergen
// \date Initial file creation: Aug 23, 2016
///

using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace LegionRuntime::Arrays;

namespace flecsi {
namespace execution {

using index_partition_t = dmp::index_partition__<size_t>;

static const size_t N = 8;

enum FieldIDs {
  FID_CELL_PART,
};

  
void mpi_task(double val) {
  int rank = 0;
  int size = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  //std::cout << "My rank: " << rank << std::endl;

  size_t part = N/size;
  size_t rem = N%size;

  size_t start = rank*(part + (rem > 0 ? 1 : 0));
  size_t end = rank < rem ? start + part+1 : start + part;

  
                      
#if 1
  std::cout << "rank: " << rank << " start: " << start <<
    " end: " << end << std::endl;
#endif

  index_partition_t ip;

  for(size_t j(0); j<N; ++j) {
    for(size_t i(start); i<end; ++i) {
      const size_t id = j*N+i;
      // exclusive
      if(i>start && i<end-1) {
        ip.exclusive.push_back(id);
        //std::cout << "rank: " << rank << " exclusive: " << id << std::endl;
      }
      else if(rank == 0 && i==start) {
        ip.exclusive.push_back(id);
        //std::cout << "rank: " << rank << " exclusive: " << id << std::endl;
      }
      else if(rank == size-1 && i==end-1) {
        ip.exclusive.push_back(id);
        //std::cout << "rank: " << rank << " exclusive: " << id << std::endl;
      }
      else if(i==start) {
        ip.shared.push_back(id);
        //std::cout << "rank: " << rank << " shared: " << id << std::endl;

        const size_t ghost_id = j*N+i-1;
        ip.ghost.push_back(ghost_id);
        //std::cout << "rank: " << rank << " ghost: " << ghost_id << std::endl;
      }
      else if(i==end-1) {
        ip.shared.push_back(id);
        //std::cout << "rank: " << rank << " shared: " << id << std::endl;

        const size_t ghost_id = j*N+i+1;
        ip.ghost.push_back(ghost_id);
        //std::cout << "rank: " << rank << " ghost: " << ghost_id << std::endl;
      } // if

    } // for
  } // for
 

  flecsi::execution::context_t & context_ =
             flecsi::execution::context_t::instance();
  context_.interop_helper_.data_storage_.push_back(
        flecsi::utils::any_t(ip));
  

} // mpi_task
  
register_task(mpi_task, mpi, single);
  

void driver(int argc, char ** argv) {
  context_t & context_ = context_t::instance();
  size_t task_key = const_string_t{"driver"}.hash();
  auto runtime = context_.runtime(task_key);
  auto context = context_.context(task_key);

  flecsi::dmp::parts partitions;
  
  // first execute mpi task to setup initial partitions 
  execute_task(mpi_task, mpi, single, 1.0);
  

  // create a field space to store my cell paritioning 
  FieldSpace fs =
    runtime->create_field_space(context);
  {
    FieldAllocator allocator = 
      runtime->create_field_allocator(context, fs);
    allocator.allocate_field(sizeof(int), FID_CELL_PART);
  }

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks); 
  int max_cells = 1000;
  Point<2> elem_rect_lo; elem_rect_lo.x[0] = 0; elem_rect_lo.x[1]=0;
  Point<2> elem_rect_hi;
  elem_rect_hi.x[0] = num_ranks;
  elem_rect_hi.x[1] = max_cells;
  Rect<2> elem_rect( elem_rect_lo, elem_rect_hi );
  
  IndexSpace is = runtime->create_index_space(context,
    Domain::from_rect<2>(elem_rect));
  
  LogicalRegion cell_parts_lr = runtime->create_logical_region(context, is, fs);
  runtime->attach_name(cell_parts_lr, "cell partitions logical region");

  // Use blockify to create num_rank slices of the data  
  Point<2> cell_part_color; cell_part_color.x[0] = 1; cell_part_color.x[1] = max_cells;
  Blockify<2> coloring(cell_part_color);
  
  IndexPartition ip = runtime->create_index_partition(context, is, coloring);
  LogicalPartition cell_parts_lp = runtime->
    get_logical_partition(context, cell_parts_lr, ip);
  runtime->attach_name(cell_parts_lp, "cell partitions logical partition");
  

  LegionRuntime::HighLevel::ArgumentMap arg_map;

  // this is the init partitions task 
  LegionRuntime::HighLevel::IndexLauncher init_cell_partitions_launcher(
    task_ids_t::instance().init_cell_partitions_task_id,
    LegionRuntime::HighLevel::Domain::from_rect<2>(context_.interop_helper_.all_processes_), 
    LegionRuntime::HighLevel::TaskArgument(0, 0),
    arg_map);
  
  
  init_cell_partitions_launcher.add_region_requirement(
    RegionRequirement(cell_parts_lp, 0/*projection ID*/,
                      WRITE_DISCARD, EXCLUSIVE, cell_parts_lr));
  init_cell_partitions_launcher.add_field(0, FID_CELL_PART);
  FutureMap fm = runtime->execute_index_space(context, init_cell_partitions_launcher);
  
  fm.wait_all_results();

  std::cout << "returned from wait_all_results()" << std::endl;

  for (int i = 0; i < num_ranks; i++) {
    std::cout << "about to call get_results" << std::endl; 
    flecsi::dmp::parts received = fm.get_result<flecsi::dmp::parts>(DomainPoint::from_point<2>(make_point(i,0)));
    std::cout << "From rank " << i << " received (exclusive, shared, ghost) "
              << "(" << received.exclusive << "," << received.shared << ","
              << received.ghost << ")" << std::endl; 
      
  }
  
  // GMS: back at the top level, we need to read this partitioning information
  RegionRequirement req(cell_parts_lr, READ_WRITE, EXCLUSIVE, cell_parts_lr);
  req.add_field(FID_CELL_PART);

  
  std::cout << "Back in driver (TTL) and checking values in LR" << std::endl;
  InlineLauncher cell_parts_launcher(req);
  PhysicalRegion cell_parts_region = runtime->map_region(context, cell_parts_launcher);
  cell_parts_region.wait_until_valid();
  RegionAccessor<AccessorType::Generic, int> acc_cell_part =
    cell_parts_region.get_field_accessor(FID_CELL_PART).typeify<int>();
  int zeros = 0;
  int non_zeros = 0; 
  int num_cells = 0; 
  for (int i = 0; i < num_ranks; i++) {
    flecsi::dmp::parts received = fm.get_result<flecsi::dmp::parts>(DomainPoint::from_point<2>(make_point(i,0)));
    int j = 0;
    num_cells += received.exclusive+received.shared; 
    for (; j < received.exclusive; j++) {
      double value =
        acc_cell_part.read(DomainPoint::from_point<2>(make_point(i,j)));
      std::cout << "partition (" << i << "," << j << ") exclusive cell id: " << value << std::endl;
    }
    for (; j < received.exclusive+received.shared; j++) {
      double value =
        acc_cell_part.read(DomainPoint::from_point<2>(make_point(i,j)));
      std::cout << "partition (" << i << "," << j << ") shared cell id: " << value << std::endl;
    }
    for (; j < received.exclusive+received.shared+received.ghost; j++) {
      double value =
        acc_cell_part.read(DomainPoint::from_point<2>(make_point(i,j)));
      std::cout << "partition (" << i << "," << j << ") ghost cell id: " << value << std::endl;
      
    }
  }    
  std::cout << "partition has " << zeros << " zeros " << non_zeros << " non zeros " << std::endl;
  

  // Now I need to create an Index space based on this data ..

#if 0 
  
  Rect<1> elem_rect_cells(Point<1>(0),Point<1>(num_cells-1));
  IndexSpace is_cells = runtime->create_index_space(context, 
                                              Domain::from_rect<2>(elem_rect_cells));
  FieldSpace input_fs = runtime->create_field_space(ctx);
  {
    FieldAllocator allocator = 
      runtime->create_field_allocator(ctx, input_fs);
    allocator.allocate_field(sizeof(double),FID_PRES);
  }

  IndexPartition ip_cells_ex;
  DomainColoring coloring_ex;

  color ... 
  

#endif
    

  
  // Then partition the cells based on the values..

  // etc.. 

  
  
  
} // driver

} // namespace execution
} // namespace flecsi

#endif // flecsi_sprint_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
