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
/*
Description:    This code solves a system corresponding to a discretization
                     of the Laplace equation -Delta u = 1 with zero boundary
                     conditions on the unit square.  The domain is split into
                     an N x N processor grid.  Thus, the given number of processors
                     should be a perfect square.  Each processor's piece of the
                     grid has n x n cells with n x n nodes connected by the
                     standard 5-point stencil. Note that the struct interface
                     assumes a cell-centered grid, and, therefore, the nodes are
                     not shared.  This example demonstrates more features than the
                     previous two struct examples (Example 1 and Example 2).  Two
                     solvers are available.
  
                     To incorporate the boundary conditions, we do the following:
                     Let x_i and x_b be the interior and boundary parts of the
                     solution vector x. We can split the matrix A as
                                     A = [A_ii A_ib; A_bi A_bb].
                     Let u_0 be the Dirichlet B.C.  We can simply say that x_b = u_0.
                     If b_i is the right-hand side, then we just need to solve in
                     the interior:
                                      A_ii x_i = b_i - A_ib u_0.
                     For this partitcular example, u_0 = 0, so we are just solving
                     A_ii x_i = b_i.
  
                     We recommend viewing examples 1 and 2 before viewing this
                     example.

*/
// system includes
#include <cinchtest.h>
#include <iostream>
#include <string>
#include <type_traits> // std::is_same
#include "legion.h"

// user includes
#include "flecsi/utils/mpi_legion_interoperability/mpi_legion_interop.h"
//#include "flecsi/utils/mpi_legion_interoperability/legion_handshake.h"
#include "flecsi/utils/mpi_legion_interoperability/mapper.h"
#include "flecsi/utils/mpi_legion_interoperability/task_ids.h"
#include "flecsi/execution/mpi_execution_policy.h"
#include "flecsi/execution/task.h"

//hypre includes
#include <math.h>
#include "_hypre_utilities.h"
#include "HYPRE_struct_ls.h"
//#include "vis.c"

using namespace flecsi::mpilegion;
using namespace LegionRuntime::HighLevel;
using namespace LegionRuntime::Accessor;
using namespace LegionRuntime::Arrays;

using execution_t = flecsi::execution_t<flecsi::mpi_execution_policy_t>;
using return_type_t = execution_t::return_type_t;

enum TaskIDs{
 TOP_LEVEL_TASK_ID         =0x00000100,
 HANDOFF_TO_MPI_TASK_ID    =0x00000300,
};


enum FieldIDs {
};

ExtLegionHandshake *handshake;


int hypre_test();

void top_level_task(const Task *task,
                    const std::vector<PhysicalRegion> &regions,
                    Context ctx, HighLevelRuntime *runtime)
{
  int num_local_procs=0;
  #ifndef SHARED_LOWLEVEL
  // Only the shared lowlevel runtime needs to iterate over all points
  // on each processor.
  int num_points = 1;
  int num_procs = 0;
  {
   std::set<Processor> all_procs;
   Realm::Machine::get_machine().get_all_processors(all_procs);
   for(std::set<Processor>::const_iterator it = all_procs.begin();
      it != all_procs.end();
      it++){
    if((*it).kind() == Processor::LOC_PROC)
      num_procs++;
   }
  }
  num_local_procs=num_procs;
#else
  int num_procs = Machine::get_machine()->get_all_processors().size();
  int num_points = rank->proc_grid_size.x[0] * rank->proc_grid_size.x[1] * rank->proc_grid_size.x[2];
#endif
  printf("Attempting to connect %d processors with %d points per processor\n",
         num_procs, num_points);
  Point<2> all_procs_lo, all_procs_hi;
  all_procs_lo.x[0] = all_procs_lo.x[1] = 0;
  all_procs_hi.x[0] = num_procs - 1;
  all_procs_hi.x[1] = num_points - 1;
  Rect<2> all_processes = Rect<2>(all_procs_lo, all_procs_hi);

  Rect<1> local_procs(0,num_local_procs);
  ArgumentMap arg_map;

  IndexLauncher connect_mpi_launcher(CONNECT_MPI_TASK_ID,
                                       Domain::from_rect<2>(all_processes),
                                       TaskArgument(0, 0),
                                       arg_map);

 //run legion_init() from each thead
  FutureMap fm1 = runtime->execute_index_space(ctx, connect_mpi_launcher);
   printf("connect_mpi finished \n");
  //run some legion task here 
   fm1.wait_all_results();
  
 //handoff to MPI
 TaskLauncher handoff_to_mpi_launcher(HANDOFF_TO_MPI_TASK_ID,
      TaskArgument(0, 0));
 std::vector<Future> future_tmp;
  future_tmp.push_back(runtime->execute_task(ctx, handoff_to_mpi_launcher));
 
}



void connect_mpi_task (const Task *task,
                      const std::vector<PhysicalRegion> &regions,
                      Context ctx, HighLevelRuntime *runtime)
{
    handshake->legion_init();
}

int handoff_to_mpi_task (const Task *legiontask,
                      const std::vector<PhysicalRegion> &regions,
                      Context ctx, HighLevelRuntime *runtime)
{
 handshake->legion_handoff_to_ext();
 return 0;
}


void complete_legion_configure(void)
{
   handshake->ext_init();
}

void run_legion_task(void)
{
   handshake->ext_handoff_to_legion();
}


void my_init_legion(){

  handshake = new ExtLegionHandshake(ExtLegionHandshake::IN_EXT, 1, 1);

  HighLevelRuntime::set_top_level_task_id(TOP_LEVEL_TASK_ID);

  HighLevelRuntime::register_legion_task<top_level_task>( TOP_LEVEL_TASK_ID,
                          Processor::LOC_PROC, true/*single*/, false/*index*/, 
                          AUTO_GENERATE_ID, TaskConfigOptions(), "top_level_task");
  HighLevelRuntime::register_legion_task< connect_mpi_task>( CONNECT_MPI_TASK_ID, 
                          Processor::LOC_PROC, false/*single*/, true/*index*/,
                          AUTO_GENERATE_ID, TaskConfigOptions(true/*leaf*/), "connect_mpi_task");

  HighLevelRuntime::register_legion_task<int,handoff_to_mpi_task>( HANDOFF_TO_MPI_TASK_ID,
                          Processor::LOC_PROC, true/*single*/, false/*index*/, 0,
                           TaskConfigOptions(), "handoff_to_mpi_task");

  const InputArgs &args = HighLevelRuntime::get_input_args();

  HighLevelRuntime::set_registration_callback(mapper_registration);

  HighLevelRuntime::start(args.argc, args.argv, true);

  complete_legion_configure();

  run_legion_task();  

  handshake->ext_wait_on_legion(); 

  int return_value;
  return_value=hypre_test();
  //copy data to HYPRE and call Hypre executables

}

#define execute(task, ...) \
  execution_t::execute_task(task, ##__VA_ARGS__)

  int hypre_test ()
  {
     int i, j, k;
  
     int myid, num_procs;
  
     int n, N, pi, pj;
     double h, h2;
     int ilower[2], iupper[2];
  
     int solver_id;
     int n_pre, n_post;
  
     HYPRE_StructGrid     grid;
     HYPRE_StructStencil  stencil;
     HYPRE_StructMatrix   A;
     HYPRE_StructVector   b;
     HYPRE_StructVector   x;
     HYPRE_StructSolver   solver;
     HYPRE_StructSolver   precond;
  
     int num_iterations;
     double final_res_norm;
  
     int vis;
  
     /* Initialize MPI */
     //MPI_Init(&argc, &argv);
     MPI_Comm_rank(MPI_COMM_WORLD, &myid);
     MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
 
    printf ("myid = %d \n", myid);
 
     /* Set defaults */
     n = 33;
     solver_id = 0;
     n_pre  = 1;
     n_post = 1;
     vis = 0;
  
 
    /* Figure out the processor grid (N x N).  The local problem
       size for the interior nodes is indicated by n (n x n).
       pi and pj indicate position in the processor grid. */
    N  = sqrt(num_procs);
    h  = 1.0 / (N*n+1); /* note that when calculating h we must
                           remember to count the boundary nodes */
    h2 = h*h;
    pj = myid / N;
    pi = myid - pj*N;
 
    /* Figure out the extents of each processor's piece of the grid. */
    ilower[0] = pi*n;
    ilower[1] = pj*n;
 
    iupper[0] = ilower[0] + n-1;
    iupper[1] = ilower[1] + n-1;
 
    /* 1. Set up a grid */
    {
       /* Create an empty 2D grid object */
       HYPRE_StructGridCreate(MPI_COMM_WORLD, 2, &grid);
 
       /* Add a new box to the grid */
       HYPRE_StructGridSetExtents(grid, ilower, iupper);
 
       /* This is a collective call finalizing the grid assembly.
          The grid is now ``ready to be used'' */
       HYPRE_StructGridAssemble(grid);
    }
 
    /* 2. Define the discretization stencil */
    {
       /* Create an empty 2D, 5-pt stencil object */
       HYPRE_StructStencilCreate(2, 5, &stencil);
 
       /* Define the geometry of the stencil */
       {
          int entry;
          int offsets[5][2] = {{0,0}, {-1,0}, {1,0}, {0,-1}, {0,1}};
 
          for (entry = 0; entry < 5; entry++)
             HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);
       }
    }
 
    /* 3. Set up a Struct Matrix */
    {
       int nentries = 5;
       int nvalues = nentries*n*n;
       double *values;
       int stencil_indices[5];
 
       /* Create an empty matrix object */
       HYPRE_StructMatrixCreate(MPI_COMM_WORLD, grid, stencil, &A);
 
       /* Indicate that the matrix coefficients are ready to be set */
       HYPRE_StructMatrixInitialize(A);
 
       values = new double[nvalues];
 
       for (j = 0; j < nentries; j++)
          stencil_indices[j] = j;
 
       /* Set the standard stencil at each grid point,
          we will fix the boundaries later */
       for (i = 0; i < nvalues; i += nentries)
       {
          values[i] = 4.0;
          for (j = 1; j < nentries; j++)
             values[i+j] = -1.0;
       }
 
       HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries,
                                      stencil_indices, values);
 
       delete [] values;
    }
 
    /* 4. Incorporate the zero boundary conditions: go along each edge of
          the domain and set the stencil entry that reaches to the boundary to
          zero.*/
    {
       int bc_ilower[2];
       int bc_iupper[2];
       int nentries = 1;
       int nvalues  = nentries*n; /*  number of stencil entries times the length
                                      of one side of my grid box */
       double *values;
       int stencil_indices[1];
 
       values = new double[nvalues];
       for (j = 0; j < nvalues; j++)
          values[j] = 0.0; 

       /* Recall: pi and pj describe position in the processor grid */
       if (pj == 0)
       {
          /* Bottom row of grid points */
          bc_ilower[0] = pi*n;
          bc_ilower[1] = pj*n;
 
          bc_iupper[0] = bc_ilower[0] + n-1;
          bc_iupper[1] = bc_ilower[1];
 
          stencil_indices[0] = 3;
 
          HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                         stencil_indices, values);
       }
 
       if (pj == N-1)
       {
          /* upper row of grid points */
          bc_ilower[0] = pi*n;
          bc_ilower[1] = pj*n + n-1;
 
          bc_iupper[0] = bc_ilower[0] + n-1;
          bc_iupper[1] = bc_ilower[1];
 
          stencil_indices[0] = 4;
 
          HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                         stencil_indices, values);  
    }
 
       if (pi == 0)
       {
          /* Left row of grid points */
          bc_ilower[0] = pi*n;
          bc_ilower[1] = pj*n;
 
          bc_iupper[0] = bc_ilower[0];
          bc_iupper[1] = bc_ilower[1] + n-1;
 
          stencil_indices[0] = 1;
 
          HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                         stencil_indices, values);
       }
 
       if (pi == N-1)
       {
          /* Right row of grid points */
          bc_ilower[0] = pi*n + n-1;
          bc_ilower[1] = pj*n;
 
          bc_iupper[0] = bc_ilower[0];
          bc_iupper[1] = bc_ilower[1] + n-1;
 
          stencil_indices[0] = 2;
 
          HYPRE_StructMatrixSetBoxValues(A, bc_ilower, bc_iupper, nentries,
                                         stencil_indices, values);
       }
 
       delete [] values;
    } 

    /* This is a collective call finalizing the matrix assembly.
       The matrix is now ``ready to be used'' */
    HYPRE_StructMatrixAssemble(A);

    /* 5. Set up Struct Vectors for b and x */
    {
       int    nvalues = n*n;
       double *values;
 
       values = new double[nvalues];
 
       /* Create an empty vector object */
       HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &b);
       HYPRE_StructVectorCreate(MPI_COMM_WORLD, grid, &x);
 
       /* Indicate that the vector coefficients are ready to be set */
       HYPRE_StructVectorInitialize(b);
       HYPRE_StructVectorInitialize(x);
 
      /* Set the values */
       for (i = 0; i < nvalues; i ++)
          values[i] = h2;
       HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
 
       for (i = 0; i < nvalues; i ++)
          values[i] = 0.0;
       HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
 
       delete [] values;
 
       /* This is a collective call finalizing the vector assembly.
          The vector is now ``ready to be used'' */
       HYPRE_StructVectorAssemble(b);
       HYPRE_StructVectorAssemble(x);
    }
 
    /* 6. Set up and use a struct solver
       (Solver options can be found in the Reference Manual.) */
    if (solver_id == 0)
    {
       HYPRE_StructPCGCreate(MPI_COMM_WORLD, &solver);
       HYPRE_StructPCGSetMaxIter(solver, 50 );
       HYPRE_StructPCGSetTol(solver, 1.0e-06 );
       HYPRE_StructPCGSetTwoNorm(solver, 1 );
       HYPRE_StructPCGSetRelChange(solver, 0 );
       HYPRE_StructPCGSetPrintLevel(solver, 2 ); /* print each CG iteration */
       HYPRE_StructPCGSetLogging(solver, 1);
 
       /* Use symmetric SMG as preconditioner */
       HYPRE_StructSMGCreate(MPI_COMM_WORLD, &precond);
       HYPRE_StructSMGSetMemoryUse(precond, 0);
       HYPRE_StructSMGSetMaxIter(precond, 1);
       HYPRE_StructSMGSetTol(precond, 0.0);
       HYPRE_StructSMGSetZeroGuess(precond);
       HYPRE_StructSMGSetNumPreRelax(precond, 1);
       HYPRE_StructSMGSetNumPostRelax(precond, 1);
 
       /* Set the preconditioner and solve */
       HYPRE_StructPCGSetPrecond(solver, HYPRE_StructSMGSolve,
                                   HYPRE_StructSMGSetup, precond);
       HYPRE_StructPCGSetup(solver, A, b, x);
       HYPRE_StructPCGSolve(solver, A, b, x);
 
       /* Get some info on the run */
       HYPRE_StructPCGGetNumIterations(solver, &num_iterations);
       HYPRE_StructPCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
 
       /* Clean up */
       HYPRE_StructPCGDestroy(solver);
    }
 
    if (solver_id == 1)
    {
       HYPRE_StructSMGCreate(MPI_COMM_WORLD, &solver);
       HYPRE_StructSMGSetMemoryUse(solver, 0);
       HYPRE_StructSMGSetMaxIter(solver, 50);
       HYPRE_StructSMGSetTol(solver, 1.0e-06);
       HYPRE_StructSMGSetRelChange(solver, 0);
       HYPRE_StructSMGSetNumPreRelax(solver, n_pre);
       HYPRE_StructSMGSetNumPostRelax(solver, n_post);
       /* Logging must be on to get iterations and residual norm info below */
       HYPRE_StructSMGSetLogging(solver, 1);
 
       /* Setup and solve */
       HYPRE_StructSMGSetup(solver, A, b, x);
       HYPRE_StructSMGSolve(solver, A, b, x);
 
       /* Get some info on the run */
       HYPRE_StructSMGGetNumIterations(solver, &num_iterations);
       HYPRE_StructSMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
 
       /* Clean up */
       HYPRE_StructSMGDestroy(solver);
    }
 
    /* Save the solution for GLVis visualization, see vis/glvis-ex3.sh */
    if (vis)
    {
       FILE *file;
       char filename[255];
 
       int nvalues = n*n;
       double *values;
       values = new double[nvalues];
 
       /* get the local solution */
       HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
 
       sprintf(filename, "%s.%06d", "vis/ex3.sol", myid);
       if ((file = fopen(filename, "w")) == NULL)
       {
          printf("Error: can't open output file %s\n", filename);
          //MPI_Finalize();
          exit(1);
       }
 
       /* save solution with global unknown numbers */
       k = 0;
       for (j = 0; j < n; j++)
          for (i = 0; i < n; i++)
             fprintf(file, "%06d %.14e\n", pj*N*n*n+pi*n+j*N*n+i, values[k++]);
 
       fflush(file);
       fclose(file);
       delete [] values;
 
       /* save global finite element mesh */
    //   if (myid == 0)
    //TOFIX      GLVis_PrintGlobalSquareMesh("vis/ex3.mesh", N*n-1);
    }
 
    if (myid == 0)
    {
       printf("\n");
       printf("Iterations = %d\n", num_iterations);
       printf("Final Relative Residual Norm = %g\n", final_res_norm);
       printf("\n");
    }
 
    /* Free memory */
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructStencilDestroy(stencil);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
 
    /* Finalize MPI */
    //MPI_Finalize();
 
    return (0);
 }


TEST(more_complex_mpi_with_legion, simple) {
 //  ASSERT_LT(execute(world_size), 1);
  
   int return_value;

   my_init_legion(); 
 
} // TEST



/*----------------------------------------------------------------------------*
 * Cinch test Macros
 *
 *  ==== I/O ====
 *  CINCH_CAPTURE()              : Insertion stream for capturing output.
 *                                 Captured output can be written or
 *                                 compared using the macros below.
 *
 *    EXAMPLE:
 *      CINCH_CAPTURE() << "My value equals: " << myvalue << endl;
 *
 *  CINCH_COMPARE_BLESSED(file); : Compare captured output with
 *                                 contents of a blessed file.
 *
 *  CINCH_WRITE(file);           : Write captured output to file.
 *
 * Google Test Macros
 *
 * Basic Assertions:
 *
 *  ==== Fatal ====             ==== Non-Fatal ====
 *  ASSERT_TRUE(condition);     EXPECT_TRUE(condition)
 *  ASSERT_FALSE(condition);    EXPECT_FALSE(condition)
 *
 *
 * Binary Comparison:
 *
 *  ==== Fatal ====             ==== Non-Fatal ====
 *  ASSERT_EQ(val1, val2);      EXPECT_EQ(val1, val2)
 *  ASSERT_NE(val1, val2);      EXPECT_NE(val1, val2)
 *  ASSERT_LT(val1, val2);      EXPECT_LT(val1, val2)
 *  ASSERT_LE(val1, val2);      EXPECT_LE(val1, val2)
 *  ASSERT_GT(val1, val2);      EXPECT_GT(val1, val2)
 *  ASSERT_GE(val1, val2);      EXPECT_GE(val1, val2)
 *
 * String Comparison:
 *
 *  ==== Fatal ====                     ==== Non-Fatal ====
 *  ASSERT_STREQ(expected, actual);     EXPECT_STREQ(expected, actual)
 *  ASSERT_STRNE(expected, actual);     EXPECT_STRNE(expected, actual)
 *  ASSERT_STRCASEEQ(expected, actual); EXPECT_STRCASEEQ(expected, actual)
 *  ASSERT_STRCASENE(expected, actual); EXPECT_STRCASENE(expected, actual)
 *----------------------------------------------------------------------------*/

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
