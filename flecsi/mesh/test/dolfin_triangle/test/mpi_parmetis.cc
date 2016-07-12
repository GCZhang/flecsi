/*~-------------------------------------------------------------------------~~*
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
 *~-------------------------------------------------------------------------~~*/

#include <cinchtest.h>
#include <parmetis.h>

#include "flecsi/execution/mpi_execution_policy.h"

using namespace flecsi;
using namespace testing;

class mpi_parmetis_2way : public Test {
protected:

  virtual void SetUp() override {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    idx_t ncon = 1;
    idx_t nparts = 2;
    idx_t wgtflag = 0;
    idx_t numflag = 0;

    // No, these three arrays can not be NULL
    real_t tpwgts[2] = {0.5, 0.5};
    real_t ubvec[1] = {1.05};
    idx_t options[1] = {0};

    idx_t edgecut;
    MPI_Comm comm = MPI_COMM_WORLD;

    auto ret = ParMETIS_V3_PartKway(
      vtxdist, xadj[rank], adjncy[rank],
      nullptr, /* vwgt */
      nullptr, /* adjwgt */
      &wgtflag,
      &numflag,
      &ncon,
      &nparts,
      tpwgts,
      ubvec,
      options,
      &edgecut,
      part,
      &comm
    );

    ASSERT_EQ(ret, METIS_OK);
  }

  idx_t vtxdist[3] = {0, 5, 10};

  idx_t xadj[2][6] = {
    {0, 2, 4, 6, 9, 11},
    {0, 2, 4, 6, 9, 11}
  };

  idx_t adjncy[2][11] = {
    {1, 9, 2, 0, 3, 1, 4, 8, 2, 5, 3},
    {6, 4, 7, 5, 8, 6, 9, 3, 7, 0, 8}
  };

  idx_t part[5];
};

TEST(mpi_parmetis, comm_size_should_be_2) {
  int comm_size = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  ASSERT_THAT(comm_size, Eq(2));
}

TEST_F(mpi_parmetis_2way, vertices_are_in_partition_0_or_1) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "rank: " << rank << " partition: ";
  for (auto i: part) {
    std::cout << i << " ";
  }
  std::cout << std::endl;
  ASSERT_THAT(part, Each(AnyOf(Eq(0), Eq(1))));
}