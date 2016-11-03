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
/*!
 *
 * \file set_intersect.h
 *
 * \brief detect set intersections.
 *
 ******************************************************************************/
#pragma once

// system includes
#include <algorithm>


namespace flecsi
{
namespace utils
{

////////////////////////////////////////////////////////////////////////////////
//! \brief Determine if two sets intersect.
//! 
//! This is very similar to the way std::includes works.  It just determines
//! whether two sets intersect and returns a boolean result. It does not
//! compute the actual intersection.
//!
//! \param [in] first1,last1  The first sorted range of elements to be examined.
//! \param [in] first2,last2  The second sorted range of elements to be 
//!                           examined.
//!
//! \remark This function has complexity O(n + m)
////////////////////////////////////////////////////////////////////////////////
template<class InputIt1, class InputIt2>
bool intersects(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2)
{
  while (first1 != last1 && first2 != last2) {
    if (*first1 < *first2) {
      ++first1;
      continue;
    } 
    if (*first2 < *first1) {
      ++first2;
      continue;
    } 
    return true;
  }
  return false;
}

////////////////////////////////////////////////////////////////////////////////
//! \brief Determine if two sets intersect.
//! 
//! This is very similar to the way std::includes works.  It just determines
//! whether two sets intersect and returns a boolean result. It does not
//! compute the actual intersection.
//!
//! \param [in] first1,last1  The first sorted range of elements to be examined.
//! \param [in] first2,last2  The second sorted range of elements to be 
//!                           examined.
//! \param [in] comp  A comparison function object which returns true if the
//!                   first argument is less than (i.e. ordered before) the 
//!                   second.
//!
//! \remark This function has complexity O(n + m)
////////////////////////////////////////////////////////////////////////////////
template<class InputIt1, class InputIt2, class Compare>
bool intersects(
  InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2,
  Compare comp
)
{
  while (first1 != last1 && first2 != last2) {
    if ( comp(*first1, *first2) ) {
      ++first1;
      continue;
    } 
    if ( comp(*first2 < *first1) ) {
      ++first2;
      continue;
    } 
    return true;
  }
  return false;
}

#if 0
////////////////////////////////////////////////////////////////////////////////
//! \brief detect intersecitions of sorted lists
//!
//! \remark  When input1 is much smaller that input2, this gives O(n * log(m)) 
//!          time.
////////////////////////////////////////////////////////////////////////////////
template<class InputIt1, class InputIt2>
bool intersects(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2) 
{
  while (first1 != last1)
    if (std::binary_search(first2, last2, *first1++))
      return true;
  return false;
}
#endif

} // namespace
} // namespace

/*~------------------------------------------------------------------------~--*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
