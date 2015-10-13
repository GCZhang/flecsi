/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flexi_iterator_h
#define flexi_iterator_h

/*!
 * \file iterator.h
 * \authors bergen
 * \date Initial file creation: Oct 09, 2015
 */

namespace flexi {

/*!
  \class iterator iterator.h
  \brief iterator provides...
 */
template<typename C>
class iterator
{
public:

  using type_t = typename C::type_t;

  //! Default constructor
  iterator(C & values, size_t index)
    : values_(values), index_(index) {}

  //! Destructor
   ~iterator() {}

  //! Copy constructor
  iterator(const iterator & it)
    : values_(it.values_), index_(it.index_) {}

  //! Assignment operator
  iterator & operator = (const iterator & it) {
    index_ = it.index_;
    values_ = it.values_;
    return *this;
  } // operator =

  //! Increment operator
  iterator & operator ++ () { ++index_; return *this; }
  
  //! Dereference operator
  type_t & operator * () { return values_[index_]; }

  //! Dereference operator
  type_t & operator -> () { return values_[index_]; }

  //! Equivalence operator
  bool operator == (const iterator & it) const
    { return index_ == it.index_; }

  //! Comparision operator
  bool operator != (const iterator & it) const
    { return index_ != it.index_; }

private:

  C & values_;
  size_t index_;

}; // class iterator

} // namespace flexi

#endif // flexi_iterator_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
