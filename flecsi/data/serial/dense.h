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

#ifndef flecsi_serial_dense_h
#define flecsi_serial_dense_h

//----------------------------------------------------------------------------//
// POLICY_NAMESPACE must be defined before including storage_type.h!!!
// Using this approach allows us to have only one storage_type_t
// definintion that can be used by all data policies -> code reuse...
#define POLICY_NAMESPACE serial
#include "flecsi/data/storage_type.h"
#undef POLICY_NAMESPACE
//----------------------------------------------------------------------------//

#include <memory>

#include "flecsi/data/data_client.h"
#include "flecsi/data/data_handle.h"
#include "flecsi/utils/const_string.h"
#include "flecsi/utils/index_space.h"

///
// \file serial/dense.h
// \authors bergen
// \date Initial file creation: Apr 7, 2016
///

namespace flecsi {
namespace data {
namespace serial {

//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//
// Helper type definitions.
//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//

//----------------------------------------------------------------------------//
// Dense accessor.
//----------------------------------------------------------------------------//

///
// \brief dense_accessor_t provides logically array-based access to data
//        variables that have been registered in the data model.
//
// \tparam T The type of the data variable. If this type is not
//           consistent with the type used to register the data, bad things
//           can happen. However, it can be useful to reinterpret the type,
//           e.g., when writing raw bytes. This class is part of the
//           low-level \e flecsi interface, so it is assumed that you
//           know what you are doing...
// \tparam MD The meta data type.
///
template<typename T, typename MD>
struct dense_accessor_t
{
  //--------------------------------------------------------------------------//
  // Type definitions.
  //--------------------------------------------------------------------------//

  using iterator_t = index_space_t::iterator_t;
  using meta_data_t = MD;
  using user_meta_data_t = typename meta_data_t::user_meta_data_t;

  //--------------------------------------------------------------------------//
  // Constructors.
  //--------------------------------------------------------------------------//

  ///
  // Default constructor.
  ///
  dense_accessor_t() {}

  ///
  // Constructor.
  //
  // \param label The c_str() version of the const_string_t used for
  //              this data variable's hash.
  // \param size The size of the associated index space.
  // \param data A pointer to the raw data.
  // \param user_meta_data A reference to the user-defined meta data.
  ///
  dense_accessor_t(
    const std::string & label,
    const size_t size,
    T * data,
    const user_meta_data_t & user_meta_data,
    bitset_t & user_attributes,
    size_t index_space
  )
  :
    label_(label),
    size_(size),
    data_(data),
    user_meta_data_(user_meta_data),
    user_attributes_(user_attributes),
    index_space_(index_space),
    is_(size)
  {}

	///
  // Copy constructor.
	///
	dense_accessor_t(
    const dense_accessor_t & a
  )
  :
    label_(a.label_),
    size_(a.size_),
    data_(a.data_),
    user_meta_data_(a.user_meta_data_),
    user_attributes_(a.user_attributes_),
    index_space_(a.index_space_),
    is_(a.is_)
  {}

  //--------------------------------------------------------------------------//
  // Member data interface.
  //--------------------------------------------------------------------------//

	///
  // \brief Return a std::string containing the label of the data variable
  //        reference by this accessor.
	///
  const std::string &
  label() const
  {
    return label_;
  } // label

	///
  // \brief Return the index space size of the data variable
  //        referenced by this accessor.
	///
  size_t
  size() const
  {
    return size_;
  } // size

	///
  // \brief Return the user meta data for this data variable.
	///
  const user_meta_data_t &
  meta_data() const
  {
    return user_meta_data_;
  } // meta_data

  ///
  //
  ///
  bitset_t &
  attributes()
  {
    return user_attributes_;
  } // attributes

  ///
  //
  ///
  size_t
  index_space() const
  {
    return index_space_;
  } // index_space

  //--------------------------------------------------------------------------//
  // Iterator interface.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  iterator_t
  begin()
  {
    return {is_, 0};
  } // begin

  ///
  //
  ///
  iterator_t
  end()
  {
    return {is_, size_};
  } // end

  //--------------------------------------------------------------------------//
  // Operators.
  //--------------------------------------------------------------------------//

	///
  // \brief Provide logical array-based access to the data for this
  //        data variable.  This is the const operator version.
  //
  // \tparam E A complex index type.
  //
  // This version of the operator is provided to support use with
  // \e flecsi mesh entity types \ref mesh_entity_base_t.
	///
  template<typename E>
  const T &
  operator [] (
    E * e
  ) const
  {
    return this->operator[](e->template id<0>());
  } // operator []

	///
  // \brief Provide logical array-based access to the data for this
  //        data variable.  This is the const operator version.
  //
  // \tparam E A complex index type.
  //
  // This version of the operator is provided to support use with
  // \e flecsi mesh entity types \ref mesh_entity_base_t.
	///
  template<typename E>
  T &
  operator [] (
    E * e
  )
  {
    return this->operator[](e->template id<0>());
  } // operator []

	///
  // \brief Provide logical array-based access to the data for this
  //        data variable.  This is the const operator version.
  //
  // \param index The index of the data variable to return.
	///
  const T &
  operator [] (
    size_t index
  ) const
  {
    assert(index < size_ && "index out of range");
    return data_[index];
  } // operator []

	///
  // \brief Provide logical array-based access to the data for this
  //        data variable.  This is the const operator version.
  //
  // \param index The index of the data variable to return.
	///
  T &
  operator [] (
    size_t index
  )
  {
    assert(index < size_ && "index out of range");
    return data_[index];
  } // operator []

	///
  // \brief Test to see if this accessor is empty
  //
  // \return true if registered.
  ///
  operator bool() const
  {
    return data_ != nullptr;
  } // operator bool

private:

  std::string label_ = "";
  size_t size_ = 0;
  T * data_ = nullptr;
  const user_meta_data_t & user_meta_data_ =
    *(std::make_unique<user_meta_data_t>());
  index_space_t is_;
  bitset_t & user_attributes_ =
    *(std::make_unique<bitset_t>());
  size_t index_space_ = 0;

}; // struct dense_accessor_t

//----------------------------------------------------------------------------//
// Dense handle.
//----------------------------------------------------------------------------//

template<typename T>
struct dense_handle_t : public data_handle_t
{
  using type = T;
}; // struct dense_handle_t

//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//
// Main type definition.
//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//

//----------------------------------------------------------------------------//
// Dense storage type.
//----------------------------------------------------------------------------//

///
// FIXME: Dense storage type.
///
template<typename DS, typename MD>
struct storage_type_t<dense, DS, MD>
{
  //--------------------------------------------------------------------------//
  // Type definitions.
  //--------------------------------------------------------------------------//

  using data_store_t = DS;
  using meta_data_t = MD;

  template<typename T>
  using accessor_t = dense_accessor_t<T, MD>;

  template<typename T>
  using handle_t = dense_handle_t<T>;

  //--------------------------------------------------------------------------//
  // Data registration.
  //--------------------------------------------------------------------------//

  ///
  // \tparam T Data type to register.
  // \tparam NS Namespace
  // \tparam Args Variadic arguments that are passed to
  //              metadata initialization.
  //
  // \param data_client Base class reference to client.
  // \param data_store A reference for accessing the low-level data.
  // \param key A const string instance containing the variable name.
  // \param versions The number of variable versions for this datum.
  // \param indices The number of indices in the index space.
  ///
  template<
    typename T,
    size_t NS,
    typename ... Args
  >
  static
  handle_t<T>
  register_data(
    data_client_t & data_client,
    data_store_t & data_store,
    const const_string_t & key,
    size_t versions,
    size_t index_space,
    Args && ... args
  )
  {
    size_t h = key.hash() ^ data_client.runtime_id();

    // Runtime assertion that this key is unique
    assert(data_store[NS].find(h) == data_store[NS].end() &&
      "key already exists");

    //------------------------------------------------------------------------//
    // Call the user meta data initialization method passing variadic
    // user arguments
    //------------------------------------------------------------------------//

    data_store[NS][h].user_data.initialize(std::forward<Args>(args) ...);

    //------------------------------------------------------------------------//
    // Set the data label
    //------------------------------------------------------------------------//
 
    data_store[NS][h].label = key.c_str();

    //------------------------------------------------------------------------//
    // Set the data size by calling the data clients indeces method.
    // This allows the user to interpret the index space argument
    // in whatever way they want.
    //------------------------------------------------------------------------//
 
    data_store[NS][h].size = data_client.indices(index_space);

    //------------------------------------------------------------------------//
    // Store the index space.
    //------------------------------------------------------------------------//

    data_store[NS][h].index_space = index_space;

    //------------------------------------------------------------------------//
    // Store the data type size information.
    //------------------------------------------------------------------------//

    data_store[NS][h].type_size = sizeof(T);

    //------------------------------------------------------------------------//
    // This allows us to set the runtime-type-information, which requires
    // a const reference.
    //------------------------------------------------------------------------//

    data_store[NS][h].rtti.reset(
      new typename meta_data_t::type_info_t(typeid(T)));

    //------------------------------------------------------------------------//
    // Store the number of versions.
    //------------------------------------------------------------------------//

    data_store[NS][h].versions = versions;

    //------------------------------------------------------------------------//
    // Allocate data for each version.
    //------------------------------------------------------------------------//

    for(size_t i=0; i<versions; ++i) {
      data_store[NS][h].attributes[i].reset();
      data_store[NS][h].data[i].resize(
        data_client.indices(index_space) * sizeof(T));
    } // for

    //------------------------------------------------------------------------//
    // num_entries is unused for this storage type.
    //------------------------------------------------------------------------//

    data_store[NS][h].num_entries = 0;

    return {};
  } // register_data

  //--------------------------------------------------------------------------//
  // Data accessors.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  template<
    typename T,
    size_t NS
  >
  static
  accessor_t<T>
  get_accessor(
    data_client_t & data_client,
    data_store_t & data_store,
    const const_string_t & key,
    size_t version
  )
  {
    const size_t h = key.hash() ^ data_client.runtime_id();
    auto search = data_store[NS].find(h);

    if(search == data_store[NS].end()) {
      return {};
    }
    else {
      auto & meta_data = search->second;

      // check that the requested version exists
      assert(meta_data.versions > version && "version out-of-range");

      return { meta_data.label, meta_data.size,
        reinterpret_cast<T *>(&meta_data.data[version][0]),
				meta_data.user_data, meta_data.attributes[version],
        meta_data.index_space };
    } // if
  } // get_accessor

  ///
  //
  ///
  template<
    typename T,
    size_t NS,
    typename P
  >
  static
  std::vector<accessor_t<T>>
  get_accessors(
    data_client_t & data_client,
    P && preficate
  )
  {
    return {};
  } // get_accessors

  //--------------------------------------------------------------------------//
  // Data handles.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  template<
    typename T,
    size_t NS
  >
  static
  handle_t<T>
  get_handle(
    data_client_t & data_client,
    data_store_t & data_store,
    const const_string_t & key,
    size_t version
  )
  {
    return {};
  } // get_handle

}; // struct storage_type_t

} // namespace serial
} // namespace data
} // namespace flecsi

#endif // flecsi_serial_dense_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
