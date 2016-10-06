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

#ifndef flecsi_serial_sparse_h
#define flecsi_serial_sparse_h

#include <map>
#include <algorithm>
#include <set>

//----------------------------------------------------------------------------//
// POLICY_NAMESPACE must be defined before including storage_type.h!!!
// Using this approach allows us to have only one storage_type_t
// definintion that can be used by all data policies -> code reuse...
#define POLICY_NAMESPACE serial
#include "flecsi/data/storage_type.h"
#undef POLICY_NAMESPACE
//----------------------------------------------------------------------------//

#include "flecsi/utils/const_string.h"

///
// \file serial/sparse.h
// \authors bergen
// \date Initial file creation: Apr 17, 2016
///

namespace flecsi {
namespace data {
namespace serial {

//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//
// Helper type definitions.
//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//

//----------------------------------------------------------------------------//
// Sparse accessor.
//----------------------------------------------------------------------------//

using index_pair_ = std::pair<size_t, size_t>;

///
//
///
template<typename T>
struct material_value__
{
  material_value__(size_t material)
  : material(material){}

  material_value__(size_t material, T value)
  : material(material),
  value(value){}

  material_value__(){}

  size_t material;
  T value;
}; // struct material_value__

static constexpr size_t INDICES_FLAG = 1UL << 63;
static constexpr size_t MATERIALS_FLAG = 1UL << 62;

///
// \brief sparse_accessor_t provides logically array-based access to data
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
struct sparse_accessor_t
{
  //--------------------------------------------------------------------------//
  // Type definitions.
  //--------------------------------------------------------------------------//

  using iterator_t = index_space_t::iterator_t;
  using meta_data_t = MD;
  using user_meta_data_t = typename meta_data_t::user_meta_data_t;

  using material_value_t = material_value__<T>;

  //--------------------------------------------------------------------------//
  // Constructors.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  sparse_accessor_t() {}

  ///
  //
  ///
  sparse_accessor_t(
    meta_data_t & meta_data,
    size_t version
  )
  :
  meta_data_(meta_data),
  version_(version),
  num_indices_(meta_data.size),
  num_materials_(meta_data_.num_materials)
  {
    auto iitr = meta_data_.data.find(version | INDICES_FLAG);
    assert(iitr != meta_data_.data.end());

    std::vector<uint8_t> & raw_indices = 
      const_cast<std::vector<uint8_t> &>(iitr->second);

    indices_ = reinterpret_cast<size_t *>(&raw_indices[0]);

    auto mitr = meta_data_.data.find(version | MATERIALS_FLAG);
    assert(mitr != meta_data_.data.end());

    std::vector<uint8_t> & raw_materials = 
      const_cast<std::vector<uint8_t> &>(mitr->second);

    materials_ = reinterpret_cast<material_value_t *>(&raw_materials[0]);      
  } // sparse_accessor_t

  //--------------------------------------------------------------------------//
  // Member data interface.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  bitset_t &
  attributes()
  {
    return meta_data_.user_attributes_;
  } // attributes

  //--------------------------------------------------------------------------//
  // Operators.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  T &
  operator () (
    size_t index,
    size_t material
  )
  {
    assert(index < num_indices_ && "sparse accessor: index out of bounds");

    material_value_t * start = materials_ + indices_[index];
    material_value_t * end = materials_ + indices_[index + 1];

    material_value_t * itr = 
      std::lower_bound(start, end, material_value_t(material),
        [](const auto & k1, const auto & k2) -> bool {
          return k1.material < k2.material;
        });

    assert(itr != end && "sparse accessor: unmapped material");

    return itr->value;
  } // operator ()

  ///
  //
  ///
  void dump()
  {
    for(size_t i = 0; i < num_indices_; ++i) {
      size_t start = indices_[i];
      size_t end = indices_[i + 1];

      std::cout << "+++++ row: " << i << std::endl;

      for(size_t j = start; j < end; ++j) {
        material_value_t & mj = materials_[j];
        std::cout << "++ material: " << mj.material << std::endl;
        std::cout << "++ value: " << mj.value << std::endl << std::endl;
      } // for
    } // for
  } // dump

  ///
  //
  ///
  T *
  data()
  {
    return nullptr;
  } // data

private:

  size_t version_;
  const meta_data_t & meta_data_ =
    *(std::make_unique<meta_data_t>());
  size_t num_indices_;
  size_t num_materials_;
  size_t * indices_;
  material_value_t * materials_;

}; // struct sparse_accessor_t

template<typename T, typename MD>
struct sparse_mutator_t {

  //--------------------------------------------------------------------------//
  // Type definitions.
  //--------------------------------------------------------------------------//
  
  using iterator_t = index_space_t::iterator_t;
  using meta_data_t = MD;
  using user_meta_data_t = typename meta_data_t::user_meta_data_t;

  using material_value_t = material_value__<T>;

  //--------------------------------------------------------------------------//
  // Constructors.
  //--------------------------------------------------------------------------//

  ///
  //
  ///
  sparse_mutator_t() {}

  ///
  //
  ///
  sparse_mutator_t(
    size_t num_slots,
    const std::string & label,
    size_t version,
    meta_data_t & meta_data,
    const user_meta_data_t & user_meta_data
  )
  :
    num_slots_(num_slots),
    label_(label),
    version_(version), 
    meta_data_(meta_data),
    num_indices_(meta_data_.size),
    num_materials_(meta_data_.num_materials),
    indices_(new index_pair_[num_indices_]),
    materials_(new material_value_t[num_indices_ * num_slots_]),
    erase_set_(nullptr)
  {}

  ///
  //
  ///
  ~sparse_mutator_t()
  {
    commit();
  } // ~sparse_mutator_t

  T &
  operator () (
    size_t index,
    size_t material
  )
  {
    assert(indices_ && "sparse mutator has alread been committed");
    assert(index < num_indices_ && material < num_materials_);

    index_pair_ & ip = indices_[index];
    
    size_t n = ip.second - ip.first;
    
    if(n >= num_slots_) {
      return spare_map_.emplace(index,
        material_value_t(material))->second.value;
    } // if

    material_value_t * start = materials_ + index * num_slots_;     
    material_value_t * end = start + n;

    material_value_t * itr = 
      std::lower_bound(start, end, material_value_t(material),
        [](const auto & k1, const auto & k2) -> bool{
          return k1.material < k2.material;
        });

    while(end != itr) {
      *(end) = *(end - 1);
      --end;
    } // while

    itr->material = material;

    ++ip.second;

    return itr->value;
  } // operator ()

  void
  erase(
    size_t index,
    size_t material
  )
  {
    assert(indices_ && "sparse mutator has alread been committed");
    assert(index < num_indices_ && material < num_materials_);

    if(!erase_set_){
      erase_set_ = new erase_set_t;
    }

    erase_set_->emplace(std::make_pair(index, material));
  }

  T *
  data()
  {
    return nullptr;
  } // data

  void
  commit()
  {
    if(!indices_) {
      return;
    } // if

    auto iitr = meta_data_.data.find(version_ | INDICES_FLAG);
    assert(iitr != meta_data_.data.end());

    std::vector<uint8_t> & raw_indices = 
      const_cast<std::vector<uint8_t> &>(iitr->second);

    size_t * indices = 
      reinterpret_cast<size_t *>(&raw_indices[0]);

    auto mitr = meta_data_.data.find(version_ | MATERIALS_FLAG);
    assert(mitr != meta_data_.data.end());

    std::vector<uint8_t> & raw_materials = 
      const_cast<std::vector<uint8_t> &>(mitr->second);

    constexpr size_t ev_bytes = sizeof(material_value_t);

    size_t s = raw_materials.size();
    size_t c = raw_materials.capacity();
    size_t d = (num_indices_ * num_slots_ + spare_map_.size()) * ev_bytes;

    if(c - s < d) {
      //raw_materials.reserve(c * ev_bytes + d);
      raw_materials.resize(c * ev_bytes + d);
    } // if

    material_value_t * materials = 
      reinterpret_cast<material_value_t *>(&raw_materials[0]);

    material_value_t * materials_end = materials + s / ev_bytes;

    for(size_t i = 0; i < num_indices_; ++i) {
      index_pair_ & ip = indices_[i];
      
      size_t n = ip.second - ip.first;

      size_t pos = indices[i];

      if(n == 0) {
        indices[i + 1] = pos;
        continue;
      } // if

      auto p = spare_map_.equal_range(i);
      
      size_t m = distance(p.first, p.second);

      material_value_t * start = materials_ + i * num_slots_;     
      material_value_t * end = start + n; 

      auto iitr = materials + pos;

      std::copy(iitr, materials_end, iitr + n + m);
      std::copy(start, end, iitr);
      materials_end += n + m;

      auto cmp = [](const auto & k1, const auto & k2) -> bool {
        return k1.material < k2.material;
      };

      auto nitr = iitr + indices[i + 1];

      std::inplace_merge(iitr, nitr, nitr + n, cmp);

      if(m == 0) {
        indices[i + 1] = pos + n;
        continue;
      } // if

      indices[i + 1] = pos + n + m;

      auto vitr = nitr + n;

      for(auto itr = p.first; itr != p.second; ++itr) {
        vitr->material = itr->second.material;
        vitr->value = itr->second.value;
        ++vitr;
      } // for

      std::inplace_merge(iitr, nitr + n, nitr + n + m, cmp);
    } // for

    delete[] indices_;
    indices_ = nullptr;

    delete[] materials_;
    materials_ = nullptr;

    if(!erase_set_){
      return;
    }

    constexpr size_t erased_marker = std::numeric_limits<size_t>::max();

    material_value_t * start; 
    material_value_t * end;

    size_t last = std::numeric_limits<size_t>::max();
  
    size_t num_erased = 0;

    for(auto p : *erase_set_){
      if(p.first != last){
        start = materials + indices[p.first];
        end = materials + indices[p.first + 1];
        last = p.first;
        indices[p.first] -= num_erased;      
      }

      material_value_t * m = 
        std::lower_bound(start, end, material_value_t(p.second),
          [](const auto & k1, const auto & k2) -> bool {
            return k1.material < k2.material;
        });

      m->material = erased_marker;
      start = m + 1;
      ++num_erased;
    }

    while(last < num_indices_){
      indices[++last] -= num_erased;
    }

    std::remove_if(materials, materials_end,
      [](const auto & k) -> bool {
        return k.material == erased_marker;
    });

    s = raw_materials.size();
    s -= num_erased * ev_bytes;
    raw_materials.resize(s);

    delete erase_set_;
    erase_set_ = nullptr;
  } // commit

private:
  using spare_map_t = std::multimap<size_t, material_value__<T>>;
  using erase_set_t = std::set<std::pair<size_t, size_t>>;

  size_t num_slots_;
  std::string label_ = "";
  size_t version_ = 0;
  const meta_data_t & meta_data_ =
    *(std::make_unique<meta_data_t>());
    

  size_t num_indices_;
  size_t num_materials_;
  index_pair_ * indices_;
  material_value__<T> * materials_;
  spare_map_t spare_map_;
  erase_set_t * erase_set_;
}; // struct sparse_accessor_t

//----------------------------------------------------------------------------//
// Sparse handle.
//----------------------------------------------------------------------------//

template<typename T>
struct sparse_handle_t {
}; // struct sparse_handle_t

//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//
// Main type definition.
//+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=//

///
// FIXME: Sparse storage type.
///
template<typename DS, typename MD>
struct storage_type_t<sparse, DS, MD> {

  //--------------------------------------------------------------------------//
  // Type definitions.
  //--------------------------------------------------------------------------//

  using data_store_t = DS;
  using meta_data_t = MD;

  template<typename T>
  using accessor_t = sparse_accessor_t<T, MD>;

  template<typename T>
  using mutator_t = sparse_mutator_t<T, MD>;

  template<typename T>
  using handle_t = sparse_handle_t<T>;

  //--------------------------------------------------------------------------//
  // Data registration.
  //--------------------------------------------------------------------------//

  template<typename T, size_t NS, typename ... Args>
  static handle_t<T> register_data(data_client_t & data_client,
    data_store_t & data_store, const const_string_t & key,
    size_t versions, size_t indices, size_t num_materials, Args && ... args) {

    size_t h = key.hash() ^ data_client.runtime_id();

    // Runtime assertion that this key is unique
    assert(data_store[NS].find(h) == data_store[NS].end() &&
      "key already exists");

    meta_data_t & md = data_store[NS][h];

    md.user_data.initialize(std::forward<Args>(args) ...);
    md.label = key.c_str();
    md.size = indices;
    md.type_size = sizeof(T);
    md.versions = versions;
    md.rtti.reset(
      new typename meta_data_t::type_info_t(typeid(T)));

    md.num_materials = num_materials;

    for(size_t version = 0; version < versions; ++version){
      auto & iv = md.data[version | INDICES_FLAG] = std::vector<uint8_t>();
      iv.resize((indices + 1) * sizeof(size_t));

      md.data[version | MATERIALS_FLAG] = std::vector<uint8_t>();
    }

    return {};
  } // register_data

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

      return { meta_data, version };
    } // if
  } // get_accessor

  ///
  //
  ///
  template<
    typename T,
    size_t NS
  >
  static
  mutator_t<T>
  get_mutator(
    data_client_t & data_client,
    data_store_t & data_store,
    const const_string_t & key,
    size_t slots,
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

      return { slots, meta_data.label, version, meta_data,
        meta_data.user_data };
    } // if
  } // get_accessor

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

#endif // flecsi_serial_sparse_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
