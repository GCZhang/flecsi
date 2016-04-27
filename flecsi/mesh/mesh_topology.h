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

#ifndef flecsi_mesh_topology_h
#define flecsi_mesh_topology_h

/*!
  \file mesh_topology.h
  \authors nickm@lanl.gov, bergen@lanl.gov
  \date Initial file creation: Sep 23, 2015
 */

#include <algorithm>
#include <iostream>
#include <array>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <functional>
#include <map>

#include "flecsi/utils/common.h"
#include "flecsi/utils/set_intersection.h"
#include "flecsi/mesh/mesh_types.h"

namespace flecsi
{
/*----------------------------------------------------------------------------*
 * class mesh_topology_t
 *----------------------------------------------------------------------------*/

/*!
  \class mesh_topology_t mesh_topology.h
  \brief mesh_topology_t is parameterized on a class (MT) which gives
    information about its entity types, connectivities and more. the mesh
    topology is responsibly for computing connectivity info between entities
    of different topological dimension, e.g: vertex -> cell,
    cell -> edge, etc. and provides methods for traversing these adjancies.
    It also holds vectors containing the entity instances.
 */
template <class MT>
class mesh_topology_t : public mesh_topology_base_t
{
 public:

  // used to find the entity type of topological dimension D and domain M
  template <size_t D, size_t M = 0>
  using entity_type = typename find_entity_<MT, D, M>::type;

  /*--------------------------------------------------------------------------*
   * class iterator
   *--------------------------------------------------------------------------*/

  /*!
   \class iterator mesh_topology.h
   \brief An iterator that returns entities for topological dimension D
    and domain M.
  */
  template <size_t D, size_t M = 0>
  class iterator
  {
   public:
    using mesh_t = mesh_topology_t;
    using entity_type = typename find_entity_<MT, D, M>::type;

    // construct a top-level iterator, e.g: cells of a mesh
    iterator(const iterator & itr)
        : mesh_(itr.mesh_), entities_(itr.entities_), index_(itr.index_)
    {
    }

    // construct a nested iterator, e.g: edges of a cell
    iterator(mesh_t & mesh, const id_vector_t & entities, size_t index)
        : mesh_(mesh), entities_(&entities), index_(index)
    {
    }

    iterator & operator++()
    {
      ++index_;
      return *this;
    } // operator ++

    iterator & operator=(const iterator & itr)
    {
      index_ = itr.index_;
      entities_ = itr.entities_;
      return *this;
    } // operator =

    domain_entity<M, entity_type> operator*()
    {
      return mesh_.get_entity<D, M>((*entities_)[index_]);
    } // operator *

    // allow the entity methods to be called into
    entity_type * operator->()
    {
      return mesh_.get_entity<D, M>((*entities_)[index_]);
    } // operator ->

    bool operator==(const iterator & itr) const
    {
      return index_ == itr.index_;
    } // operator ==

    bool operator!=(const iterator & itr) const
    {
      return index_ != itr.index_;
    } // operator !=

   private:
    mesh_t & mesh_;
    const id_vector_t * entities_;
    size_t index_;

  }; // class iterator

  /*--------------------------------------------------------------------------*
   * class const_iterator
   *--------------------------------------------------------------------------*/

  /*!
   \class const_iterator mesh_topology.h
   \brief A const iterator that returns entities for topological dimension D
   and domain M.
  */
  template <size_t D, size_t M = 0>
  class const_iterator
  {
   public:
    using mesh_t = const mesh_topology_t;
    using entity_type = typename find_entity_<MT, D, M>::type;

    // construct a top-level iterator, e.g: cells of a mesh
    const_iterator(const const_iterator & itr)
        : mesh_(itr.mesh_), entities_(itr.entities_), index_(itr.index_)
    {
    }

    // construct a nested iterator, e.g: edges of a cell
    const_iterator(mesh_t & mesh, const id_vector_t & entities,
        size_t index)
        : mesh_(mesh), entities_(&entities), index_(index)
    {
    }

    const_iterator & operator++()
    {
      ++index_;
      return *this;
    } // operator ++

    const_iterator & operator=(const const_iterator & itr)
    {
      index_ = itr.index_;
      entities_ = itr.entities_;
      return *this;
    } // operator =

    domain_entity<M, entity_type> operator*()
    {
      return mesh_.get_entity<D, M>((*entities_)[index_]);
    } // operator *

    const entity_type * operator->() const
    {
      return mesh_.get_entity<D, M>((*entities_)[index_]);
    } // operator ->

    bool operator==(const const_iterator & itr) const
    {
      return index_ == itr.index_;
    } // operator ==

    bool operator!=(const const_iterator & itr) const
    {
      return index_ != itr.index_;
    } // operator !=

   private:
    mesh_t & mesh_;
    const id_vector_t * entities_;
    size_t index_;

  }; // class iterator

  /*--------------------------------------------------------------------------*
   * class entity_set
   *--------------------------------------------------------------------------*/

  template <size_t D, class I, size_t M = 0>
  class entity_set
  {
  public:
  using iterator_t = I;
  using mesh_t = typename iterator_t::mesh_t;
  using entity_type = typename iterator_t::entity_type;
  using domain_entity_t = domain_entity<M, entity_type>;
  using domain_entity_vector_t = std::vector<domain_entity_t>;
  
  using filter_function = std::function<bool(domain_entity_t)>;
  
  using apply_function = std::function<void(domain_entity_t)>;

  template<typename T>
  using map_function = std::function<T(domain_entity_t)>;

  template<typename T>
  using reduce_function = std::function<void(domain_entity_t,T&)>;

  // default constructor
  entity_set() = default;

  // Top-level constructor, e.g: cells of a mesh
  entity_set(mesh_t & mesh, const id_vector_t & v,
            bool sorted = false)
  : mesh_(&mesh), v_(&v), begin_(0), end_(v_->size()),
  owned_(false), sorted_(sorted) { }

     // Nested constructor, e.g: edges of a cell
  entity_set(
            mesh_t & mesh, const id_vector_t & v, size_t begin,
            size_t end, bool sorted = false)
  : mesh_(&mesh), v_(&v), begin_(begin), end_(end),
  owned_(false), sorted_(sorted) { }

  entity_set(const entity_set & r)
  : mesh_(r.mesh_), begin_(0), owned_(r.owned_), sorted_(r.sorted_)
  {
    if (owned_) {
      v_ = new id_vector_t(*r.v_);
    }
    else {
      v_ = r.v_;
    }

    end_ = v_->size();
  }

  // Top-level constructor, e.g: cells of a mesh
  entity_set(mesh_t & mesh, id_vector_t && v, bool sorted)
  : mesh_(&mesh), v_(new id_vector_t(std::move(v))),
  begin_(0), end_(v_->size()), owned_(true), sorted_(sorted) { }

  ~entity_set(){
    if(owned_){
      delete v_;
    }
  }

  entity_set & operator=(const entity_set & r) = default;

  iterator_t begin() const { return iterator_t(*mesh_, *v_, begin_); } // begin
  
  iterator_t end() const { return iterator_t(*mesh_, *v_, end_); } // end
   
  /*!
   convert this range to a vector
  */
  domain_entity_vector_t to_vec() const
  {
    domain_entity_vector_t ret;
    for (size_t i = begin_; i < end_; ++i) {
      ret.push_back(mesh_->template get_entity<D, M>((*v_)[i]));
    } // for
       
    return ret;
  } // to_vec
     
  domain_entity<M, entity_type> operator[](size_t i) const
  {
    return mesh_->template get_entity<D, M>((*v_)[begin_ + i]);
  } // []
 
  domain_entity<M, entity_type> at(size_t i) const
  {
    assert( i >= begin_ && i < end_ );
    return mesh_->get_entity<D, M>((*v_)[begin_ + i]);
  } // at
     
   
  domain_entity<M, entity_type> front() const
  {
    return mesh_->template get_entity<D, M>((*v_)[begin_]);
  } // first
   
  domain_entity<M, entity_type> back() const
  {
    return mesh_->template get_entity<D, M>((*v_)[end_ - 1]);
  } // last
   
  size_t size() const { return end_ - begin_; } // size
   
  entity_set filter(filter_function f) const {
    id_vector_t v;

    for (auto ent : *this) {
      if (f(ent)) {
        v.push_back(ent.id());
      }
    }

    return entity_set(*mesh_, std::move(v), sorted_);
  }

  template<typename T>
  std::vector<entity_set> scatter(map_function<T> f) const {

    std::unordered_map<T, id_vector_t> id_map;
    for (auto ent : *this)
      id_map[f(ent)].push_back(ent.id());

    std::vector<entity_set> ent_map;
    ent_map.reserve( id_map.size() );
    for ( auto entry : id_map )
      ent_map.emplace_back( 
        std::move( entity_set(*mesh_, std::move(entry.second), sorted_) )
      );

    return ent_map;
  }
   
  void apply(apply_function f) const {
    for (auto ent : *this) {
      f(ent);
    }
  }
   
  template<typename T>
  std::vector<T> map(map_function<T> f) const {
    std::vector<T> ret;
    ret.reserve(v_->size());
    
    for (auto ent : *this) {
      ret.push_back(f(ent));
    }
    return ret;
  }

  template<typename T>
  T reduce(T start, reduce_function<T> f) const {
    T r = start;
    
    for (auto ent : *this) {
      f(ent, r);
    }

    return r;
  }
   
  void prepare_(){
    if(!owned_){
      v_ = new id_vector_t(*v_);
      owned_ = true;
    }

    if(!sorted_){
      auto vc = const_cast<id_vector_t*>(v_);
      std::sort(vc->begin(), vc->end());
      sorted_ = true;
    }
  }
   
  entity_set& operator&=(const entity_set& r){
    prepare_();

    id_vector_t ret;

    if(r.sorted_){
      ret.resize(std::min(v_->size(), r.v_->size()));

      auto itr = std::set_intersection(v_->begin(), v_->end(),
                                       r.v_->begin(), r.v_->end(), ret.begin());

      ret.resize(itr - ret.begin());
    }
    else{
      id_vector_t v2(*r.v_);
      std::sort(v2.begin(), v2.end());

      ret.resize(std::min(v_->size(), v2.size()));

      auto itr = std::set_intersection(v_->begin(), v_->end(),
                                       v2.begin(), v2.end(), ret.begin());

      ret.resize(itr - ret.begin());
    }

    delete v_;
    v_ = new id_vector_t(std::move(ret));

    begin_ = 0;
    end_ = v_->size();

    return *this;
  }
   
  entity_set operator&(const entity_set& r) const{
    entity_set ret(*this);
    ret &= r;
    return ret;
  }
   
  entity_set& operator|=(const entity_set& r){
    prepare_();

    id_vector_t ret;

    if(r.sorted_){
      ret.resize(v_->size() + r.v_->size());

      auto itr = std::set_union(v_->begin(), v_->end(),
                                r.v_->begin(), r.v_->end(), ret.begin());

     ret.resize(itr - ret.begin());
    }
    else{
      id_vector_t v2(*r.v_);

      std::sort(v2.begin(), v2.end());

      ret.resize(v_->size() + v2.size());

      auto itr = std::set_union(v_->begin(), v_->end(),
                                v2.begin(), v2.end(), ret.begin());

      ret.resize(itr - ret.begin());
    }

    delete v_;
    v_ = new id_vector_t(std::move(ret));

    begin_ = 0;
    end_ = v_->size();

    return *this;
  }
   
  entity_set operator|(const entity_set& r) const{
    entity_set ret(*this);
    ret |= r;
    return ret;
  } 
   
  entity_set& operator-=(const entity_set& r){
    prepare_();

    id_vector_t ret(v_->size());

    if(r.sorted_){
      auto itr = std::set_difference(v_->begin(), v_->end(),
                                     r.v_->begin(), r.v_->end(), ret.begin());

      ret.resize(itr - ret.begin());
    }
    else{
      id_vector_t v2(*r.v_);

      std::sort(v2.begin(), v2.end());

      auto itr = std::set_difference(v_->begin(), v_->end(),
                                     v2.begin(), v2.end(), ret.begin());

      ret.resize(itr - ret.begin());
    }

    delete v_;
    v_ = new id_vector_t(std::move(ret));

    begin_ = 0;
    end_ = v_->size();

    return *this;
  }
   
  entity_set operator-(const entity_set& r) const{
    entity_set ret(*this);
    ret -= r;
    return ret;
  }
   
  void add(const domain_entity<M, entity_type>& ent){
    if(!owned_){
      v_ = new id_vector_t(*v_);
      owned_ = true;
    }

    auto vc = const_cast<id_vector_t*>(v_);

    if(sorted_){
      auto id = ent.id();
      auto itr = std::upper_bound(vc->begin(), vc->end(), id);
      vc->insert(itr, id);
    }
    else{
      vc->push_back(ent.id());
    }
  }
   
  entity_set& operator<<(const domain_entity<M, entity_type>& ent){
    add(ent);
    return *this;
  }
     
  private:
    mesh_t * mesh_ = nullptr;
    const id_vector_t * v_ = nullptr;
    size_t begin_ = 0;
    size_t end_ = 0;
    bool owned_ = false;
    bool sorted_ = true;
  }; // class entity_set

  template<size_t D, size_t M = 0>
  using entity_set_t = entity_set<D, iterator<D, M>, M>;

  template<size_t D, size_t M = 0>
  using const_entity_set_t = entity_set<D, const_iterator<D, M>, M>;

  /*--------------------------------------------------------------------------*
   * class iterator
   *--------------------------------------------------------------------------*/

  /*!
   \class id_iterator mesh_topology.h
   \brief An iterator that returns only the id's of entities
    (not the entities themselves) for performance reasons.
  */
  class id_iterator
  {
   public:
    id_iterator(const id_iterator & itr)
        : entities_(itr.entities_), index_(itr.index_)
    {
    }

    id_iterator(const id_vector_t & entities, size_t index)
        : entities_(&entities), index_(index)
    {
    }

    id_iterator & operator++()
    {
      ++index_;
      return *this;
    } // operator ++

    id_iterator & operator=(const id_iterator & itr)
    {
      index_ = itr.index_;
      entities_ = itr.entities_;
      return *this;
    } // oerator =

    id_t operator*() { return (*entities_)[index_]; }
    bool operator==(const id_iterator & itr) const
    {
      return index_ == itr.index_;
    } // operator ==

    bool operator!=(const id_iterator & itr) const
    {
      return index_ != itr.index_;
    } // operator !=

   private:
    const id_vector_t * entities_;
    size_t index_;

  }; // class iterator

  /*--------------------------------------------------------------------------*
   * class id_range
   *--------------------------------------------------------------------------*/

  /*!
    \class id_range mesh_topology.h
    \brief Used to implement range-based for iteration for id iterators.
   */

  class id_range
  {
   public:
    id_range(const id_vector_t & v) : v_(v), begin_(0), end_(v_.size()) {}
    id_range(const id_vector_t & v, id_t begin, id_t end)
        : v_(v), begin_(begin), end_(end)
    {
    }

    id_range(const id_range & r) : v_(r.v_), begin_(0), end_(v_.size()) {}
    id_iterator begin() const { return id_iterator(v_, begin_); }
    id_iterator end() const { return id_iterator(v_, end_); }
    // Convert this range into a vector which can then be indexed
    id_vector_t to_vec() const
    {
      id_vector_t ret;

      for (size_t i = begin_; i < end_; ++i) {
        ret.push_back(v_[i]);
      } // for

      return ret;
    } // to_vec

    auto operator[](size_t i) const
    { return v_[i]; }

    size_t size() const { return end_ - begin_; } // size
   private:
    const id_vector_t & v_;
    id_t begin_;
    id_t end_;

  }; // class id_range

  // Don't allow the mesh to be copied or copy constructed

  mesh_topology_t(const mesh_topology_t &) = delete;

  mesh_topology_t & operator=(const mesh_topology_t &) = delete;

  // Allow move operations
  mesh_topology_t(mesh_topology_t &&) = default;

  //! override default move assignement
  mesh_topology_t & operator=(mesh_topology_t && o) = default;

  //! Constructor
  mesh_topology_t()
  {
    // initialize all lower connectivities because the user might 
    // specify different combinations of connections
    for (size_t i = 1; i < MT::dimension+1; ++i)
      for (size_t j = 0; j < i; ++j)
        get_connectivity_(0, i, j).init();
  } // mesh_topology_t()

  // The mesh retains ownership of the entities and deletes them
  // upon mesh destruction
  virtual ~mesh_topology_t()
  {
    for (size_t d = 0; d < MT::num_domains; ++d) {
      for (auto & ev : ms_.entities[d]) {
        for (auto ent : ev) {
          delete ent;
        }
      }
    }
  }

  // Add and entity to a mesh domain and assign its id per domain
  template <size_t D, size_t M = 0>
  void add_entity(mesh_entity_base_t<MT::num_domains> * ent,
                  size_t partition_id=0)
  {
    auto & ents = ms_.entities[M][D];

    id_t global_id = id_t::make<D, M>(ents.size(), partition_id);

    ent->ids_[M] = global_id;
    ents.push_back(ent);

    auto & id_vec = ms_.id_vecs[M][D];
    id_vec.push_back(global_id);
  } // add_entity

  // A mesh is constructed by creating cells and vertices and associating
  // vertices with cells as in this method.
  template <size_t M, class C, typename V>
  void init_cell(C * cell, V && verts)
  {
    init_cell_<M>(cell, std::forward<V>(verts) );
  } // init_cell

  template <size_t M, class C, typename V>
  void init_cell(C * cell, std::initializer_list<V *> verts)
  {
    init_cell_<M>(cell, verts );
  } // init_cell

  template < size_t M, typename V >
  void init_cell_(entity_type<MT::dimension, M> * cell, V && verts)
  {
    auto & c = get_connectivity_(M, MT::dimension, 0);

    assert(cell->template id<M>() == c.from_size() && "id mismatch");

    for (entity_type<0, M> * v : std::forward<V>(verts) ) {
      c.push(v->template global_id<M>());
    } // for

    c.end_from();
  } // init_cell


  // Initialize an entities connectivity with a subset of another
  template < size_t M, size_t D1, size_t D2, class E1, class E2 >
  void init_entity(E1 * super, E2 && subs)
  {
    init_entity_<M,D1,D2>( super, std::forward<E2>(subs) );
  } // init_entity

  template < size_t M, size_t D1, size_t D2, class E1, class E2 >
  void init_entity(E1 * super, std::initializer_list<E2*> subs )
  {
    init_entity_<M,D1,D2>( super, subs );
  } // init_entity

  template < size_t M, size_t D1, size_t D2, class E2 >
  void init_entity_(
    entity_type<D1, M> * super, E2 && subs)
  {
    auto & c = get_connectivity_(M, D1, D2);

    assert(super->template id<M>() == c.from_size() && "id mismatch");

    for ( auto e : subs ) {
      c.push( e->template global_id<M>() );
    } // for

    c.end_from();
  } // init_entity

  // Get the number of entities in a given domain and topological dimension
  size_t num_entities_(size_t dim, size_t domain=0) const
  {
    return ms_.entities[domain][dim].size();
  } // num_entities_

  // Virtual method of num_entities_()
  size_t num_entities(size_t dim, size_t domain=0) const override
  {
    return num_entities_(dim, domain);
  } // num_entities

  /*!
    Build connectivity informaiton and add entities to the mesh for the
    given dimension.
   */
  template <size_t M, size_t D>
  void build_connectivity()
  {
    // std::cerr << "build: " << D << std::endl;

    // Sanity check
    assert(D <= MT::dimension);

    // Reference to storage from cells to the entity (to be created here).
    connectivity_t & cell_to_entity = get_connectivity_(M, MT::dimension, D);

    // Storage for entity-to-vertex connectivity information.
    connection_vector_t entity_vertex_conn;

    // Helper variables
    size_t entity_id = 0;
    size_t max_cell_entity_conns = 1;

    // Get connectivity for cells to vertices.
    connectivity_t & cell_to_vertex = get_connectivity_(M, MT::dimension, 0);
    assert(!cell_to_vertex.empty());

    const size_t _num_cells = num_entities<MT::dimension, M>();

    // Storage for cell-to-entity connectivity information.
    connection_vector_t cell_entity_conn(_num_cells);

    // This map is primarily used to make sure that entities are not
    // created multiple times, i.e., that they are unique.  The
    // emplace method of the map is used to only define a new entity
    // if it does not already exist in the map.
    id_vector_map_t entity_vertices_map;

    // This buffer should be large enough to hold all entities
    // vertices that potentially need to be created
    std::array<id_t, 4096> entity_vertices;

    for (size_t c = 0; c < _num_cells; ++c) {
      // Get the cell object
      auto cell = static_cast<entity_type<MT::dimension, M> *>(
          ms_.entities[M][MT::dimension][c]);
      
      id_t cell_id = cell->template global_id<M>();

      // Get storage reference.
      id_vector_t & conns = cell_entity_conn[c];

      // Try to optimize storage.
      conns.reserve(max_cell_entity_conns);

      // Get the vertices that define the current cell
      size_t end_index;
      id_t * vertices = cell_to_vertex.get_entities(c, end_index);

      // This call allows the users specialization to create
      // whatever entities are needed to complete the mesh.
      //
      // p.first:   The number of entities per cell.
      // p.second:  A std::vector of id_t containing the ids of the
      //            vertices that define the entity.
      auto sv =
          cell->create_entities(D, entity_vertices.data(), vertices, end_index);

      size_t n = sv.size();

      // iterate over the newly-defined entities
      for (size_t i = 0; i < n; ++i) {
        size_t m = sv[i];

        // Get the vertices that define this entity by getting
        // a pointer to the vector-of-vector data and then constructing
        // a vector of ids for only this entity.
        id_t * a = &entity_vertices[i * m];
        id_vector_t ev(a, a + m);

        // Sort the ids for the current entity so that they are
        // monotonically increasing. This ensures that entities are
        // created uniquely (using emplace_back below) because the ids
        // will always occur in the same order for the same entity.
        std::sort(ev.begin(), ev.end());

        // Emplace the sorted vertices into the entity map
        auto itr = entity_vertices_map.emplace(
            std::move(ev), id_t::make<D, M>(entity_id, cell_id.partition()));

        // Add this id to the cell to entity connections
        conns.push_back(itr.first->second);

        // If the insertion took place
        if (itr.second) {
          // what does this do?
          id_vector_t ev2 = id_vector_t(a, a + m);
          entity_vertex_conn.emplace_back(std::move(ev2));

          max_cell_entity_conns = 
            std::max(max_cell_entity_conns, conns.size());

          // A new entity was added, so we advance the id counter.
          ++entity_id;
        } // if
      } // for
    } // for

    // This call will create the entity objects in the mesh (The above
    // logic only defines the indices and connectivity.)
    cell_to_entity.init_create<MT, M, D>(
        ms_.id_vecs[M][D], ms_.entities[M][D], cell_entity_conn, *this);

    // Set the connectivity information from the created entities to
    // the vertices.
    connectivity_t & entity_to_vertex = get_connectivity_(M, D, 0);
    entity_to_vertex.init(entity_vertex_conn);
  } // build_connectivity

  /*!
     used internally to compute connectivity information for
     topological dimension
       FD -> TD where FD < TD
   */
  template <size_t FM, size_t TM, size_t FD, size_t TD>
  void transpose()
  {
    // std::cerr << "transpose: " << FD << " -> " << TD << std::endl;

    // The connectivity we will be populating
    auto & out_conn = get_connectivity_(FM, TM, FD, TD);
    if (!out_conn.empty()) {
      return;
    } // if

    index_vector_t pos(num_entities_(FD, FM), 0);

    for (auto to_entity : entities<TD, TM>()) {
      for (id_t from_id : entity_ids<FD, TM, FM>(to_entity)) {
        ++pos[from_id.entity()];
      }
    }

    out_conn.resize(pos);

    std::fill(pos.begin(), pos.end(), 0);

    for (auto to_entity : entities<TD, TM>()) {
      for (id_t from_id : entity_ids<FD, TM, FM>(to_entity)) {
        out_conn.set(from_id.entity(), to_entity->template global_id<TM>(),
            pos[from_id.entity()]++);
      }
    }
  } // transpose

  /*!
     Used internally to compute connectivity information for
     topological dimension
       FD -> TD using FD -> D' and D' -> TD
   */
  template <size_t FM, size_t TM, size_t FD, size_t TD, size_t D>
  void intersect()
  {
    // std::cerr << "intersect: " << FD << " -> " << TD << std::endl;

    // The connectivity we will be populating
    connectivity_t & out_conn = get_connectivity_(FM, TM, FD, TD);
    if (!out_conn.empty()) {
      return;
    } // if

    // the number of each entity type
    auto num_from_ent = num_entities_(FD, FM);
    auto num_to_ent = num_entities_(TD, FM);

    // Temporary storage for connection id's
    connection_vector_t conns(num_from_ent);

    // Keep track of which to id's we have visited
    using visited_vec = std::vector<bool>;
    visited_vec visited(num_to_ent);

    size_t max_size = 1;

    // Read connectivities
    connectivity_t & c = get_connectivity_(FM, FD, D);
    assert(!c.empty());

    connectivity_t & c2 = get_connectivity_(TM, TD, D);
    assert(!c2.empty());

    // Iterate through entities in from topological dimension
    for (auto from_entity : entities<FD, FM>()) {
      id_t from_id = from_entity->template global_id<FM>();
      id_vector_t & ents = conns[from_id.entity()];
      ents.reserve(max_size);

      size_t count;
      id_t * ep = c.get_entities(from_id.entity(), count);

      // Create a copy of to vertices so they can be sorted
      id_vector_t from_verts(ep, ep+count);
      // sort so we have a unique key for from vertices
      std::sort(from_verts.begin(), from_verts.end());

      // initially set all to id's to unvisited
      for (auto from_ent2 : entities<D, FM>(from_entity)) {
        for (id_t to_id : entity_ids<TD, TM>(from_ent2)) {
          visited[to_id.entity()] = false;
        }
      }

      // Loop through each from entity again
      for (auto from_ent2 : entities<D, FM>(from_entity)) {
        for (id_t to_id : entity_ids<TD, TM>(from_ent2)) {

          // If we have already visited, skip
          if (visited[to_id.entity()]) {
            continue;
          } // if

          visited[to_id.entity()] = true;

          // If the topological dimensions are the same, always add to id
          if (FD == TD) {
            if (from_id != to_id) {
              ents.push_back(to_id);
            } // if
          } else {
            size_t count;
            id_t * ep = c2.get_entities(to_id.entity(), count);

            // Create a copy of to vertices so they can be sorted
            id_vector_t to_verts(ep, ep + count);
            // Sort to verts so we can do an inclusion check
            std::sort(to_verts.begin(), to_verts.end());

            // If from vertices contains the to vertices add to id
            // to this connection set
            if ( D < TD ) {
              if ( std::includes( from_verts.begin(), from_verts.end(),
                                  to_verts.begin(), to_verts.end()) )
                ents.emplace_back(to_id); 
            } 
            // If we are going through a higher level, then set
            // intersection is sufficient. i.e. one set does not need to 
            // be a subset of the other
            else {
              if ( utils::intersects( from_verts.begin(), from_verts.end(),
                                      to_verts.begin(), to_verts.end() ) )
                ents.emplace_back(to_id); 
            } // if

          } // if
        } // for
      } // for

      max_size = std::max(ents.size(), max_size);
    } // for

    // Finally create the connection from the temporary conns
    out_conn.init(conns);
  } // intersect

  /*!
     Used to compute connectivity information for topological dimension
       D1 -> D2
   */
  template <size_t M, size_t FD, size_t TD>
  void compute_connectivity()
  {
    // std::cerr << "compute: " << FD << " -> " << TD << std::endl;

    // Get the output connectivity
    connectivity_t & out_conn = get_connectivity_(M, FD, TD);

    // Check if we have already computed it
    if (!out_conn.empty()) {
      return;
    } // if

    // if we don't have cell -> vertex connectivities, then
    // try building cell -> vertex connectivity through the
    // faces (3d) or edges(2d)
    static_assert( MT::dimension <= 3, 
                   "this needs to be re-thought for higher dimensions" );

    if ( get_connectivity_(M, MT::dimension, 0).empty() ) {
      assert( !get_connectivity_(M, MT::dimension-1, 0).empty() && 
              " need at least edges(2d)/faces(3) -> vertex connectivity" );
      // assume we have cell -> faces, so invert it to get faces -> cells
      transpose<M, M, MT::dimension-1, MT::dimension>();
      // invert faces -> vertices to get vertices -> faces
      transpose<M, M, 0, MT::dimension-1>();
      // build cells -> vertices via intersections with faces
      intersect<M, M, MT::dimension, 0, MT::dimension-1>();
    }

    // Check if we need to build entities, e.g: edges or faces
    if (num_entities_(FD, M) == 0) {
      build_connectivity<M, FD>();
    } // if

    if (num_entities_(TD, M) == 0) {
      build_connectivity<M, TD>();
    } // if

    if (num_entities_(FD, M) == 0 && num_entities_(TD, M) == 0) {
      return;
    } // if

    // Depending on the corresponding topological dimensions, call transpose
    // or intersect as need
     if (FD < TD) {
      compute_connectivity<M, TD, FD>();
      transpose<M, M, FD, TD>();
    } else {
       if (FD == 0 && TD == 0) {
         // compute vertex to vertex connectivities through shared cells.
         compute_connectivity<M, FD, MT::dimension>();
         compute_connectivity<M, MT::dimension, TD>();
         intersect<M, M, FD, TD, MT::dimension>();
       } else {
         // computer connectivities through shared vertices.
         compute_connectivity<M, FD, 0>();
         compute_connectivity<M, 0, TD>();
         intersect<M, M, FD, TD, 0>();
       }
    } // if
  } // compute_connectivity

  template <size_t FM, size_t TM, size_t FD, size_t TD>
  void compute_bindings()
  {
    connectivity_t & out_conn = get_connectivity_(FM, TM, FD, TD);

    if (!out_conn.empty()) {
      return;
    } // if

    if (FD <= TD) {
      connectivity_t & trans_conn = get_connectivity_(TM, FM, TD, FD);
      if(!trans_conn.empty() && 
         !ms_.entities[FD].empty() &&
         !ms_.entities[TD].empty()){
        transpose<FM, TM, FD, TD>();
        return;        
      }
    }

    if (FD == TD) {
      connectivity_t & out_conn = get_connectivity_(FM, TM, FD, TD);
    }

    if (ms_.entities[TM][TD].empty()) {
      build_bindings<FM, TM, TD>();
    }

  } // compute_bindings

  /*!
    Build bindings associated with a from/to domain and topological dimension.
    compute_bindings will call this on each binding found in the tuple of
    bindings specified in the mesh type/traits mesh specialization.
   */
  template <size_t FM, size_t TM, size_t TD>
  void build_bindings()
  {
    // Sanity check
    static_assert(TD <= MT::dimension, "invalid dimension");

    // Helper variables
    size_t entity_id = 0;
    size_t max_cell_conns = 1;
    const size_t _num_cells = num_entities<MT::dimension, FM>();

    // Storage for cell connectivity information
    connection_vector_t cell_conn(_num_cells);

    // Get cell definitions from domain 0
    using ent_vec_t = entity_vector_t<MT::num_domains>;
    ent_vec_t & cells = ms_.entities[FM][MT::dimension];

    static constexpr size_t M0 = 0;

    for (size_t i = 0; i < MT::dimension; ++i) {
      get_connectivity_<TM, FM, TD>(i).init();
    }

    std::array<id_t *, MT::dimension> primal_ids;
    std::array<size_t, MT::dimension> num_primal_ids;

    // This buffer should be large enough to hold all entities
    // that potentially need to be created
    std::array<id_t, 4096> entity_ids;

    // Iterate over cells
    for (auto c : cells) {
      // Map used to ensure unique entity creation
      id_vector_map_t entity_ids_map;

      // Get a cell object.
      auto cell = static_cast<entity_type<MT::dimension, M0> *>(c);
      id_t cell_id = cell->template global_id<FM>();

      // Get ids of entities with at least this dimension
      for (size_t dim = 0; dim < MT::dimension; ++dim) {
        auto & c = get_connectivity_<FM, FM, MT::dimension>(dim);
        primal_ids[dim] = c.get_entities( cell_id.entity(), num_primal_ids[dim] );
      } // for

      // p.first:   The number of entities per cell.
      // p.second:  A std::vector of id_t containing the ids of the
      //            entities that define the bound entity.

      auto sv = cell->create_bound_entities(
        FM, TM, TD, primal_ids.data(), num_primal_ids.data(), entity_ids.data() );

      size_t n = sv.size();

      // Iterate over the newly-defined entities
      id_vector_t & conns = cell_conn[cell_id.entity()];

      conns.reserve(n);

      size_t pos = 0;

      for (size_t i = 0; i < n; ++i) {
        size_t m = sv[i];

        id_t create_id = id_t::make<TD, TM>(entity_id, cell_id.partition());

        // Add this id to the cell entity connections
        conns.push_back(create_id);

        uint32_t dim_flags = 0;

        for (size_t k = 0; k < m; ++k) {
          id_t global_id = entity_ids[pos + k];

          get_connectivity_<TM, FM, TD>(global_id.dimension()).push(global_id);
          dim_flags |= 1U << global_id.dimension();
        }

        for (size_t i = 0; i < MT::dimension; ++i) {
          if (dim_flags & (1U << i)) {
            get_connectivity_<TM, FM, TD>(i).end_from();
          }
        }

        ++entity_id;

        pos += m;
      } // for
    } // for

    // Reference to storage from cells to the entity (to be created here).
    connectivity_t & cell_out = get_connectivity_(FM, TM, MT::dimension, TD);

    if (ms_.entities[TM][TD].empty()) {
      // Create the entity objects
      cell_out.init_create<MT, TM, TD>(
          ms_.id_vecs[TM][TD], ms_.entities[TM][TD], cell_conn, *this);
    } else {
      cell_out.init(cell_conn);
    }

  } // build_bindings

  /*!
    The init method builds entities as edges/faces and computes adjacencies
    and bindings.
   */
  template <size_t M = 0>
  void init()
  {
    // Compute mesh connectivity
    using TP = typename MT::connectivities;
    compute_connectivity_<M, std::tuple_size<TP>::value, TP>::compute(*this);

    using BT = typename MT::bindings;
    compute_bindings_<M, std::tuple_size<BT>::value, BT>::compute(*this);
  } // init

  /*!
    Similar to init(), but only compute bindings. This method should be called
    when a domain is sparse, i.e: missing certain entity types such as cells
    and it is not possible to compute connectivities.
   */
  template <size_t M = 0>
  void init_bindings()
  {
    using BT = typename MT::bindings;
    compute_bindings_<M, std::tuple_size<BT>::value, BT>::compute(*this);
  } // init

  /*!
   Return the number of entities contained in specified topological dimension
   and domain.
   */
  template <size_t D, size_t M = 0>
  decltype(auto) num_entities() const
  {
    return ms_.entities[M][D].size();
  } // num_entities

  /*!
   Get the connectivity of the specified from/to domain and from/to topological
   dimensions.
   */
  const connectivity_t & get_connectivity(size_t from_domain, size_t to_domain,
      size_t from_dim, size_t to_dim) const override
  {
    return get_connectivity_(from_domain, to_domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified from/to domain and from/to topological
   dimensions.
   */
  connectivity_t & get_connectivity(size_t from_domain, size_t to_domain,
      size_t from_dim, size_t to_dim) override
  {
    return get_connectivity_(from_domain, to_domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified domain and from/to topological
   dimensions.
   */
  const connectivity_t & get_connectivity(
      size_t domain, size_t from_dim, size_t to_dim) const override
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Get the connectivity of the specified domain and from/to topological
   dimensions.
   */
  connectivity_t & get_connectivity(
      size_t domain, size_t from_dim, size_t to_dim) override
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  const connectivity_t & get_connectivity_(size_t from_domain, size_t to_domain,
      size_t from_dim, size_t to_dim) const
  {
    assert(from_dim < ms_.topology[from_domain][to_domain].size() &&
        "invalid from_dim");
    auto & t = ms_.topology[from_domain][to_domain][from_dim];
    assert(to_dim < t.size() && "invalid to_dim");
    return t[to_dim];
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  connectivity_t & get_connectivity_(
      size_t from_domain, size_t to_domain, size_t from_dim, size_t to_dim)
  {
    assert(from_dim < ms_.topology[from_domain][to_domain].size() &&
        "invalid from_dim");
    auto & t = ms_.topology[from_domain][to_domain][from_dim];
    assert(to_dim < t.size() && "invalid to_dim");
    return t[to_dim];
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  template <size_t FM, size_t TM, size_t FD>
  connectivity_t & get_connectivity_(size_t to_dim)
  {
    assert(FD < ms_.topology[FM][TM].size() && "invalid from_dim");
    auto & t = ms_.topology[FM][TM][FD];
    assert(to_dim < t.size() && "invalid to_dim");
    return t[to_dim];
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  template <size_t FM, size_t TM, size_t FD, size_t TD>
  connectivity_t & get_connectivity_()
  {
    return ms_.topology[FM][TM][FD][TD];
  } // get_connectivity

  /*!
   Implementation of get_connectivity for various get_connectivity convenience
   methods.
   */
  const connectivity_t & get_connectivity_(
      size_t domain, size_t from_dim, size_t to_dim) const
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  connectivity_t & get_connectivity_(
      size_t domain, size_t from_dim, size_t to_dim)
  {
    return get_connectivity_(domain, domain, from_dim, to_dim);
  } // get_connectivity

  size_t topological_dimension() const override { return MT::dimension; }
  /*!
    This method should be called to construct and entity rather than
    calling the constructor directly. This way, the ability to have
    extra initialization behavior is reserved.
  */
  template <class T, class... S>
  T * make(S &&... args)
  {
    T * entity = new T(std::forward<S>(args)...);
    return entity;
  } // make

  template <size_t M = 0>
  const entity_vector_t<MT::num_domains> & get_entities_(size_t dim) const
  {
    return ms_.entities[M][dim];
  } // get_entities_

  template <size_t M = 0>
  const id_vector_t & get_id_vec_(size_t dim) const
  {
    return ms_.id_vecs[M][dim];
  } // get_id_vec_

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template <size_t D, size_t M = 0>
  auto get_entity(id_t global_id) const
  {
    using entity_type = typename find_entity_<MT, D, M>::type;
    return static_cast<entity_type *>(
        ms_.entities[M][D][global_id.entity()]);
  } // get_entity

  /*!
    Get an entity in domain M of topological dimension D with specified id.
  */
  template <size_t M = 0>
  auto get_entity(size_t dim, id_t global_id)
  {
    return ms_.entities[M][dim][global_id.entity()];
  } // get_entity

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM, size_t TM = FM, class E>
  const_entity_set_t<D, TM> entities(const E * e) const
  {
    const connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    const index_vector_t & fv = c.get_from_index_vec();
    return const_entity_set_t<D, TM>(*this, c.get_entities(),
        fv[e->template id<FM>()], fv[e->template id<FM>() + 1]);
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM, size_t TM = FM, class E>
  entity_set_t<D, TM> entities(E * e)
  {
    const connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    const index_vector_t & fv = c.get_from_index_vec();
    return entity_set_t<D, TM>(*this, c.get_entities(),
        fv[e->template id<FM>()], fv[e->template id<FM>() + 1]);
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM = 0, size_t TM = FM, class E>
  decltype(auto) entities(domain_entity<FM, E> & e) const
  {
    return entities<D, FM, TM>(e.entity()());
  } // entities

  /*!
    Get the entities of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM = 0, size_t TM = FM, class E>
  decltype(auto) entities(domain_entity<FM, E> & e)
  {
    return entities<D, FM, TM>(e.entity());
  } // entities

  /*!
    Get the top-level entities of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template <size_t D, size_t M = 0>
  const_entity_set_t<D, M> entities() const
  {
    assert(!ms_.id_vecs[M][D].empty());
    return const_entity_set_t<D>(*this, ms_.id_vecs[M][D], true);
  } // entities

  /*!
    Get the top-level entities of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template <size_t D, size_t M = 0>
  entity_set_t<D, M> entities()
  {
    assert(!ms_.id_vecs[M][D].empty());
    return entity_set_t<D, M>(*this, ms_.id_vecs[M][D], true);
  } // entities

  /*!
    Get the top-level entity id's of topological dimension D of the specified
    domain M. e.g: cells of the mesh.
  */
  template <size_t D, size_t M = 0>
  id_range entity_ids() const
  {
    assert(!ms_.id_vecs[M][D].empty());
    return id_range(ms_.id_vecs[M][D]);
  } // entity_ids

  /*!
    Get the entity id's of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM = 0, size_t TM = FM, class E>
  decltype(auto) entity_ids(domain_entity<FM, E> & e)
  {
    return entity_ids<D, FM, TM>(e.entity());
  } // entities

  /*!
    Get the entity id's of topological dimension D connected to another entity
    by specified connectivity from domain FM and to domain TM.
  */
  template <size_t D, size_t FM = 0, size_t TM = FM, class E>
  id_range entity_ids(const E * e) const
  {
    const connectivity_t & c = get_connectivity(FM, TM, E::dimension, D);
    assert(!c.empty() && "empty connectivity");
    const index_vector_t & fv = c.get_from_index_vec();
    return id_range(c.get_entities(), fv[e->template id<FM>()],
        fv[e->template id<FM>() + 1]);
  } // entities

  template<typename I>
  void compute_graph_partition(
    size_t domain,
    size_t dim,
    const std::vector<I>& partition_sizes,
    std::vector<mesh_graph_partition<I>>& partitions){

    using int_t = I;

    partitions.reserve(partition_sizes.size());

    int_t total_size = 0;
    for(auto pi : partition_sizes){
      total_size += pi;
    }

    size_t n = num_entities_(dim, domain);
    size_t pn = n / total_size;

    size_t to_dim;

    if (dim == 0) {
      // vertex -> vertex via shared edge.
      to_dim = 1;
    } else {
      // edge -> edge via shared vertex, cell -> cell via shared edge/face etc.
      to_dim = dim - 1;
    }

    const connectivity_t& c1 = get_connectivity(domain, dim, to_dim);
    assert(!c1.empty() && "empty connectivity c1");
    const index_vector_t& fv1 = c1.get_from_index_vec();

    const connectivity_t& c2 = get_connectivity(domain, to_dim, dim);
    assert(!c2.empty() && "empty connectivity c2");
    const index_vector_t& fv2 = c2.get_from_index_vec();

    mesh_graph_partition<int_t> cp;
    cp.offset.reserve(pn);

    size_t offset = 0;
    size_t pi = 0;

    std::vector<int_t> partition;
    partition.push_back(0);

    for(size_t from_id = 0; from_id < n; ++from_id){
      auto to_ids = id_range(c1.get_entities(), fv1[from_id], fv1[from_id + 1]);
      cp.offset.push_back(offset);
      
      for(auto to_id : to_ids){
        auto ret_ids = id_range(c2.get_entities(), fv2[to_id], fv2[to_id + 1]);
        
        for(auto ret_id : ret_ids){
          if(ret_id != from_id){
            cp.index.push_back(ret_id);
            ++offset;
          }
        }
      }

      size_t m = cp.offset.size();

      if(m >= pn * partition_sizes[pi]){
        partitions.emplace_back(std::move(cp));
        partition.push_back(m + partition.back());
        offset = 0;
        ++pi;
      }
    }

    for(auto& pi : partitions){
      pi.partition = partition;
    }
  }

  /*!
    Debug method to dump the connectivity of the mesh over all domains and
    topological dimensions.
  */
  void dump()
  {
    for (size_t from_domain = 0; from_domain < MT::num_domains; ++from_domain) {
      std::cout << "=================== from domain: " << from_domain
                << std::endl;
      for (size_t to_domain = 0; to_domain < MT::num_domains; ++to_domain) {
        std::cout << "========== to domain: " << to_domain << std::endl;
        size_t n = ms_.topology[from_domain][to_domain].size();
        for (size_t i = 0; i < n; ++i) {
          auto & ci = ms_.topology[from_domain][to_domain][i];
          for (size_t j = 0; j < ci.size(); ++j) {
            auto & cj = ci[j];
            std::cout << "------------- " << i << " -> " << j << std::endl;
            cj.dump();
          }
        }
      }
    }
  } // dump

 private:
  mesh_storage_t<MT::dimension, MT::num_domains> ms_;

}; // class mesh_topology_t

} // flecsi

#endif // flecsi_mesh_topology_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
