// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file bfl_tree_iterators.hpp
 * 
 * This library provides detail-implementations of BFL-tree iterators.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2013
 */

#ifndef BOOST_BFL_TREE_ITERATORS_HPP
#define BOOST_BFL_TREE_ITERATORS_HPP

#include <boost/config.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <boost/variant.hpp>

#include <boost/graph/detail/boost_container_generators.hpp>

#include <iterator>

namespace boost {

namespace graph { 

namespace detail {




struct bfltree_edge_desc {
  std::size_t source_vertex;
  std::size_t edge_index;
  bfltree_edge_desc(std::size_t aSrc = 0, std::size_t aEdgeId = 0) : 
                    source_vertex(aSrc), edge_index(aEdgeId) { };
  
};

inline bool operator==( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) { 
  return ((lhs.source_vertex == rhs.source_vertex) && (lhs.edge_index == rhs.edge_index)); 
};
inline bool operator!=( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) { 
  return ((lhs.source_vertex != rhs.source_vertex) || (lhs.edge_index != rhs.edge_index)); 
};
inline bool operator <( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) {
  if( lhs.source_vertex == rhs.source_vertex )
    return ( lhs.edge_index < rhs.edge_index );
  else 
    return ( lhs.source_vertex < rhs.source_vertex );
};
inline bool operator<=( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) {
  if( lhs.source_vertex == rhs.source_vertex )
    return ( lhs.edge_index <= rhs.edge_index );
  else 
    return ( lhs.source_vertex < rhs.source_vertex );
};
inline bool operator >( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) {
  if( lhs.source_vertex == rhs.source_vertex )
    return ( lhs.edge_index > rhs.edge_index );
  else 
    return ( lhs.source_vertex > rhs.source_vertex );
};
inline bool operator>=( const bfltree_edge_desc& lhs, const bfltree_edge_desc& rhs) {
  if( lhs.source_vertex == rhs.source_vertex )
    return ( lhs.edge_index >= rhs.edge_index );
  else 
    return ( lhs.source_vertex > rhs.source_vertex );
};




template <typename VProp>
bool bfltree_is_vertex_valid(const VProp& vp) {
  return (vp.out_degree != std::numeric_limits<std::size_t>::max());
};

template <typename Container>
struct bfltree_vertex_validity {
  const Container* p_cont;
  explicit bfltree_vertex_validity(const Container* aPCont = NULL) : p_cont(aPCont) { };
  bool operator()(std::size_t d) {
    return ((d < p_cont->size()) && (bfltree_is_vertex_valid((*p_cont)[d])));
  };
};

template <typename Container, std::size_t Arity>
struct bfltree_edge_validity {
  const Container* p_cont;
  explicit bfltree_edge_validity(const Container* aPCont = NULL) : p_cont(aPCont) { };
  bool operator()(bfltree_edge_desc d) {
    std::size_t d_child = Arity * d.source_vertex + 1 + d.edge_index;
    return ((d.source_vertex < p_cont->size()) && bfltree_is_vertex_valid((*p_cont)[d.source_vertex])) &&
           ((d_child < p_cont->size()) && bfltree_is_vertex_valid((*p_cont)[d_child]));
  };
};





struct bfltree_oeiter : 
  public iterator_facade<
    bfltree_oeiter,
    const bfltree_edge_desc,
    std::random_access_iterator_tag
  > {
  public:
    
    typedef bfltree_oeiter self;
    
    explicit bfltree_oeiter(bfltree_edge_desc aE = bfltree_edge_desc()) : e(aE) { };
    
  public: // private:
    friend class iterator_core_access;
    
    void increment() { ++e.edge_index; };
    void decrement() { --e.edge_index; };
    bool equal(const self& rhs) const { return (this->e == rhs.e); };
    const bfltree_edge_desc& dereference() const { return e; };
    
    void advance(std::ptrdiff_t i) { e.edge_index += i; };
    std::ptrdiff_t distance_to(const self& rhs) const { return rhs.e.edge_index - this->e.edge_index; }; 
    
    bfltree_edge_desc e;
};


struct bfltree_ieiter : 
  public iterator_facade<
    bfltree_ieiter,
    const bfltree_edge_desc,
    std::random_access_iterator_tag
  > {
  public:
    
    typedef bfltree_ieiter self;
    
    explicit bfltree_ieiter(bfltree_edge_desc aE = bfltree_edge_desc()) : e(aE) { };
    
  public: // private:
    friend class iterator_core_access;
    
    void increment() { ++e.edge_index; };
    void decrement() { --e.edge_index; };
    bool equal(const self& rhs) const { return (this->e == rhs.e); };
    const bfltree_edge_desc& dereference() const { return e; };
    
    void advance(std::ptrdiff_t i) { e.edge_index += i; };
    std::ptrdiff_t distance_to(const self& rhs) const { return rhs.e.edge_index - this->e.edge_index; }; 
    
    bfltree_edge_desc e;
};






template <std::size_t Arity>
struct bfltree_eiter : 
  public iterator_facade<
    bfltree_eiter<Arity>,
    const bfltree_edge_desc,
    std::bidirectional_iterator_tag
  > {
  public:
    
    typedef bfltree_eiter<Arity> self;
    
    explicit bfltree_eiter(bfltree_edge_desc aE = bfltree_edge_desc()) : e(aE) { };
    
  public: // private:
    friend class iterator_core_access;
    
    void increment() { 
      ++e.edge_index; 
      if(e.edge_index == Arity) {
        ++e.source_vertex;
        e.edge_index = 0;
      };
    };
    void decrement() { 
      if(e.edge_index == 0) {
        --e.source_vertex;
        e.edge_index = Arity;
      };
      --e.edge_index; 
    };
    bool equal(const self& rhs) const { return (this->e == rhs.e); };
    const bfltree_edge_desc& dereference() const { return e; };
    
    bfltree_edge_desc e;
};



struct bfltree_viter : 
  public iterator_facade<
    bfltree_viter,
    const std::size_t,
    std::random_access_iterator_tag
  > {
  public:
    
    typedef bfltree_viter self;
    
    explicit bfltree_viter(std::size_t aVIt = 0) : v_it(aVIt) { };
    
  public: // private:
    friend class iterator_core_access;
    
    void increment() { ++v_it; };
    void decrement() { --v_it; };
    bool equal(const self& rhs) const { return (this->v_it == rhs.v_it); };
    const std::size_t& dereference() const { return v_it; };
    
    void advance(std::ptrdiff_t i) { v_it += i; };
    std::ptrdiff_t distance_to(const self& rhs) const { return rhs.v_it - this->v_it; }; 
    
    std::size_t v_it;
};




};

};

};


#endif





