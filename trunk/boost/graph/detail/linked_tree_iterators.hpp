// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file linked_tree_iterators.hpp
 * 
 * This library provides detail-implementations of (doubly-)linked tree iterators.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_LINKED_TREE_ITERATORS_HPP
#define BOOST_LINKED_TREE_ITERATORS_HPP

#include <boost/config.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include <boost/variant.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/int.hpp>

#include <iterator>

namespace boost {


namespace graph { 

namespace detail {
  
  
  template <typename ContIterator>
  class ltree_viter_from_iter : 
    public iterator_facade<
      ltree_viter_from_iter<ContIterator>,
      const ContIterator,
      typename std::iterator_traits<ContIterator>::iterator_category
    > {
    public:
      typedef ltree_viter_from_iter<ContIterator> self;
      
      explicit ltree_viter_from_iter(ContIterator aVIt = ContIterator()) : v_it(aVIt) { };
      
      /* Utility: */
      template <typename VertexContainer>
      static self begin(VertexContainer& c) { return self(c.begin()); };
      template <typename VertexContainer>
      static self end(VertexContainer& c) { return self(c.end()); };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++v_it; };
      void decrement() { --v_it; };
      bool equal(const self& rhs) const { return this->v_it == rhs.v_it; };
      const ContIterator& dereference() const { return v_it; };
      
      ContIterator v_it;
  };
  
  
  struct ltree_viter_from_index : 
    public iterator_facade<
      ltree_viter_from_index,
      const std::size_t,
      std::random_access_iterator_tag
    > {
    public:
      
      typedef ltree_viter_from_index self;
      
      explicit ltree_viter_from_index(std::size_t aVIt = 0) : v_it(aVIt) { };
      
      /* Utility: */
      template <typename VertexContainer>
      static self begin(VertexContainer&) { return self(0); };
      template <typename VertexContainer>
      static self end(VertexContainer& c) { return self(c.size()); };
      
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
  
  
  template <typename VDesc, typename ContIterator>
  class ltree_child_viter_from_iter : 
    public iterator_facade<
      ltree_child_viter_from_iter<VDesc, ContIterator>,
      const VDesc,
      typename std::iterator_traits<ContIterator>::iterator_category
    > {
    public:
      typedef ltree_child_viter_from_iter<VDesc, ContIterator> self;
      
      explicit ltree_child_viter_from_iter(ContIterator aEdgeIt = ContIterator()) : e_it(aEdgeIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++e_it; };
      void decrement() { --e_it; };
      bool equal(const self& rhs) const { return (this->e_it == rhs.e_it); };
      const VDesc& dereference() const { return e_it->target; };
      
      void advance(std::ptrdiff_t i) { std::advance(this->v_it, i); };
      std::ptrdiff_t distance_to(const self& rhs) const { return std::distance(this->v_it, rhs.v_it); }; 
      
      ContIterator e_it;
  };
  
  
  
  template <typename EDesc>
  class ltree_oeiter_from_iter : 
    public iterator_facade<
      ltree_oeiter_from_iter<EDesc>,
      const EDesc,
      std::bidirectional_iterator_tag
    > {
    public:
      typedef ltree_oeiter_from_iter<EDesc> self;
      
      explicit ltree_oeiter_from_iter(EDesc aEdge = EDesc()) : e_it(aEdge) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++(e_it.edge_id); };
      void decrement() { --(e_it.edge_id); };
      bool equal(const self& rhs) const { return (this->e_it.source == rhs.e_it.source) && (this->e_it.edge_id == rhs.e_it.edge_id); };
      const EDesc& dereference() const { return e_it; };
      
      EDesc e_it;
  };
  
  
  template <typename EDesc>
  class ltree_oeiter_from_index : 
    public iterator_facade<
      ltree_oeiter_from_index<EDesc>,
      const EDesc,
      std::random_access_iterator_tag
    > {
    public:
      typedef ltree_oeiter_from_index<EDesc> self;
      
      explicit ltree_oeiter_from_index(EDesc aEdge = EDesc()) : e_it(aEdge) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++(e_it.edge_id); };
      void decrement() { --(e_it.edge_id); };
      bool equal(const self& rhs) const { return (this->e_it.source == rhs.e_it.source) && (this->e_it.edge_id == rhs.e_it.edge_id); };
      const EDesc& dereference() const { return e_it; };
      
      void advance(std::ptrdiff_t i) { this->e_it.edge_id += i; };
      std::ptrdiff_t distance_to(const self& rhs) const { return rhs.e_it.edge_id - this->e_it.edge_id; }; 
      
      EDesc e_it;
  };
  
  
  template <typename EDesc, typename ContIterator>
  class ltree_eiter_from_iter : 
    public iterator_facade<
      ltree_eiter_from_iter<EDesc, ContIterator>,
      const EDesc,
      std::bidirectional_iterator_tag
    > {
    public:
      typedef ltree_eiter_from_iter<EDesc, ContIterator> self;
      
      explicit ltree_eiter_from_iter(ContIterator aVIt = ContIterator()) : v_it(aVIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++v_it; };
      void decrement() { --v_it; };
      bool equal(const self& rhs) const { return (this->v_it == rhs.v_it); };
      const EDesc& dereference() const { return v_it->in_edge; };
      
      ContIterator v_it;
  };
  
  template <typename EDesc, typename ContIterator>
  class ltree_eiter_from_pooliter : 
    public iterator_facade<
      ltree_eiter_from_pooliter<EDesc, ContIterator>,
      const EDesc,
      std::random_access_iterator_tag
    > {
    public:
      typedef ltree_eiter_from_pooliter<EDesc, ContIterator> self;
      
      explicit ltree_eiter_from_pooliter(ContIterator aVIt = ContIterator()) : v_it(aVIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++v_it; };
      void decrement() { --v_it; };
      bool equal(const self& rhs) const { return (this->v_it == rhs.v_it); };
      const EDesc& dereference() const { 
        typedef typename std::iterator_traits<ContIterator>::value_type var_type;
        typedef typename mpl::at< typename var_type::types, mpl::int_<0> >::type value_type;
        if(v_it->which() == 0)
          return get<value_type>(*v_it).in_edge;
        else
          return EDesc::null_value;
      };
      
      void advance(std::ptrdiff_t i) { this->v_it += i; };
      std::ptrdiff_t distance_to(const self& rhs) const { return rhs.v_it - this->v_it; }; 
      
      ContIterator v_it;
  };
  
  
  
}; 


};



};


#endif


















