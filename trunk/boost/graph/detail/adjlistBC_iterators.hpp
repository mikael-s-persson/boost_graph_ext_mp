// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file adjlistBC_iterators.hpp
 * 
 * This library provides detail-implementations of Boost.Container Adjacency-list iterators.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
 */

#ifndef BOOST_ADJLISTBC_ITERATORS_HPP
#define BOOST_ADJLISTBC_ITERATORS_HPP

#include <boost/config.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <boost/variant.hpp>

#include <boost/graph/detail/boost_container_generators.hpp>

#include <iterator>

namespace boost {

namespace graph { 

namespace detail {
  
  
  template <typename ContIterator>
  class adjlistBC_viter_from_iter : 
    public iterator_facade<
      adjlistBC_viter_from_iter<ContIterator>,
      const ContIterator,
      typename std::iterator_traits<ContIterator>::iterator_category
    > {
    public:
      typedef adjlistBC_viter_from_iter<ContIterator> self;
      
      explicit adjlistBC_viter_from_iter(ContIterator aVIt = ContIterator()) : v_it(aVIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++v_it; };
      void decrement() { --v_it; };
      bool equal(const self& rhs) const { return this->v_it == rhs.v_it; };
      const ContIterator& dereference() const { return v_it; };
      
      ContIterator v_it;
  };
  
  
  struct adjlistBC_viter_from_index : 
    public iterator_facade<
      adjlistBC_viter_from_index,
      const std::size_t,
      std::random_access_iterator_tag
    > {
    public:
      
      typedef adjlistBC_viter_from_index self;
      
      explicit adjlistBC_viter_from_index(std::size_t aVIt = 0) : v_it(aVIt) { };
      
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
  
  
  template <typename ListS, typename Container>
  struct adjlistBC_select_vertex_iterator {
    typedef adjlistBC_viter_from_iter< typename Container::iterator > type;
    
    static std::pair< type, type > create_range(Container& cont) {
      return std::pair< type, type >(type(cont.begin()), type(cont.end()));
    };
  };
  
  template <typename Container>
  struct adjlistBC_select_vertex_iterator<vecBC, Container> {
    typedef adjlistBC_viter_from_index type;
    
    static std::pair< type, type > create_range(Container& cont) {
      return std::pair< type, type >(type(0), type(cont.size()));
    };
  };
  
  template <typename Container>
  struct adjlistBC_select_vertex_iterator<poolBC, Container> {
    typedef typename Container::container_type RawContainer;
    typedef filter_iterator< BC_is_not_hole<RawContainer>, adjlistBC_viter_from_index > type;
    
    static std::pair< type, type > create_range(Container& cont) {
      adjlistBC_viter_from_index v_beg(0);
      adjlistBC_viter_from_index v_end(cont.m_data.size());
      return std::pair< type, type >(
        type( BC_is_not_hole<RawContainer>(&cont.m_data), v_beg, v_end ), 
        type( BC_is_not_hole<RawContainer>(&cont.m_data), v_end, v_end ));
    };
  };
  
  
  template <typename EDesc>
  class adjlistBC_oeiter_from_iter : 
    public iterator_facade<
      adjlistBC_oeiter_from_iter<EDesc>,
      const EDesc,
      std::bidirectional_iterator_tag
    > {
    public:
      typedef adjlistBC_oeiter_from_iter<EDesc> self;
      
      explicit adjlistBC_oeiter_from_iter(EDesc aEdge = EDesc()) : e_it(aEdge) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++(e_it.edge_id); };
      void decrement() { --(e_it.edge_id); };
      bool equal(const self& rhs) const { return (this->e_it.source == rhs.e_it.source) && (this->e_it.edge_id == rhs.e_it.edge_id); };
      const EDesc& dereference() const { return e_it; };
      
      EDesc e_it;
  };
  
  
  template <typename EDesc>
  class adjlistBC_oeiter_from_index : 
    public iterator_facade<
      adjlistBC_oeiter_from_index<EDesc>,
      const EDesc,
      std::random_access_iterator_tag
    > {
    public:
      typedef adjlistBC_oeiter_from_index<EDesc> self;
      
      explicit adjlistBC_oeiter_from_index(EDesc aEdge = EDesc()) : e_it(aEdge) { };
      
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
  
  
  
  template <typename ListS, typename Container, typename EDesc>
  struct adjlistBC_select_out_edge_iterator {
    typedef adjlistBC_oeiter_from_iter< EDesc > type;
    
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container& cont) {
      return std::pair< type, type >(type(EDesc(u,cont.begin())), type(EDesc(u,cont.end())));
    };
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container* cont) { return create_range(u, *cont); };
  };
  
  template <typename Container, typename EDesc>
  struct adjlistBC_select_out_edge_iterator<vecBC, Container, EDesc> {
    typedef adjlistBC_oeiter_from_index< EDesc > type;
    
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container& cont) {
      return std::pair< type, type >(type(EDesc(u,0)), type(EDesc(u,cont.size())));
    };
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container* cont) { return create_range(u, *cont); };
  };
  
  template <typename Container, typename EDesc>
  struct adjlistBC_select_out_edge_iterator<poolBC, Container, EDesc> {
    typedef typename Container::container_type RawContainer;
    typedef filter_iterator< BC_edge_is_not_hole<RawContainer>, adjlistBC_oeiter_from_index< EDesc > > type;
    
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container& cont) {
      adjlistBC_oeiter_from_index< EDesc > e_beg(EDesc(u,0));
      adjlistBC_oeiter_from_index< EDesc > e_end(EDesc(u,cont.m_data.size()));
      return std::pair< type, type >(
        type( BC_edge_is_not_hole<RawContainer>(&cont.m_data), e_beg, e_end ), 
        type( BC_edge_is_not_hole<RawContainer>(&cont.m_data), e_end, e_end ));
    };
    template <typename Vertex>
    static std::pair< type, type > create_range(Vertex u, Container* cont) { return create_range(u, *cont); };
  };
  
  
  
  
  template <typename VDesc, typename OEdgeContainer, typename OEdgeIterator>
  class adjlistBC_adjacent_viter : 
    public iterator_facade<
      adjlistBC_adjacent_viter<VDesc, OEdgeContainer, OEdgeIterator>,
      const VDesc,
      typename std::iterator_traits<OEdgeIterator>::iterator_category
    > {
    public:
      typedef adjlistBC_adjacent_viter<VDesc, OEdgeContainer, OEdgeIterator> self;
      
      explicit adjlistBC_adjacent_viter(const OEdgeContainer* aPCont = NULL, 
                                        OEdgeIterator aEdgeIt = OEdgeIterator()) : p_cont(aPCont), e_it(aEdgeIt) { };
      explicit adjlistBC_adjacent_viter(const OEdgeContainer& aPCont, 
                                        OEdgeIterator aEdgeIt) : p_cont(&aPCont), e_it(aEdgeIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++e_it; };
      void decrement() { --e_it; };
      bool equal(const self& rhs) const { return (this->e_it == rhs.e_it); };
      const VDesc& dereference() const { return BC_desc_to_iterator(*p_cont, e_it->edge_id)->target; };
      
      void advance(std::ptrdiff_t i) { std::advance(this->e_it, i); };
      std::ptrdiff_t distance_to(const self& rhs) const { return std::distance(this->e_it, rhs.e_it); }; 
      
      const OEdgeContainer* p_cont;
      OEdgeIterator e_it;
  };
  
  
  template <typename VDesc, typename ContIterator>
  class adjlistBC_inv_adjacent_viter : 
    public iterator_facade<
      adjlistBC_inv_adjacent_viter<VDesc, ContIterator>,
      const VDesc,
      typename std::iterator_traits<ContIterator>::iterator_category
    > {
    public:
      typedef adjlistBC_inv_adjacent_viter<VDesc, ContIterator> self;
      
      explicit adjlistBC_inv_adjacent_viter(ContIterator aEdgeIt = ContIterator()) : e_it(aEdgeIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++e_it; };
      void decrement() { --e_it; };
      bool equal(const self& rhs) const { return (this->e_it == rhs.e_it); };
      const VDesc& dereference() const { return e_it->source; };
      
      void advance(std::ptrdiff_t i) { std::advance(this->e_it, i); };
      std::ptrdiff_t distance_to(const self& rhs) const { return std::distance(this->e_it, rhs.e_it); }; 
      
      ContIterator e_it;
  };
  
  
  
  
  
  
  
  template <typename Container, typename Desc>
  int adjlistBC_is_desc_at_begin(Container& cont, Desc d) {
    if(d == cont.begin())
      return 1;
    return 0;
  };
  
  template <typename Container, typename Desc>
  int adjlistBC_is_desc_at_begin(Container* cont, Desc d) {
    return adjlistBC_is_desc_at_begin(*cont, d);
  };
  
  template <typename ValueType, typename Desc>
  int adjlistBC_is_desc_at_begin( ::boost::container::vector<ValueType>& cont, Desc d) {
    if(d == 0)
      return 1;
    return 0;
  };
  
  template <typename ValueType, typename Desc>
  int adjlistBC_is_desc_at_begin( BC_pooled_vector<ValueType>& cont, Desc d) {
    if(d == 0)
      return 1;
    return 0;
  };
  
  
  template <typename Container, typename Desc>
  int adjlistBC_check_desc_validity(Container& cont, Desc d) {
    if(d == cont.end())
      return 1;
    return 0;
  };
  
  template <typename Container, typename Desc>
  int adjlistBC_check_desc_validity(Container* cont, Desc d) {
    return adjlistBC_check_desc_validity(*cont, d);
  };
  
  template <typename ValueType, typename Desc>
  int adjlistBC_check_desc_validity( ::boost::container::vector<ValueType>& cont, Desc d) {
    if(d >= cont.size())
      return 1;
    return 0;
  };
  
  template <typename ValueType, typename Desc>
  int adjlistBC_check_desc_validity( BC_pooled_vector<ValueType>& cont, Desc d) {
    if(d >= cont.size())
      return 1;
    if(cont.m_data[d].which() == 1)
      return 2;
    return 0;
  };
  
  
  template <typename VContainer, typename EDesc>
  class adjlistBC_edge_iterator :
    public iterator_facade<
      adjlistBC_edge_iterator<VContainer, EDesc>,
      const EDesc,
      std::bidirectional_iterator_tag
    > {
    public:
      typedef adjlistBC_edge_iterator<VContainer, EDesc> self;
      
      adjlistBC_edge_iterator() : p_cont(NULL), e() { };
      
      static self begin(VContainer& vcont) {
        typedef typename EDesc::edge_id_type RawEDesc;
        self result;
        result.p_cont = &vcont;
        result.e.source = BC_get_begin_desc(vcont);
        while(true) {
          while(adjlistBC_check_desc_validity(*p_cont, result.e.source) == 2)
            ++(result.e.source);
          if(adjlistBC_check_desc_validity(*p_cont, result.e.source) == 1) {
            result.e.edge_id = RawEDesc(); // end iterator (also begin iterator, vertex range is empty).
            break;
          };
          result.e.edge_id = BC_get_begin_desc(BC_desc_to_value(*p_cont, result.e.source).out_edges);
          while(adjlistBC_check_desc_validity(BC_desc_to_value(*p_cont, result.e.source).out_edges, result.e.edge_id) == 2)
            ++(result.e.edge_id);
          if(adjlistBC_check_desc_validity(BC_desc_to_value(*p_cont, result.e.source).out_edges, result.e.edge_id) == 1) {
            ++(result.e.source);
            continue;
          };
          break;
        };
        return result;
      };
      
      static self end(VContainer& vcont) {
        typedef typename EDesc::edge_id_type RawEDesc;
        self result;
        result.p_cont = &vcont;
        result.e.source = BC_get_end_desc(vcont);
        result.e.edge_id = RawEDesc();
        return result;
      };
      
    public:
      friend class iterator_core_access;
      
      void increment() {
        typedef typename EDesc::edge_id_type RawEDesc;
        if(adjlistBC_check_desc_validity(*p_cont, e.source))
          return; // iterator at the end.
        
        ++(e.edge_id);
        while(true) {
          while(adjlistBC_check_desc_validity(BC_desc_to_value(*p_cont, e.source).out_edges, e.edge_id) == 2)
            ++(e.edge_id);
          if(adjlistBC_check_desc_validity(BC_desc_to_value(*p_cont, e.source).out_edges, e.edge_id) == 0)
            return;
          
          ++(e.source);
          while(adjlistBC_check_desc_validity(*p_cont, e.source) == 2)
            ++(e.source);
          if(adjlistBC_check_desc_validity(*p_cont, e.source) == 1) {
            e.edge_id = RawEDesc(); // end iterator.
            return;
          };
          e.edge_id = BC_get_begin_desc(BC_desc_to_value(*p_cont, e.source).out_edges);
        };
      };
      void decrement() { 
        typedef typename EDesc::edge_id_type RawEDesc;
        typedef typename EDesc::source_descriptor SrcDesc;
        
        EDesc tmp_e = e; // apply decrement to temporary descriptor, in case e is at the beginning.
        SrcDesc last_source = e.source;
        
        while(true) {
          while(adjlistBC_check_desc_validity(*p_cont, tmp_e.source)) { // while at the end vertex.
            if(adjlistBC_is_desc_at_begin(*p_cont, tmp_e.source))
              return; // the vertex container is empty.
            --(tmp_e.source);
          };
          if(tmp_e.source != last_source) // if we moved back the vertex
            tmp_e.edge_id = BC_get_end_desc(BC_desc_to_value(*p_cont, tmp_e.source).out_edges);
          
          bool edge_range_was_empty = false;
          while(adjlistBC_check_desc_validity(BC_desc_to_value(*p_cont, tmp_e.source).out_edges, tmp_e.edge_id)) {
            // check if out-edge range is empty:
            if(adjlistBC_is_desc_at_begin(BC_desc_to_value(*p_cont, tmp_e.source).out_edges, tmp_e.edge_id)) {
              if(adjlistBC_is_desc_at_begin(*p_cont, tmp_e.source))
                return; // reached the beginning without finding a valid edge.
              --(tmp_e.source);
              edge_range_was_empty = true;
              break; // go back to scan for first valid source vertex.
            };
            --(tmp_e.edge_id);
          };
          if(!edge_range_was_empty)
            break;
        };
        
        e = tmp_e; // we have found a valid descriptor.
      };
      bool equal(const self& rhs) const { 
        return (this->e == rhs.e);
      };
      const EDesc& dereference() const { 
        return e;
      };
      
      VContainer* p_cont;
      EDesc e;
  };
  
  
  
  
}; // namespace detail


}; // namespace graph


}; // namespace boost


#endif


















