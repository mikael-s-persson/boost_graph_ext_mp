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

#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/detail/boost_container_generators.hpp>

#include <iterator>

namespace boost {

namespace graph { 

namespace detail {
  
  
  
  template <typename DirectedS>
  struct adjlistBC_traversal_tag :
    public virtual vertex_list_graph_tag,
    public virtual incidence_graph_tag,
    public virtual adjacency_graph_tag,
    public virtual edge_list_graph_tag { };
  
  template <>
  struct adjlistBC_traversal_tag< undirectedS > :
    public virtual vertex_list_graph_tag,
    public virtual incidence_graph_tag,
    public virtual adjacency_graph_tag,
    public virtual edge_list_graph_tag,
    public virtual bidirectional_graph_tag { };
  
  template <>
  struct adjlistBC_traversal_tag< bidirectionalS > :
    public virtual vertex_list_graph_tag,
    public virtual incidence_graph_tag,
    public virtual adjacency_graph_tag,
    public virtual edge_list_graph_tag,
    public virtual bidirectional_graph_tag { };
  
  
  
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
  
  
  template <typename EDesc, typename IEIter, typename OEIter>
  class adjlistBC_undir_ioeiter : 
    public iterator_facade<
      adjlistBC_undir_ioeiter<EDesc, IEIter, OEIter>,
      EDesc,
      std::bidirectional_iterator_tag,
      EDesc
    > {
    public:
      typedef adjlistBC_undir_ioeiter<EDesc, IEIter, OEIter> self;
      
      adjlistBC_undir_ioeiter() : view_as_out_edges(true), cur_ie(), ie_end(), oe_beg(), cur_oe(), oe_end() { };
      
      adjlistBC_undir_ioeiter(bool aViewAsOutEdges, 
                              IEIter aCurIE, IEIter aIEEnd, 
                              OEIter aOEBeg, OEIter aOEEnd) : 
                              view_as_out_edges(aViewAsOutEdges),
                              cur_ie(aCurIE), ie_end(aIEEnd),
                              oe_beg(aOEBeg), cur_oe(aOEBeg), oe_end(aOEEnd) { };
      
      adjlistBC_undir_ioeiter(bool aViewAsOutEdges, 
                              IEIter aIEEnd, OEIter aOEBeg, OEIter aCurOE, OEIter aOEEnd) : 
                              view_as_out_edges(aViewAsOutEdges),
                              cur_ie(aIEEnd), ie_end(aIEEnd),
                              oe_beg(aOEBeg), cur_oe(aCurOE), oe_end(aOEEnd) { };
      
      adjlistBC_undir_ioeiter(bool aViewAsOutEdges, 
                              IEIter aIEEnd, OEIter aOEBeg, OEIter aOEEnd) : 
                              view_as_out_edges(aViewAsOutEdges),
                              cur_ie(aIEEnd), ie_end(aIEEnd),
                              oe_beg(aOEBeg), cur_oe(aOEEnd), oe_end(aOEEnd) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() {
        if(cur_ie != ie_end)
          ++cur_ie;
        else if(cur_oe != oe_end)
          ++cur_oe;
      };
      void decrement() { 
        if( cur_oe == oe_beg )
          --cur_ie;
        else
          --cur_oe;
      };
      bool equal(const self& rhs) const { return (this->cur_ie == rhs.cur_ie) && (this->cur_oe == rhs.cur_oe); };
      EDesc dereference() const {
        if(view_as_out_edges) {
          if( cur_ie != ie_end )
            return EDesc(*cur_ie,true);
          return EDesc(*cur_oe,false);
        } else {
          if( cur_ie != ie_end )
            return EDesc(*cur_ie,false);
          return EDesc(*cur_oe,true);
        };
      };
      
      bool view_as_out_edges;
      IEIter cur_ie;
      IEIter ie_end;
      OEIter oe_beg;
      OEIter cur_oe;
      OEIter oe_end;
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
    if(d >= cont.m_data.size())
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
        result.e.source = BC_iterator_to_desc(vcont, BC_get_begin_iter(vcont));
        while(true) {
          while(adjlistBC_check_desc_validity(*result.p_cont, result.e.source) == 2)
            ++(result.e.source);
          if(adjlistBC_check_desc_validity(*result.p_cont, result.e.source) == 1) {
            result.e.edge_id = RawEDesc(); // end iterator (also begin iterator, vertex range is empty).
            break;
          };
          result.e.edge_id = BC_iterator_to_desc(BC_desc_to_value(*result.p_cont, result.e.source).out_edges, BC_get_begin_iter(BC_desc_to_value(*result.p_cont, result.e.source).out_edges));
          while(adjlistBC_check_desc_validity(BC_desc_to_value(*result.p_cont, result.e.source).out_edges, result.e.edge_id) == 2)
            ++(result.e.edge_id);
          if(adjlistBC_check_desc_validity(BC_desc_to_value(*result.p_cont, result.e.source).out_edges, result.e.edge_id) == 1) {
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
        result.e.source = BC_iterator_to_desc(vcont, BC_get_end_iter(vcont));
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
          e.edge_id = BC_iterator_to_desc(BC_desc_to_value(*p_cont, e.source).out_edges, BC_get_begin_iter(BC_desc_to_value(*p_cont, e.source).out_edges));
        };
      };
      void decrement() { 
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
            tmp_e.edge_id = BC_iterator_to_desc(BC_desc_to_value(*p_cont, tmp_e.source).out_edges, BC_get_end_iter(BC_desc_to_value(*p_cont, tmp_e.source).out_edges));
          
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
  
  
  
  
  template <typename EDesc, typename EIter>
  class adjlistBC_undir_eiter : 
    public iterator_facade<
      adjlistBC_undir_eiter<EDesc, EIter>,
      EDesc,
      std::bidirectional_iterator_tag,
      EDesc
    > {
    public:
      typedef adjlistBC_undir_eiter<EDesc, EIter> self;
      
      adjlistBC_undir_eiter(EIter aCurIt = EIter()) : cur_it(aCurIt) { };
      
    public: // private:
      friend class iterator_core_access;
      
      void increment() { ++cur_it; };
      void decrement() { --cur_it; };
      bool equal(const self& rhs) const { return (this->cur_it == rhs.cur_it); };
      EDesc dereference() const { return EDesc(*cur_it); };
      
      EIter cur_it;
  };
  
  
  
  
}; // namespace detail


}; // namespace graph


}; // namespace boost


#endif


















