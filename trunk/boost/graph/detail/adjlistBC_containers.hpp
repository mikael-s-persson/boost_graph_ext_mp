// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file adjlistBC_containers.hpp
 * 
 * This library implements adjacency-list data structures that underpin the adjacency_list_v2 template. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_ADJLISTBC_CONTAINERS_HPP
#define BOOST_ADJLISTBC_CONTAINERS_HPP

#include <boost/config.hpp>

#include <boost/graph/detail/boost_container_generators.hpp>
#include <boost/graph/detail/adjlistBC_iterators.hpp>

#include <boost/graph/detail/container_generators.hpp>  // for vecBC, poolBC, listBC, ... and BC_container_gen.
#include <boost/graph/graph_selectors.hpp>              // for directedS, undirectedS, bidirectionalS.

#include <iterator>
#include <utility>

#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/not.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>


namespace boost {

namespace graph { 

namespace detail {
  
  
  
/*************************************************************************
 *        value-types for vertices and edges
 * **********************************************************************/
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_vertex_stored_type;  // forward declaration.
  
  
  // this is fine because boost-containers allow incomplete types:
  template <typename VertexListS, typename OutEdgeListS,  typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_vertex_config {
    typedef adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> stored_type;
    typedef typename BC_select_value_type<VertexListS, stored_type >::type value_type;
    typedef typename BC_container_gen<VertexListS, value_type >::type container;
    typedef typename BC_select_descriptor<VertexListS, typename container::iterator >::type descriptor;
    
    static descriptor null_vertex() { return BC_null_desc<descriptor>::value(); };
  };
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_edge_stored_type {
    typedef adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> self;
    typedef adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::descriptor vertex_descriptor;
    
    vertex_descriptor target;
    EdgeProperties data;
    
    adjlistBC_edge_stored_type(vertex_descriptor aTarget = VConfig::null_vertex()) : target(aTarget), data() { };
    adjlistBC_edge_stored_type(vertex_descriptor aTarget, const EdgeProperties& aData) : target(aTarget), data(aData) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    adjlistBC_edge_stored_type(vertex_descriptor aTarget, EdgeProperties&& aData) : target(aTarget), data(std::move(aData)) { };
#endif
    
    bool operator<(const self& rhs) const { return BC_desc_less_than(this->target, rhs.target); };
    bool operator<=(const self& rhs) const { return !BC_desc_less_than(rhs.target, this->target); };
    bool operator>(const self& rhs) const { return BC_desc_less_than(rhs.target, this->target); };
    bool operator>=(const self& rhs) const { return !BC_desc_less_than(this->target, rhs.target); };
    bool operator==(const self& rhs) const { return (this->target == rhs.target); };
    bool operator!=(const self& rhs) const { return (this->target != rhs.target); };
  };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  std::size_t hash_value(const adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& ep) {
    typedef adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>::vertex_descriptor Vertex;
    ::boost::hash<Vertex> hasher;
    return hasher(ep.target);
  };
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_edge_config {
    typedef typename adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>::descriptor vertex_descriptor;
    
    typedef adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> stored_type;
    typedef typename BC_select_value_type<OutEdgeListS, stored_type >::type value_type;
    typedef typename BC_container_gen<OutEdgeListS, value_type >::type container;
    typedef typename BC_select_descriptor<OutEdgeListS, typename container::iterator >::type raw_descriptor;
    typedef typename mpl::if_< 
      mpl::and_<
        is_same< vertex_descriptor, std::size_t >,
        mpl::not_< is_same< raw_descriptor, std::size_t > > 
      >,
      container*, container >::type container_ptr;
    typedef BC_edge_desc<vertex_descriptor, raw_descriptor> descriptor;
    
    static descriptor null_edge() { return BC_null_desc<descriptor>::value(); };
  };
  
  
/*************************************************************************
 *        vertex values, as stored in the vertex containers.
 * **********************************************************************/
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_vertex_stored_type {
    typedef adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> self;
    typedef DirectedS directed_tag;
    typedef adjlistBC_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> Config;
    typedef typename Config::container edge_container;
    typedef typename Config::container_ptr edge_container_ptr;
    typedef typename Config::descriptor edge_descriptor;
    typedef ::boost::container::vector<edge_descriptor> in_edge_container;
    
    VertexProperties data;
    edge_container_ptr out_edges;
    in_edge_container in_edges;
    
    adjlistBC_vertex_stored_type() : data(), out_edges(), in_edges() { };
    adjlistBC_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges(), in_edges() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    adjlistBC_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges(), in_edges() { };
#endif
    
    bool operator<(const self& rhs) const { return (this->data < rhs.data); };
    bool operator<=(const self& rhs) const { return !(rhs.data < this->data); };
    bool operator>(const self& rhs) const { return (rhs.data < this->data); };
    bool operator>=(const self& rhs) const { return !(this->data < rhs.data); };
    bool operator==(const self& rhs) const { return (this->data == rhs.data); };
    bool operator!=(const self& rhs) const { return !(this->data == rhs.data); };
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> {
    typedef adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> self;
    typedef directedS directed_tag;
    typedef adjlistBC_edge_config<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> Config;
    typedef typename Config::container edge_container;
    typedef typename Config::container_ptr edge_container_ptr;
    typedef typename Config::descriptor edge_descriptor;
    
    VertexProperties data;
    edge_container_ptr out_edges;
    
    adjlistBC_vertex_stored_type() : data(), out_edges() { };
    adjlistBC_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    adjlistBC_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges() { };
#endif
    
    bool operator<(const self& rhs) const { return (this->data < rhs.data); };
    bool operator<=(const self& rhs) const { return !(rhs.data < this->data); };
    bool operator>(const self& rhs) const { return (rhs.data < this->data); };
    bool operator>=(const self& rhs) const { return !(this->data < rhs.data); };
    bool operator==(const self& rhs) const { return (this->data == rhs.data); };
    bool operator!=(const self& rhs) const { return !(this->data == rhs.data); };
  };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
            typename VertexProperties, typename EdgeProperties>
  std::size_t hash_value(const adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp) {
    ::boost::hash<VertexProperties> hasher;
    return hasher(vp.data);
  };
  
  
  
/*************************************************************************
 *        out-edges container pointers (or not): create and destroy
 * **********************************************************************/
  
  template <typename OutEdgeCont>
  struct adjlistBC_out_edges_factory {
    template <typename VertexValue>
    static void create_out_edges(VertexValue&) { };
    template <typename VertexValue>
    static void destroy_out_edges(VertexValue&) { };
    template <typename VertexContainer>
    static void destroy_all_out_edges(VertexContainer&) { };
  };
  
  template <typename OutEdgeCont>
  struct adjlistBC_out_edges_factory<OutEdgeCont*> {
    template <typename VertexValue>
    static void create_out_edges(VertexValue& vp) {
      vp.out_edges = new OutEdgeCont();
    };
    template <typename VertexValue>
    static void destroy_out_edges(VertexValue& vp) {
      delete vp.out_edges;
    };
    template <typename VertexValue>
    static void destroy_all_out_edges( ::boost::container::vector<VertexValue>& vcont) {
      typedef ::boost::container::vector<VertexValue> VertexContainer;
      typedef typename VertexContainer::iterator Iter;
      for(Iter vi = vcont.begin(); vi != vcont.end(); ++vi)
        delete vi->out_edges;
    };
    template <typename VertexValue>
    static void destroy_all_out_edges(BC_pooled_vector<VertexValue>& vcont) {
      typedef typename BC_pooled_vector<VertexValue>::container_type VertexContainer;
      typedef typename VertexContainer::iterator Iter;
      for(Iter vi = vcont.m_data.begin(); vi != vcont.m_data.end(); ++vi) {
        if(vi->which() == 0)
          delete get<VertexValue>(*vi).out_edges;
      };
    };
  };
  
  template <typename OutEdgeCont>
  struct adjlistBC_out_edges_range {
    typedef typename OutEdgeCont::iterator Iter;
    static Iter begin(OutEdgeCont& econt) { return econt.begin(); };
    static Iter end(OutEdgeCont& econt) { return econt.end(); };
    template <typename Desc>
    static Iter from_desc(OutEdgeCont& econt, Desc d) { return BC_desc_to_iterator(econt, d); };
  };
  
  template <typename OutEdgeCont>
  struct adjlistBC_out_edges_range<OutEdgeCont*> {
    typedef typename OutEdgeCont::iterator Iter;
    static Iter begin(OutEdgeCont* econt) { return econt->begin(); };
    static Iter end(OutEdgeCont* econt) { return econt->end(); };
    template <typename Desc>
    static Iter from_desc(OutEdgeCont* econt, Desc d) { return BC_desc_to_iterator(*econt, d); };
  };
  
  
  
  
  
/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:   O(E/V)      O(E/V)      O(E/V)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)

  
  template <typename VertexDesc>
  struct adjlistBC_is_in_edge_of {
    VertexDesc v;
    explicit adjlistBC_is_in_edge_of(VertexDesc aV) : v(aV) { };
    template <typename EdgeProp>
    bool operator()(const EdgeProp& ep) const {
      return (BC_get_value(ep).target == v);
    };
  };
  template <typename VertexDesc>
  adjlistBC_is_in_edge_of<VertexDesc> adjlistBC_make_is_in_edge_of(VertexDesc aV) {
    return adjlistBC_is_in_edge_of<VertexDesc>(aV);
  };
  
  // for OutEdgeListS = listS, setS, ...
  template <typename Container, typename Predicate, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge_if(Container& cont, Predicate pred, VertexCont&, VertexDesc) {
    typedef typename Container::iterator Iter;
    for(Iter it = cont.begin(); it != cont.end(); ) {
      if( pred(*it) )
        it = cont.erase(it);
      else
        ++it;
    };
  };
  
  template <typename Container, typename Predicate, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge_if(Container* cont, Predicate pred, VertexCont& vcont, VertexDesc v) {
    adjlistBC_erase_edge_if(*cont, pred, vcont, v);
  };
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename VertexDesc>
  void adjlistBC_update_in_edge_id(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                                   VertexDesc v, std::size_t old_id, std::size_t new_id) { };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename VertexDesc>
  void adjlistBC_update_in_edge_id(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                                   VertexDesc v, std::size_t old_id, std::size_t new_id) {
    typedef adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VertexValue;
    typedef typename VertexValue::in_edge_container InEdgeCont;
    typedef typename InEdgeCont::iterator InEdgeIter;
    
    for(InEdgeIter ei = vp.in_edges.begin(); ei != vp.in_edges.end(); ++ei) {
      if( (ei->source == v) && (ei->edge_id == old_id) ) {
        ei->edge_id = new_id;
        break;
      };
    };
  };
  
  // for OutEdgeListS = vecS
  template <typename ValueType, typename Predicate, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge_if( ::boost::container::vector<ValueType>& cont, Predicate pred, VertexCont& vertex_cont, VertexDesc v) {
    typedef ::boost::container::vector<ValueType> Container;
    typedef typename Container::iterator Iter;
    using std::swap;
    
    Iter it_last = cont.end();
    for(Iter it = cont.begin(); it != it_last; ) {
      if( pred(*it) ) {
        --it_last;
        if(it != it_last) {
          swap(*it, *it_last);
          // If this graph has in-edge references, then they must be updated.
          adjlistBC_update_in_edge_id(BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), 
                                      v, it_last - cont.begin(), it - cont.begin());
        };
      } else
        ++it;
    };
    cont.erase(it_last, cont.end());
  };
  
  // for OutEdgeListS = poolS
  template <typename ValueType, typename Predicate, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge_if( BC_pooled_vector<ValueType>& cont, Predicate pred, VertexCont&, VertexDesc) {
    typedef typename BC_pooled_vector<ValueType>::container_type Container;
    typedef typename Container::iterator Iter;
    
    for(Iter it = cont.m_data.begin(); it != cont.m_data.end(); ++it) {
      if( (it->which() == 0) && pred(get<ValueType>(*it)) ) {
        *it = cont.m_first_hole;
        cont.m_first_hole = BC_hole_desc(it - cont.m_data.begin());
        --(cont.m_num_elements);
      };
    };
  };
  
  
  
/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))
  
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void adjlistBC_erase_in_edge(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                               EdgeDesc e) { };
  
  // for vector of edge-desc (in-edges)
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void adjlistBC_erase_in_edge(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                               EdgeDesc e) {
    vp.in_edges.erase( std::find(vp.in_edges.begin(), vp.in_edges.end(), e) );
  };
  
  
  // for OutEdgeListS = listS, setS, ...
  template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge(Container& cont, EdgeDesc e, VertexCont&, VertexDesc) {
    cont.erase(e.edge_id);
  };
  template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge(Container* cont, EdgeDesc e, VertexCont& vcont, VertexDesc v) {
    adjlistBC_erase_edge(*cont, e, vcont, v);
  };
  
  // for OutEdgeListS = vecS
  template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge( ::boost::container::vector<ValueType>& cont, EdgeDesc e, VertexCont& vertex_cont, VertexDesc v) {
    typedef ::boost::container::vector<ValueType> Container;
    typedef typename Container::iterator Iter;
    using std::swap;
    
    Iter it = BC_desc_to_iterator(cont, e.edge_id);
    Iter it_last = cont.end(); --it_last;
    if(it != it_last) {
      swap(*it, *it_last);
      // If this graph has in-edge references, then they must be updated.
      adjlistBC_update_in_edge_id(BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), 
                                  v, it_last - cont.begin(), it - cont.begin());
    };
    cont.erase(it_last, cont.end());
  };
  
  // for OutEdgeListS = poolS
  template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void adjlistBC_erase_edge( BC_pooled_vector<ValueType>& cont, EdgeDesc e, VertexCont&, VertexDesc) {
    cont.m_data[e.edge_id] = cont.m_first_hole;
    cont.m_first_hole = BC_hole_desc(e.edge_id);
    --(cont.m_num_elements);
  };
  
  
  
/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)
  
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void adjlistBC_add_in_edge(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) { };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void adjlistBC_add_in_edge(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) {
    vp.in_edges.push_back(e);
  };
  
  
  // for OutEdgeListS = listS, multisetS, ...
  template <typename Container, typename EdgeProperties, typename VertexDesc>
  std::pair<typename Container::iterator, bool> adjlistBC_add_edge(Container& cont, const EdgeProperties& ep, VertexDesc v) {
    typedef typename Container::value_type ValueType;
    return std::pair<typename Container::iterator, bool>(cont.insert(cont.end(), ValueType(v, ep)), true);
  };
  
  // for OutEdgeListS = setS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<typename ::boost::container::set<ValueType>::iterator, bool> adjlistBC_add_edge( ::boost::container::set<ValueType>& cont, const EdgeProperties& ep, VertexDesc v) {
    return cont.insert(ValueType(v, ep));
  };
  
  // for OutEdgeListS = unordered_setS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<typename ::boost::unordered_set<ValueType>::iterator, bool> adjlistBC_add_edge( ::boost::unordered_set<ValueType>& cont, const EdgeProperties& ep, VertexDesc v) {
    return cont.insert(ValueType(v, ep));
  };
  
  template <typename Container, typename EdgeProperties, typename VertexDesc>
  std::pair<typename Container::iterator, bool> adjlistBC_add_edge(Container* cont, const EdgeProperties& ep, VertexDesc v) {
    return adjlistBC_add_edge(*cont, ep, v);
  };
  
  // for OutEdgeListS = vecS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<std::size_t, bool> adjlistBC_add_edge( ::boost::container::vector<ValueType>& cont, const EdgeProperties& ep, VertexDesc v) {
    cont.push_back(ValueType(v, ep));
    return std::pair<std::size_t, bool>(cont.size() - 1, true);
  };
  
  // for OutEdgeListS = poolS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<std::size_t, bool> adjlistBC_add_edge( BC_pooled_vector<ValueType>& cont, const EdgeProperties& ep, VertexDesc v) {
    typedef typename BC_pooled_vector<ValueType>::iterator Iter;
    
    if(cont.m_first_hole == BC_hole_desc()) {
      Iter it = cont.m_data.insert(cont.m_data.end(), ValueType(v, ep));
      ++(cont.m_num_elements);
      return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
    } else {
      Iter it = cont.m_data.begin() + cont.m_first_hole.value;
      cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole]);
      *it = ValueType(v, ep);
      ++(cont.m_num_elements);
      return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
    };
  };
  
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  
  // for OutEdgeListS = listS, multisetS, ...
  template <typename Container, typename EdgeProperties, typename VertexDesc>
  std::pair<typename Container::iterator, bool> adjlistBC_add_edge(Container& cont, EdgeProperties&& ep, VertexDesc v) {
    typedef typename Container::value_type ValueType;
    return std::pair<typename Container::iterator, bool>(cont.insert(cont.end(), ValueType(v, std::move(ep))), true);
  };
  
  // for OutEdgeListS = setS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<typename ::boost::container::set<ValueType>::iterator, bool> adjlistBC_add_edge( ::boost::container::set<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
    return cont.insert(ValueType(v, std::move(ep)));
  };
  
  // for OutEdgeListS = unordered_setS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<typename ::boost::unordered_set<ValueType>::iterator, bool> adjlistBC_add_edge( ::boost::unordered_set<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
    return cont.insert(ValueType(v, std::move(ep)));
  };
  
  template <typename Container, typename EdgeProperties, typename VertexDesc>
  std::pair<typename Container::iterator, bool> adjlistBC_add_edge(Container* cont, EdgeProperties&& ep, VertexDesc v) {
    return adjlistBC_add_edge(*cont, std::move(ep), v);
  };
  
  // for OutEdgeListS = vecS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<std::size_t, bool> adjlistBC_add_edge( ::boost::container::vector<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
    cont.push_back(ValueType(v, std::move(ep)));
    return std::pair<std::size_t, bool>(cont.size() - 1, true);
  };
  
  // for OutEdgeListS = poolS
  template <typename ValueType, typename EdgeProperties, typename VertexDesc>
  std::pair<std::size_t, bool> adjlistBC_add_edge( BC_pooled_vector<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
    typedef typename BC_pooled_vector<ValueType>::iterator Iter;
    
    if(cont.m_first_hole == BC_hole_desc()) {
      Iter it = cont.m_data.insert(cont.m_data.end(), ValueType(v, std::move(ep)));
      ++(cont.m_num_elements);
      return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
    } else {
      Iter it = cont.m_data.begin() + cont.m_first_hole.value;
      cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole]);
      *it = ValueType(v, std::move(ep));
      ++(cont.m_num_elements);
      return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
    };
  };
  
#endif
  
  
  
  
/*************************************************************************
 *        functions for clearing the edges of a vertex
 * **********************************************************************/
  // NOTE: The function 'adjlistBC_clear_vertex' works for all graph types.
  
// NOTE: Time complexities: 
//                 vecBC         poolBC         listBC
// directedS:   O(V*(E/V))     O(V*(E/V))     O(V*(E/V))
// bidir:       O((E/V)^2)     O((E/V)^2)     O((E/V)^2)
  
  
  template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type adjlistBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v) {
    typedef typename VertexCont::iterator VertexIter;
    
    // first, just clear the out-going edges. No need to synchronize in-edges (there are none).
    vp.out_edges.clear();
    VertexIter vi = BC_desc_to_iterator(cont, v);
    
    // now, the stupid part... we have to traverse all other vertices (and their edges) 
    // to look for in-edges that lead to "v", and erase them.
    for(VertexIter ui = cont.begin(); ui != cont.end(); ++ui) {
      if((ui == vi) || (!adjlistBC_is_elem_valid(*ui)))
        continue;
      VertexValue& up = BC_get_value(*ui);
      adjlistBC_erase_edge_if(up.out_edges, adjlistBC_make_is_in_edge_of(v), cont, u);
    };
  };
  
  template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type adjlistBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v) {
    typedef typename VertexValue::edge_container OutEdgeCont;
    typedef typename OutEdgeCont::iterator OutEdgeIter;
    typedef typename VertexValue::in_edge_container InEdgeCont;
    typedef typename InEdgeCont::iterator InEdgeIter;
    
    // first, remove the in-edge references from the adjacent vertices of v.
    for(OutEdgeIter ei = vp.out_edges.begin(); ei != vp.out_edges.end(); ++ei) {
      if( !adjlistBC_is_elem_valid(*ei) )
        continue;
      VertexValue& wp = BC_get_value( *BC_desc_to_iterator(cont, BC_get_value(*ei).target) );
      for(InEdgeIter iei = wp.in_edges.begin(); iei != wp.in_edges.end(); ++iei) {
        if( (iei->source == v) && (BC_desc_to_iterator(vp.out_edges, iei->edge_id) == ei) ) {
          wp.in_edges.erase(iei);
          break;
        };
      };
    };
    
    // then, clear the out-going edges.
    vp.out_edges.clear();
    
    // finally, remove the required out-edges of the "parent" vertices of v.
    for(InEdgeIter iei = vp.in_edges.begin(); iei != vp.in_edges.end(); ++iei) {
      VertexDesc u = iei->source;
      VertexValue& up = BC_get_value( *BC_desc_to_iterator(cont, iei->source) );
      adjlistBC_erase_edge_if(up.out_edges, adjlistBC_make_is_in_edge_of(v), cont, u);
    };
    
  };
  
  
  
  
  
  
/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
  // NOTE: The function 'adjlistBC_erase_vertex' works for all graph types.
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

  
  template <typename Container, typename VertexDesc>
  bool adjlistBC_erase_vertex(Container& cont, VertexDesc v) {
    typedef typename Container::iterator Iter;
    typedef typename Container::value_type VertexValue;
    typedef adjlistBC_out_edges_factory< typename VertexValue::edge_container_ptr > OutEdgeFactory;
    OutEdgeFactory::destroy_out_edges(*v);
    cont.erase(v);
  };
  
  // O(E)
  template <typename DirectedS, typename ValueType>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type adjlistBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    typedef ::boost::container::vector<ValueType> VContainer;
    typedef typename VContainer::iterator VIter;
    typedef typename ValueType::edge_container OutEdgeCont;
    typedef adjlistBC_out_edges_range< typename ValueType::edge_container_ptr > OERange;
    typedef typename OutEdgeCont::iterator OEIter;
    typedef typename ValueType::edge_descriptor EdgeDesc;
    
    for(VIter vi = cont.begin(); vi != cont.end(); ++vi)
      for(OEIter ei = OERange::begin(vi->out_edges); ei != OERange::end(vi->out_edges); ++ei)
        if(BC_get_value(*ei).target == old_v_id)
          BC_get_value(*ei).target = new_v_id;
    
  };
  
  // O((E/V)^2)
  template <typename DirectedS, typename ValueType>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type adjlistBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    typedef ::boost::container::vector<ValueType> VContainer;
    typedef typename VContainer::iterator VIter;
    typedef typename ValueType::edge_container OutEdgeCont;
    typedef adjlistBC_out_edges_range< typename ValueType::edge_container_ptr > OERange;
    typedef typename OutEdgeCont::iterator OEIter;
    typedef typename ValueType::edge_descriptor EdgeDesc;
    typedef typename ValueType::in_edge_container InEdgeCont;
    typedef typename InEdgeCont::iterator InEdgeIter;
    
    // first, update in-edge vertices
    for(InEdgeIter iei = cont[new_v_id].in_edges.begin(); iei != cont[new_v_id].in_edges.end(); ++iei) {
      ValueType& up = cont[iei->source];
      EdgeDesc ed = OERange::from_desc(up.out_edges, iei->edge_id);
      for(OEIter ei = OERange::begin(up.out_edges); ei != OERange::end(up.out_edges); ++ei) {
        if((BC_get_value(*ei).target == old_v_id) && (ed == ei)) {
          BC_get_value(*ei).target = new_v_id;
          break;
        };
      };
    };
    
    // second, update out-edge vertices
    for(OEIter ei = OERange::begin(cont[new_v_id].out_edges); ei != OERange::end(cont[new_v_id].out_edges); ++ei) {
      ValueType& wp = cont[BC_get_value(*ei).target];
      for(InEdgeIter iei = wp.in_edges.begin(); iei != wp.in_edges.end(); ++iei) {
        if((iei->source == old_v_id) && 
           (ei == OERange::from_desc(cont[new_v_id].out_edges, iei->edge_id))) {
          iei->source = new_v_id;
          break;
        };
      };
    };
    
  };
  
  template <typename ValueType>
  void adjlistBC_erase_vertex( ::boost::container::vector<ValueType>& cont, std::size_t v) {
    typedef ::boost::container::vector<ValueType> Container;
    typedef typename Container::iterator Iter;
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OutEdgeFactory;
    typedef typename ValueType::directed_tag DirectedS;
    using std::swap;
    
    Iter it = BC_desc_to_iterator(cont, v);
    OutEdgeFactory::destroy_out_edges(*it);
    Iter it_last = cont.end(); --it_last;
    if(it == it_last) {
      cont.erase(it_last);
      return;
    };
    swap(*it, *it_last);
    std::size_t old_id = it_last - cont.begin();
    std::size_t new_id = it - cont.begin();
    cont.erase(it_last);
    adjlistBC_update_out_edges<DirectedS>(cont, old_id, new_id);
  };
  
  template <typename ValueType>
  void adjlistBC_erase_vertex(BC_pooled_vector<ValueType>& cont, std::size_t v) {
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OutEdgeFactory;
    // the first_hole will become v and v will be referring to first_hole (i.e., like a simple linked-list):
    OutEdgeFactory::destroy_out_edges(get<ValueType>(cont.m_data[v]));
    cont.m_data[v] = cont.m_first_hole;
    cont.m_first_hole = BC_hole_desc(v);
    --(cont.m_num_elements);
  };
  
  
  
  
/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

  
  template <typename Container, typename VertexProperties>
  typename Container::iterator adjlistBC_add_vertex(Container& cont, const VertexProperties& vp) {
    typedef typename Container::value_type ValueType;
    typedef typename Container::iterator Iter;
    Iter it = cont.insert(cont.end(), ValueType(vp));
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OEFactory;
    OEFactory::create_out_edges(*it);
    return it;
  };
  
  template <typename Container, typename VertexProperties>
  typename Container::iterator adjlistBC_add_vertex(Container* cont, const VertexProperties& vp) {
    return adjlistBC_add_vertex(*cont, vp);
  };
  
  template <typename ValueType, typename VertexProperties>
  std::size_t adjlistBC_add_vertex( ::boost::container::vector<ValueType>& cont, const VertexProperties& vp) {
    typedef typename ::boost::container::vector<ValueType>::iterator Iter;
    Iter it = cont.insert(cont.end(), ValueType(vp));
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OEFactory;
    OEFactory::create_out_edges(*it);
    return it - cont.begin();
  };
  
  template <typename ValueType, typename VertexProperties>
  std::size_t adjlistBC_add_vertex( BC_pooled_vector<ValueType>& cont, const VertexProperties& vp) {
    typedef typename BC_pooled_vector<ValueType>::iterator Iter;
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OEFactory;
    
    if(cont.m_first_hole == BC_hole_desc()) {
      Iter it = cont.m_data.insert(cont.m_data.end(), ValueType(vp));
      ++(cont.m_num_elements);
      OEFactory::create_out_edges(get<ValueType>(*it));
      return it - cont.m_data.begin();
    } else {
      Iter it = cont.m_data.begin() + cont.m_first_hole.value;
      cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole]);
      *it = ValueType(vp);
      ++(cont.m_num_elements);
      OEFactory::create_out_edges(get<ValueType>(*it));
      return it - cont.m_data.begin();
    };
  };
  
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  
  template <typename Container, typename VertexProperties>
  typename Container::iterator adjlistBC_add_vertex(Container& cont, VertexProperties&& vp) {
    typedef typename Container::value_type ValueType;
    return cont.insert(cont.end(), ValueType(std::move(vp)));
  };
  
  template <typename Container, typename VertexProperties>
  typename Container::iterator adjlistBC_add_vertex(Container* cont, VertexProperties&& vp) {
    return adjlistBC_add_vertex(*cont, std::move(vp));
  };
  
  template <typename ValueType, typename VertexProperties>
  std::size_t adjlistBC_add_vertex( ::boost::container::vector<ValueType>& cont, VertexProperties&& vp) {
    return cont.insert(cont.end(), ValueType(std::move(vp)));
  };
  
  template <typename ValueType, typename VertexProperties>
  std::size_t adjlistBC_add_vertex( BC_pooled_vector<ValueType>& cont, VertexProperties&& vp) {
    typedef typename BC_pooled_vector<ValueType>::iterator Iter;
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OEFactory;
    
    if(cont.m_first_hole == BC_hole_desc()) {
      Iter it = cont.m_data.insert(cont.m_data.end(), ValueType(std::move(vp)));
      ++(cont.m_num_elements);
      OEFactory::create_out_edges(get<ValueType>(*it));
      return it - cont.m_data.begin();
    } else {
      Iter it = cont.m_data.begin() + cont.m_first_hole.value;
      cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole]);
      *it = ValueType(std::move(vp));
      ++(cont.m_num_elements);
      OEFactory::create_out_edges(get<ValueType>(*it));
      return it - cont.m_data.begin();
    };
  };
  
#endif
  
  
  
  
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
            typename VertexProperties, typename EdgeProperties>
  struct adjlistBC_vertex_container {
    
    typedef adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::container vertex_container;
    typedef typename vertex_container::size_type vertices_size_type;
    typedef typename VConfig::descriptor vertex_descriptor;
    typedef typename VConfig::stored_type vertex_stored_type;
    typedef typename VConfig::value_type vertex_value_type;
    
    typedef adjlistBC_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> EConfig;
    typedef typename EConfig::container edge_container;
    typedef typename edge_container::size_type edges_size_type;
    typedef typename EConfig::descriptor edge_descriptor;
    typedef typename EConfig::stored_type edge_stored_type;
    typedef typename EConfig::value_type edge_value_type;
    
    mutable vertex_container m_vertices;
    
    adjlistBC_vertex_container() : m_vertices() { };
    
    std::size_t size() const { return m_vertices.size(); };
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    vertex_stored_type& get_stored_vertex(vertex_descriptor v) const { 
      return BC_desc_to_value(m_vertices, v);
    };
    
    edge_stored_type& get_stored_edge(edge_descriptor e) const { 
      return BC_desc_to_value(BC_desc_to_value(m_vertices, e.source).out_edges, e.edge_id);
    };
    
    
    // NOTE: this operation does not invalidate anything.
    // NOTE: This WORKS for ALL vertex container types.
    vertex_descriptor add_vertex(const VertexProperties& vp) {
      return adjlistBC_add_vertex(m_vertices, vp);
    };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_vertex(VertexProperties&& vp) {
      return adjlistBC_add_vertex(m_vertices, std::move(vp));
    };
#endif
    
    // NOTE: this operation only invalidates existing vertex-iterators, 
    // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vecBC).
    // NOTE: This WORKS for ALL vertex container types.
    void clear_vertex(vertex_descriptor v) {
      adjlistBC_clear_vertex<DirectedS>(m_vertices, BC_desc_to_value(m_vertices, v), v);
    };
    
    // NOTE: this operation only invalidates existing vertex-descriptors, 
    // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vecBC).
    // NOTE: This WORKS for ALL vertex container types.
    void remove_vertex(vertex_descriptor v) {
      clear_vertex(v);
      adjlistBC_erase_vertex(m_vertices, v);
    };
    
    // NOTE: this operation does not invalidate anything.
    // NOTE: This WORKS for ALL vertex container types.
    std::pair< edge_descriptor, bool > add_edge(vertex_descriptor u, vertex_descriptor v, const EdgeProperties& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(BC_desc_to_value(m_vertices, u).out_edges, ep, v);
      
      if( raw_result.second ) {
        adjlistBC_add_in_edge(BC_desc_to_value(m_vertices, v), edge_descriptor(u, raw_result.first));
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      } else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< edge_descriptor, bool > add_edge(vertex_descriptor u, vertex_descriptor v, EdgeProperties&& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(BC_desc_to_value(m_vertices, u).out_edges, std::move(ep), v);
      
      if( raw_result.second ) {
        adjlistBC_add_in_edge(BC_desc_to_value(m_vertices, v), edge_descriptor(u, raw_result.first));
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      } else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
#endif
    
    // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
    // NOTE: This WORKS for ALL vertex container types.
    void remove_edge(edge_descriptor e) {
      adjlistBC_erase_in_edge(BC_desc_to_value(m_vertices, get_stored_edge(e).target), e);
      adjlistBC_erase_edge(BC_desc_to_value(m_vertices, e.source).out_edges, e, m_vertices, e.source);
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    void clear() { 
      typedef adjlistBC_out_edges_factory< typename vertex_value_type::edge_container_ptr > OEFactory;
      OEFactory::destroy_all_out_edges(m_vertices);
      m_vertices.clear();
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef adjlistBC_select_vertex_iterator<VertexListS, vertex_container> VIterSelect;
    typedef typename VIterSelect::type vertex_iterator;
    
    std::pair< vertex_iterator, vertex_iterator > vertices() const {
      return VIterSelect::create_range(m_vertices);
    };
    
    
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    typedef adjlistBC_select_out_edge_iterator<OutEdgeListS, edge_container> OEIterSelect;
    typedef typename OEIterSelect::type out_edge_iterator;
    
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor u) const {
      return OEIterSelect::create_range(u, BC_desc_to_value(m_vertices, u).out_edges);
    };
    
    std::pair<edge_descriptor,bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
      std::pair< out_edge_iterator, out_edge_iterator > oe_pair = out_edges(u);
      for(; oe_pair.first != oe_pair.second; ++(oe_pair.first)) {
        if(get_stored_edge(*(oe_pair.first)).target == v)
          return std::pair<edge_descriptor,bool>(*(oe_pair.first), true);
      };
      return std::pair<edge_descriptor,bool>(edge_descriptor(), false);
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    typedef adjlistBC_adjacent_viter<vertex_descriptor, edge_container, out_edge_iterator> adjacency_iterator;
    
    std::pair< adjacency_iterator, adjacency_iterator > adjacent_vertices(vertex_descriptor u) const {
      std::pair< out_edge_iterator, out_edge_iterator > oe_pair = out_edges(u);
      return std::pair< adjacency_iterator, adjacency_iterator >(
        adjacency_iterator(BC_desc_to_value(m_vertices, u).out_edges, oe_pair.first),
        adjacency_iterator(BC_desc_to_value(m_vertices, u).out_edges, oe_pair.second)
      );
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef typename mpl::if_< is_same< DirectedS, directedS >,
      void, 
      typename vertex_stored_type::in_edge_container::iterator >::type in_edge_iterator;
    
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) const {
      return std::pair< in_edge_iterator, in_edge_iterator >(BC_desc_to_value(m_vertices, v).in_edges.begin(), 
                                                             BC_desc_to_value(m_vertices, v).in_edges.end());
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef typename mpl::if_< is_same< DirectedS, directedS >,
      void, 
      adjlistBC_inv_adjacent_viter<vertex_descriptor, in_edge_iterator> >::type inv_adjacency_iterator;
    
    std::pair< inv_adjacency_iterator, inv_adjacency_iterator > inv_adjacent_vertices(vertex_descriptor v) const {
      return std::pair< inv_adjacency_iterator, inv_adjacency_iterator >(
        inv_adjacency_iterator(BC_desc_to_value(m_vertices, v).in_edges.begin()), 
        inv_adjacency_iterator(BC_desc_to_value(m_vertices, v).in_edges.end())
      );
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef adjlistBC_edge_iterator<vertex_container, edge_descriptor> edge_iterator;
    
    std::pair< edge_iterator, edge_iterator > edges() const {
      return std::pair< edge_iterator, edge_iterator >(edge_iterator::begin(m_vertices), edge_iterator::end(m_vertices));
    };
    
  };
  
  
  
  
}; // namespace detail


}; // namespace graph


}; // namespace boost


#endif


















