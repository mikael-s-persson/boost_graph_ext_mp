// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file ltreeBC_containers.hpp
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
  struct ltreeBC_vertex_stored_type;  // forward declaration.
  
  
  // this is fine because boost-containers allow incomplete types:
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltreeBC_vertex_config {
    typedef ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> stored_type;
    typedef typename BC_container_gen<VertexListS, stored_type >::type container;
    typedef typename container::value_type value_type;
    typedef typename BC_select_descriptor< container >::type descriptor;
    
    static descriptor null_vertex() { return BC_null_desc<descriptor>::value(); };
  };
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltreeBC_edge_stored_type {
    typedef ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> self;
    typedef ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::descriptor vertex_descriptor;
    
    vertex_descriptor target;
    mutable EdgeProperties data;
    
    ltreeBC_edge_stored_type(vertex_descriptor aTarget = VConfig::null_vertex()) : target(aTarget), data() { };
    ltreeBC_edge_stored_type(vertex_descriptor aTarget, const EdgeProperties& aData) : target(aTarget), data(aData) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltreeBC_edge_stored_type(vertex_descriptor aTarget, EdgeProperties&& aData) : target(aTarget), data(std::move(aData)) { };
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
  std::size_t hash_value(const ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& ep) {
    return BC_desc_get_hash(ep.target);
  };
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltreeBC_edge_config {
    typedef typename ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>::descriptor vertex_descriptor;
    
    typedef ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> stored_type;
    typedef typename BC_container_gen<OutEdgeListS, stored_type >::type container;
    typedef typename container::value_type value_type;
    typedef typename BC_select_descriptor< container >::type raw_descriptor;
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
  struct ltreeBC_vertex_stored_type {
    typedef ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> self;
    typedef DirectedS directed_tag;
    typedef ltreeBC_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> Config;
    typedef typename Config::container edge_container;
    typedef typename Config::container_ptr edge_container_ptr;
    typedef typename Config::descriptor edge_descriptor;
    typedef edge_descriptor* in_edge_iterator;
    
    VertexProperties data;
    edge_container_ptr out_edges;
    edge_descriptor in_edge;
    
    ltreeBC_vertex_stored_type() : data(), out_edges(), in_edge() { };
    ltreeBC_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges(), in_edge() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltreeBC_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges(), in_edge() { };
#endif
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> {
    typedef ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> self;
    typedef directedS directed_tag;
    typedef ltreeBC_edge_config<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties> Config;
    typedef typename Config::container edge_container;
    typedef typename Config::container_ptr edge_container_ptr;
    typedef typename Config::descriptor edge_descriptor;
    typedef void in_edge_iterator;
    
    VertexProperties data;
    edge_container_ptr out_edges;
    
    ltreeBC_vertex_stored_type() : data(), out_edges() { };
    ltreeBC_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltreeBC_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges() { };
#endif
  };
  
  
  
  //NOTE: ltreeBC_out_edges_factory == adjlistBC_out_edges_factory
  //NOTE: ltreeBC_add_edge = adjlistBC_add_edge
  //NOTE: ltreeBC_find_edge_to = adjlistBC_find_edge_to
  
  
  
  
/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))
  
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_erase_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                               EdgeDesc e) { };
  
  // for vector of edge-desc (in-edges)
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_erase_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                               EdgeDesc e) {
    if(e == vp.in_edge)
      vp.in_edge = EdgeDesc::null_value();
  };
  
  
  // for OutEdgeListS = listS, setS, ...
  template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void ltreeBC_erase_edge(Container& cont, EdgeDesc e, VertexCont& vcont, VertexDesc v) {
    adjlistBC_erase_edge(cont, e, vcont, v);
  };
  template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void ltreeBC_erase_edge(Container* cont, EdgeDesc e, VertexCont& vcont, VertexDesc v) {
    adjlistBC_erase_edge(*cont, e, vcont, v);
  };
  
  // for OutEdgeListS = vecS
  template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
  void ltreeBC_erase_edge( ::boost::container::vector<ValueType>& cont, EdgeDesc e, VertexCont& vertex_cont, VertexDesc v) {
    typedef ::boost::container::vector<ValueType> Container;
    typedef typename Container::iterator Iter;
    using std::swap;
    
    Iter it = BC_desc_to_iterator(cont, e.edge_id);
    Iter it_last = cont.end(); --it_last;
    if(it != it_last) {
      swap(*it, *it_last);
      // If this graph has in-edge references, then they must be updated.
      ltreeBC_update_in_edge_id(BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), 
                                  v, it_last - cont.begin(), it - cont.begin());
    };
    cont.erase(it_last, cont.end());
  };
  
  
  
/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)
  
  
  /* THESE FUNCTIONS SHOULD NOT BE NEEDED
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_add_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) { };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_add_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) {
    vp.in_edge = e;
  };
  */
  
  
  
  
  
/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC    (multi)setBC      unordered_(multi)setBC
// directedS:   O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)    O(log(E/V))       O(1)
  
  //NOTE: ltreeBC_erase_edges_to = adjlistBC_erase_edges_to
  // As long as these special versions of update_in_edge_id are used (by ADL):
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename VertexDesc>
  void adjlistBC_update_in_edge_id(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                                   VertexDesc v, std::size_t old_id, std::size_t new_id) { };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename VertexDesc>
  void adjlistBC_update_in_edge_id(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                                   VertexDesc v, std::size_t old_id, std::size_t new_id) {
    typedef ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VertexValue;
    typedef typename VertexValue::edge_descriptor EdgeDesc;
    
    if(vp.in_edge == EdgeDesc::null_value())
      return;
    if((vp.in_edge.source == v) && (vp.in_edge.edge_id == old_id))
      vp.in_edge.edge_id = new_id;
  };
  
  
/*************************************************************************
 *        functions for clearing the edges of a vertex
 * **********************************************************************/
  // NOTE: The function 'ltreeBC_clear_vertex' works for all graph types.
  
// NOTE: Time complexities: 
//                 vecBC         poolBC         listBC
// directedS:   O(V*(E/V))     O(V*(E/V))     O(V*(E/V))
// bidir:       O((E/V)^2)     O((E/V)^2)     O((E/V)^2)
  
  
  template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v, std::size_t& e_count) {
    adjlistBC_clear_vertex<DirectedS>(cont, vp, v, e_count);
  };
  
  template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v, std::size_t& e_count) {
    typedef typename VertexValue::edge_container OutEdgeCont;
    typedef typename OutEdgeCont::iterator OutEdgeIter;
    typedef typename VertexValue::edge_descriptor EdgeDesc;
    
    // first, remove the in-edge references from the adjacent vertices of v.
    for(OutEdgeIter ei = BC_get_begin_iter(vp.out_edges); ei != BC_get_end_iter(vp.out_edges); ++ei) {
      if( !BC_is_elem_valid(*ei) )
        continue;
      VertexValue& wp = BC_get_value( *BC_desc_to_iterator(cont, BC_get_value(*ei).target) );
      if((wp.in_edge != EdgeDesc::null_value()) && (wp.in_edge.source == v) && 
         (BC_desc_to_iterator(vp.out_edges, wp.in_edge.edge_id) == ei))
        wp.in_edge = EdgeDesc::null_value();
    };
    
    // then, clear the out-going edges.
    e_count -= BC_get_size(vp.out_edges);
    BC_clear_all(vp.out_edges);
    
    // finally, remove the required out-edges of the "parent" vertex of v.
    if(vp.in_edge != EdgeDesc::null_value()) {
      VertexDesc u = vp.in_edge.source;
      VertexValue& up = BC_get_value( *BC_desc_to_iterator(cont, vp.in_edge.source) );
      adjlistBC_erase_edges_to(up.out_edges, cont, u, v, e_count);
      vp.in_edge = EdgeDesc::null_value();
    };
    
  };
  
  
  
  
  
  
/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
  // NOTE: The function 'ltreeBC_erase_vertex' works for all graph types.
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

  
  template <typename Container, typename VertexDesc>
  void ltreeBC_erase_vertex(Container& cont, VertexDesc v) {
    adjlistBC_erase_vertex(cont, v);
  };
  
  
  //NOTE: ltreeBC_update_out_edges_impl = adjlistBC_update_out_edges_impl
  
  
  // O(E)
  template <typename DirectedS, typename ValueType>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    adjlistBC_update_out_edges(cont, old_v_id, new_v_id);
  };
  
  // O(E/V)
  template <typename DirectedS, typename ValueType>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    typedef typename ValueType::edge_container OutEdgeCont;
    typedef typename OutEdgeCont::iterator OEIter;
    typedef typename ValueType::edge_descriptor EdgeDesc;
    typedef typename InEdgeCont::iterator InEdgeIter;
    
    // first, update in-edge vertices
    if(cont[new_v_id].in_edge != EdgeDesc::null_value()) {
      EdgeDesc& e_in = cont[new_v_id].in_edge;
      ValueType& up = cont[e_in.source];
      adjlistBC_update_out_edges_impl(up.out_edges, old_v_id, new_v_id, BC_desc_to_iterator(up.out_edges, e_in.edge_id));
    };
    
    // second, update out-edge vertices
    for(OEIter ei = BC_get_begin_iter(cont[new_v_id].out_edges); ei != BC_get_end_iter(cont[new_v_id].out_edges); ++ei) {
      if(!BC_is_elem_valid(*ei))
        continue;
      ValueType& wp = cont[BC_get_value(*ei).target];
      EdgeDesc& e_in = wp.in_edge;
      if(e_in != EdgeDesc::null_value()) {
        if((e_in.source == old_v_id) && 
           (ei == BC_desc_to_iterator(cont[new_v_id].out_edges, e_in.edge_id))) {
          e_in.source = new_v_id;
          break;
        };
      };
    };
    
  };
  
  template <typename ValueType>
  void ltreeBC_erase_vertex( ::boost::container::vector<ValueType>& cont, std::size_t v) {
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
    ltreeBC_update_out_edges<DirectedS>(cont, old_id, new_id);
  };
  
  template <typename ValueType>
  void ltreeBC_erase_vertex(BC_pooled_vector<ValueType>& cont, std::size_t v) {
    adjlistBC_erase_vertex(cont, v);
  };
  
  
  
  //NOTE: ltreeBC_add_vertex = adjlistBC_add_vertex
  
  
  
  
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
            typename VertexProperties, typename EdgeProperties>
  struct ltreeBC_vertex_container {
    
    typedef ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::container vertex_container;
    typedef typename vertex_container::size_type vertices_size_type;
    typedef typename VConfig::descriptor vertex_descriptor;
    typedef typename VConfig::stored_type vertex_stored_type;
    typedef typename VConfig::value_type vertex_value_type;
    
    typedef ltreeBC_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> EConfig;
    typedef typename EConfig::container edge_container;
    typedef typename edge_container::size_type edges_size_type;
    typedef typename EConfig::descriptor edge_descriptor;
    typedef typename EConfig::stored_type edge_stored_type;
    typedef typename EConfig::value_type edge_value_type;
    
    mutable vertex_container m_vertices;
    std::size_t m_num_edges;
    
    ltreeBC_vertex_container() : m_vertices(), m_num_edges(0) { };
    
    
/********************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE SAME AS ADJ-LIST-BC - NOTE NOTE NOTE *
 * ******************************************************************************/
    
    std::size_t size() const { return m_vertices.size(); };
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    std::size_t num_edges() const { return m_num_edges; };
    
    vertex_stored_type& get_stored_vertex(vertex_descriptor v) const { 
      return BC_get_value(*BC_desc_to_iterator(m_vertices, v));
    };
    
    const edge_stored_type& get_stored_edge(edge_descriptor e) const { 
      return BC_get_value(*BC_desc_to_iterator(get_stored_vertex(e.source).out_edges, e.edge_id));
    };
    
    std::size_t get_out_degree(vertex_descriptor v) const {
      return BC_get_size(get_stored_vertex(v).out_edges);
    };
    
    std::size_t get_in_degree(vertex_descriptor v) const {
      return BC_get_size(get_stored_vertex(v).in_edges);
    };
    
    
    
    
    
/********************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE SAME AS ADJ-LIST-BC - NOTE NOTE NOTE *
 * ******************************************************************************/
    
    
    // NOTE: This WORKS for ALL vertex container types.
    void clear() { 
      typedef adjlistBC_out_edges_factory< typename vertex_value_type::edge_container_ptr > OEFactory;
      OEFactory::destroy_all_out_edges(m_vertices);
      m_vertices.clear();
      m_num_edges = 0;
    };
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef adjlistBC_select_vertex_iterator<VertexListS, vertex_container> VIterSelect;
    typedef typename VIterSelect::type vertex_iterator;
    
    std::pair< vertex_iterator, vertex_iterator > vertices() const {
      return VIterSelect::create_range(m_vertices);
    };
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    typedef adjlistBC_select_out_edge_iterator<OutEdgeListS, edge_container, edge_descriptor> OEIterSelect;
    typedef typename OEIterSelect::type out_edge_iterator;
    
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor u) const {
      return OEIterSelect::create_range(u, get_stored_vertex(u).out_edges);
    };
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    typedef ltreeBC_adjacent_viter<vertex_descriptor, edge_container, out_edge_iterator> adjacency_iterator;
    
    std::pair< adjacency_iterator, adjacency_iterator > adjacent_vertices(vertex_descriptor u) const {
      std::pair< out_edge_iterator, out_edge_iterator > oe_pair = out_edges(u);
      return std::pair< adjacency_iterator, adjacency_iterator >(
        adjacency_iterator(get_stored_vertex(u).out_edges, oe_pair.first),
        adjacency_iterator(get_stored_vertex(u).out_edges, oe_pair.second)
      );
    };
    
    
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef adjlistBC_edge_iterator<vertex_container, edge_descriptor> edge_iterator;
    
    std::pair< edge_iterator, edge_iterator > edges() const {
      return std::pair< edge_iterator, edge_iterator >(edge_iterator::begin(m_vertices), edge_iterator::end(m_vertices));
    };
    
    
    
/***************************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE DIFFERENT FROM ADJ-LIST-BC - NOTE NOTE NOTE *
 * *************************************************************************************/
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef typename vertex_stored_type::in_edge_iterator in_edge_iterator;
    
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) const {
      if(get_stored_vertex(u).in_edge == EConfig::null_edge())
        return std::pair< in_edge_iterator, in_edge_iterator >(NULL,NULL);
      else
        return std::pair< in_edge_iterator, in_edge_iterator >(&(get_stored_vertex(v).in_edge), 
                                                               &(get_stored_vertex(v).in_edge) + 1);
    };
    
    // NOTE: This WORKS for ALL vertex container types.
    typedef typename mpl::if_< is_same< DirectedS, directedS >,
      void, 
      adjlistBC_inv_adjacent_viter<vertex_descriptor, in_edge_iterator> >::type inv_adjacency_iterator;
    
    std::pair< inv_adjacency_iterator, inv_adjacency_iterator > inv_adjacent_vertices(vertex_descriptor v) const {
      std::pair< in_edge_iterator, in_edge_iterator > result_ie = in_edges(v);
      return std::pair< inv_adjacency_iterator, inv_adjacency_iterator >(
        inv_adjacency_iterator(result_ie.first), 
        inv_adjacency_iterator(result_ie.second)
      );
    };
    
    
    
  
  
/*********************************************************
 * NOTE NOTE NOTE - VERIFIED UP TO HERE - NOTE NOTE NOTE *
 * *******************************************************/
  
    
/*************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS WILL DISAPPEAR FOR L-TREE - NOTE NOTE NOTE *
 * ***********************************************************************/
    
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
      ltreeBC_clear_vertex<DirectedS>(m_vertices, get_stored_vertex(v), v, m_num_edges);
    };
    
    // NOTE: this operation only invalidates existing vertex-descriptors, 
    // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vecBC).
    // NOTE: This WORKS for ALL vertex container types.
    void remove_vertex(vertex_descriptor v) {
      clear_vertex(v);
      ltreeBC_erase_vertex(m_vertices, v);
    };
    
    // NOTE: this operation does not invalidate anything.
    // NOTE: This WORKS for ALL vertex container types.
    std::pair< edge_descriptor, bool > add_edge(vertex_descriptor u, vertex_descriptor v, const EdgeProperties& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(get_stored_vertex(u).out_edges, ep, v);
      
      if( raw_result.second ) {
        ++m_num_edges;
        ltreeBC_add_in_edge(get_stored_vertex(v), edge_descriptor(u, raw_result.first));
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      } else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< edge_descriptor, bool > add_edge(vertex_descriptor u, vertex_descriptor v, EdgeProperties&& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(get_stored_vertex(u).out_edges, std::move(ep), v);
      
      if( raw_result.second ) {
        ++m_num_edges;
        ltreeBC_add_in_edge(get_stored_vertex(v), edge_descriptor(u, raw_result.first));
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      } else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
#endif
    
    // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
    // NOTE: This WORKS for ALL vertex container types.
    void remove_edge(edge_descriptor e) {
      ltreeBC_erase_in_edge(get_stored_vertex(get_stored_edge(e).target), e);
      ltreeBC_erase_edge(get_stored_vertex(e.source).out_edges, e, m_vertices, e.source);
      --m_num_edges;
    };
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    std::pair<edge_descriptor,bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      
      std::pair< RawEDesc, bool > raw_result = adjlistBC_find_edge_to(get_stored_vertex(u).out_edges, v);
      
      if( raw_result.second )
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
    
    
    
    
  };
  
  
  
  
}; // namespace detail


}; // namespace graph


}; // namespace boost


#endif


















