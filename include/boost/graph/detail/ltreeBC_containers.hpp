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

#ifndef BOOST_LTREEBC_CONTAINERS_HPP
#define BOOST_LTREEBC_CONTAINERS_HPP

#include <boost/config.hpp>

#include <boost/graph/detail/boost_container_generators.hpp>
#include <boost/graph/detail/adjlistBC_iterators.hpp>
#include <boost/graph/detail/adjlistBC_containers.hpp>

#include <boost/graph/graph_selectors.hpp>              // for directedS, undirectedS, bidirectionalS.

#include <iterator>
#include <utility>
#include <stack>
#include <queue>
#include <vector>

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
  void swap(ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
            ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& rhs) {
    using std::swap;
    swap(lhs.target, rhs.target);
    swap(lhs.data, rhs.data);
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
    void swap(self& rhs) {
      using std::swap;
      swap(data, rhs.data);
      swap(out_edges, rhs.out_edges);
      swap(in_edge, rhs.in_edge);
    };
    self& operator=(self rhs) {
      this->swap(rhs);
      return *this;
    };
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
    typedef int* in_edge_iterator;
    
    VertexProperties data;
    edge_container_ptr out_edges;
    
    ltreeBC_vertex_stored_type() : data(), out_edges() { };
    ltreeBC_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltreeBC_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges() { };
#endif
    void swap(self& rhs) {
      using std::swap;
      swap(data, rhs.data);
      swap(out_edges, rhs.out_edges);
    };
    self& operator=(self rhs) {
      this->swap(rhs);
      return *this;
    };
  };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, 
            typename VertexProperties, typename EdgeProperties>
  void swap(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
            ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& rhs) {
    lhs.swap(rhs);
  };
  
  
  
  //NOTE: ltreeBC_out_edges_factory == adjlistBC_out_edges_factory
  //NOTE: ltreeBC_add_edge = adjlistBC_add_edge
  //NOTE: ltreeBC_find_edge_to = adjlistBC_find_edge_to
  //NOTE: ltreeBC_add_vertex = adjlistBC_add_vertex
  
  
  
  
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
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))
  
  
  
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
      adjlistBC_update_in_edge_id(BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), 
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
  
  
  template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_add_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) { };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
  void ltreeBC_add_in_edge(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
                             EdgeDesc e) {
    vp.in_edge = e;
  };
  
  
  
/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
  // NOTE: The function 'ltreeBC_erase_vertex' works for all graph types.
// NOTE: Time complexities: 
//               vecBC      poolBC      listBC
// directedS:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

  
  
  template <typename Container, typename VDescContainer>
  void ltreeBC_erase_vertices(Container& cont, VDescContainer& v_list) {
    for(typename VDescContainer::iterator it = v_list.begin(); it != v_list.end(); ++it)
      adjlistBC_erase_vertex(cont, *it);
  };
  
  
  //NOTE: ltreeBC_update_out_edges_impl = adjlistBC_update_out_edges_impl
  
  
  // O(E)
  template <typename DirectedS, typename ValueType>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    adjlistBC_update_out_edges<DirectedS>(cont, old_v_id, new_v_id);
  };
  
  // O(E/V)
  template <typename DirectedS, typename ValueType>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont, 
                                          std::size_t old_v_id, std::size_t new_v_id) {
    typedef typename ValueType::edge_container OutEdgeCont;
    typedef typename OutEdgeCont::iterator OEIter;
    typedef typename ValueType::edge_descriptor EdgeDesc;
    // first, update in-edge vertices
    if(cont[new_v_id].in_edge != EdgeDesc::null_value()) {
      EdgeDesc& e_in = cont[new_v_id].in_edge;
      ValueType& up = cont[e_in.source];
      // std::cout << "address 3: " << cont[new_v_id].data << " " << up.data << " " << e_in.edge_id << " " << up.out_edges.size() << std::endl;
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
  
  
  
  // O(1)
  template <typename DirectedS, typename ValueType>
  typename enable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_invalidate_in_edges(::boost::container::vector<ValueType>&, std::size_t) { };
  
  // O(1)
  template <typename DirectedS, typename ValueType>
  typename disable_if< is_same< DirectedS, directedS >,
  void >::type ltreeBC_invalidate_in_edges(::boost::container::vector<ValueType>& cont, std::size_t u) {
    typedef typename ValueType::edge_container OutEdgeCont;
    typedef typename OutEdgeCont::iterator OEIter;
    typedef typename ValueType::edge_descriptor EdgeDesc;
    
    for(OEIter ei = BC_get_begin_iter(cont[u].out_edges); ei != BC_get_end_iter(cont[u].out_edges); ++ei) {
      if(!BC_is_elem_valid(*ei))
        continue;
      cont[BC_get_value(*ei).target].in_edge = EdgeDesc::null_value();
    };
  };
  
  
  
  template <typename ValueType>
  void ltreeBC_erase_vertices( ::boost::container::vector<ValueType>& cont, std::vector< std::size_t >& v_list) {
    typedef ::boost::container::vector<ValueType> Container;
    typedef typename Container::iterator Iter;
    typedef typename std::vector< std::size_t >::iterator VIter;
    typedef adjlistBC_out_edges_factory< typename ValueType::edge_container_ptr > OutEdgeFactory;
    typedef typename ValueType::directed_tag DirectedS;
    using std::swap;
    
    for(VIter v_it = v_list.begin(); v_it != v_list.end(); ++v_it) {
      // NOTE: First, invalidate the in-edge references of the children, then destroy out-edge list.
      ltreeBC_invalidate_in_edges<DirectedS>(cont, *v_it);
      Iter it = BC_desc_to_iterator(cont, *v_it);
      OutEdgeFactory::destroy_out_edges(*it);
      Iter it_last = cont.end(); --it_last;
      if(it == it_last) {
        cont.erase(it_last);
        continue;
      };
      swap(*it, *it_last);
      std::size_t old_id = it_last - cont.begin();
      std::size_t new_id = it - cont.begin();
      cont.erase(it_last); 
      ltreeBC_update_out_edges<DirectedS>(cont, old_id, new_id);
      for(VIter v_it2 = v_it + 1; v_it2 != v_list.end(); ++v_it2)
        if(*v_it2 == old_id)
          *v_it2 = new_id;
    };
  };
  
  template <typename ValueType>
  void ltreeBC_erase_vertices(BC_pooled_vector<ValueType>& cont, std::vector< std::size_t >& v_list) {
    for(typename std::vector< std::size_t >::iterator it = v_list.begin(); it != v_list.end(); ++it)
      adjlistBC_erase_vertex(cont, *it);
  };
  
  
  
  
  
  // O(E)
  template <typename DirectedS, typename Container, typename VertexValue, typename Vertex>
  typename enable_if< is_same< DirectedS, directedS >,
  VertexValue >::type::edge_descriptor ltreeBC_get_in_edge(Container& vcont, VertexValue& vp, Vertex v) {
    typedef typename Container::iterator VIter;
    typedef typename VertexValue::edge_container EdgeCont;
    typedef typename EdgeCont::iterator OEIter;
    typedef typename VertexValue::edge_descriptor EdgeDesc;
    
    for(VIter it = vcont.begin(); it != vcont.end(); ++it) {
      if(BC_is_elem_valid(*it)) {
        VertexValue& up = BC_get_value(*it);
        for(OEIter ei = BC_get_begin_iter(up.out_edges); ei != BC_get_end_iter(up.out_edges); ++ei) 
          if(BC_is_elem_valid(*ei) && (BC_get_value(*ei).target == v)) 
            return EdgeDesc(BC_iterator_to_desc(vcont, it), BC_iterator_to_desc(up.out_edges, ei));
      };
    };
    
    return EdgeDesc::null_value();
  };
  
  // O(1)
  template <typename DirectedS, typename Container, typename VertexValue, typename Vertex>
  typename disable_if< is_same< DirectedS, directedS >,
  VertexValue >::type::edge_descriptor ltreeBC_get_in_edge(Container& vcont, VertexValue& vp, Vertex v) {
    return vp.in_edge;
  };
  
  
  
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
    vertex_descriptor m_root;
    
    ltreeBC_vertex_container() : m_vertices(), m_root(BC_null_desc<vertex_descriptor>::value()) { };
    ~ltreeBC_vertex_container() { clear(); };
    
  private:
    ltreeBC_vertex_container(const ltreeBC_vertex_container&);
    ltreeBC_vertex_container& operator=(const ltreeBC_vertex_container&);
  public:
    
    void swap(ltreeBC_vertex_container& rhs) {
      using std::swap;
      m_vertices.swap(rhs.m_vertices);
      // swap(m_vertices, rhs.m_vertices);
      swap(m_root, rhs.m_root);
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltreeBC_vertex_container(ltreeBC_vertex_container&& rhs) : m_vertices(), m_root(BC_null_desc<vertex_descriptor>::value()) {
      swap(rhs);
    };
    ltreeBC_vertex_container& operator=(ltreeBC_vertex_container&& rhs) {
      swap(rhs);
      return *this;
    };
#endif
    
/********************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE SAME AS ADJ-LIST-BC - NOTE NOTE NOTE *
 * ******************************************************************************/
    
    std::size_t size() const { return m_vertices.size(); };
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    std::size_t num_edges() const { return (m_vertices.size() == 0 ? 0 : m_vertices.size() - 1); };
    
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
      if( (v != BC_null_desc<vertex_descriptor>::value()) &&
          (v != m_root) &&
          (BC_is_elem_valid(*BC_desc_to_iterator(m_vertices, v))) )
        return 1;
      else
        return 0;
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    void clear() { 
      typedef adjlistBC_out_edges_factory< typename vertex_value_type::edge_container_ptr > OEFactory;
      OEFactory::destroy_all_out_edges(m_vertices);
      m_vertices.clear();
      m_root = BC_null_desc<vertex_descriptor>::value();
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
      if(get_stored_vertex(v).in_edge == EConfig::null_edge())
        return std::pair< in_edge_iterator, in_edge_iterator >(NULL,NULL);
      else
        return std::pair< in_edge_iterator, in_edge_iterator >(&(get_stored_vertex(v).in_edge), 
                                                               &(get_stored_vertex(v).in_edge) + 1);
    };
    
    
    // NOTE: This WORKS for ALL vertex container types.
    vertex_descriptor get_parent(vertex_descriptor v) const {
      return ltreeBC_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v).source;
    };
    
    
    
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
    std::size_t get_depth(const vertex_descriptor& u) const {
      typedef typename edge_container::iterator OEIter;
      typedef std::pair< std::size_t, vertex_descriptor > TaskType;
      
      std::size_t max_depth = 0;
      std::stack< TaskType > tasks;
      tasks.push( TaskType(0, u) );
      while( ! tasks.empty() ) {
        TaskType cur = tasks.top(); tasks.pop();
        ++(cur.first);
        if(cur.first > max_depth)
          max_depth = cur.first;
        vertex_stored_type& vp = get_stored_vertex(cur.second);
        for(OEIter ei = BC_get_begin_iter(vp.out_edges); ei != BC_get_end_iter(vp.out_edges); ++ei)
          tasks.push( TaskType(cur.first, BC_get_value(*ei).target) );
      };
      return max_depth;
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
    
    
    
    // NOTE: This WORKS for ALL vertex container types.
    // NOTE: This WORKS for ALL edge container types.
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_new_vertex(const VertexProperties& vp) {
      return adjlistBC_add_vertex(m_vertices, vp);
    };
    void add_root_vertex(const VertexProperties& vp) { m_root = add_new_vertex(vp); };
#else
    template <typename VP>
    vertex_descriptor add_new_vertex(VP&& vp) {
      return adjlistBC_add_vertex(m_vertices, std::forward<VP>(vp));
    };
    template <typename VP>
    void add_root_vertex(VP&& vp) { m_root = add_new_vertex(std::forward<VP>(vp)); };
#endif
    
    
    // NOTE: this operation does not invalidate anything.
    // NOTE: This WORKS for ALL vertex container types.
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< edge_descriptor, bool > add_new_edge(vertex_descriptor u, vertex_descriptor v, const EdgeProperties& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(get_stored_vertex(u).out_edges, ep, v);
#else
    template <typename EP>
    std::pair< edge_descriptor, bool > add_new_edge(vertex_descriptor u, vertex_descriptor v, EP&& ep) {
      typedef typename edge_descriptor::edge_id_type RawEDesc;
      std::pair< RawEDesc, bool > raw_result = adjlistBC_add_edge(get_stored_vertex(u).out_edges, std::forward<EP>(ep), v);
#endif
      if( raw_result.second ) {
        ltreeBC_add_in_edge(get_stored_vertex(v), edge_descriptor(u, raw_result.first));
        return std::pair< edge_descriptor, bool >(edge_descriptor(u, raw_result.first), true);
      } else 
        return std::pair< edge_descriptor, bool >(edge_descriptor(), false);
    };
    
    
    // NOTE: this operation does not invalidate anything.
    // NOTE: This WORKS for ALL vertex container types.
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, const VertexProperties& vp, const EdgeProperties& ep) {
      vertex_descriptor new_node = add_new_vertex(vp);
      std::pair< edge_descriptor, bool > new_edge = add_new_edge(v, new_node, ep);
      return std::make_pair(new_node, new_edge.first);
    };
#else
    template <typename VP, typename EP>
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, VP&& vp, EP&& ep) {
      vertex_descriptor new_node = add_new_vertex(std::forward<VP>(vp));
      std::pair< edge_descriptor, bool > new_edge = add_new_edge(v, new_node, std::forward<EP>(ep));
      return std::make_pair(new_node, new_edge.first);
    };
#endif
    
    
    
    
    template <typename Vertex_OIter, typename Edge_OIter>
    std::pair<Vertex_OIter, Edge_OIter> clear_children_impl(vertex_descriptor v, Vertex_OIter vit_out, Edge_OIter eit_out) {
      typedef typename edge_container::iterator OEIter;
      
      std::vector< vertex_descriptor > death_row;
      std::queue< vertex_descriptor > bft_queue;
      bft_queue.push(v);
      // Put all children on death-row:
      while( ! bft_queue.empty() ) {
        vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
        bft_queue.pop(); 
        for(OEIter ei = BC_get_begin_iter(v_value.out_edges); ei != BC_get_end_iter(v_value.out_edges); ++ei) {
          if(!BC_is_elem_valid(*ei))
            continue;
          death_row.push_back(BC_get_value(*ei).target);
          bft_queue.push(BC_get_value(*ei).target);
          vertex_stored_type& u_value = get_stored_vertex(BC_get_value(*ei).target);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          *(vit_out++) = std::move(u_value.data);
          *(eit_out++) = std::move(BC_get_value(*ei).data);
#else
          *(vit_out++) = u_value.data;
          *(eit_out++) = BC_get_value(*ei).data;
#endif
        };
      };
      
      // Check if we removed the root:
      if(v == m_root) { // v is the root vertex.
        clear();
      } else {
        // remove the out-edges.
        BC_clear_all(get_stored_vertex(v).out_edges);
        // Execute the death-row vertices!
        ltreeBC_erase_vertices(m_vertices, death_row);
      };
      
      return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
    };
    
    template <typename OutputIter>
    OutputIter clear_children_impl(vertex_descriptor v, OutputIter it_out) {
      return clear_children_impl(v, it_out, ignore_output_iter()).first;
    };
    
    void clear_children_impl(vertex_descriptor v) { clear_children_impl(v, ignore_output_iter(), ignore_output_iter()); };
    
    
    
    
    template <typename Vertex_OIter, typename Edge_OIter>
    std::pair<Vertex_OIter, Edge_OIter> remove_branch_impl(vertex_descriptor v, Vertex_OIter vit_out, Edge_OIter eit_out) {
      typedef typename edge_container::iterator OEIter;
      
      edge_descriptor p_edge = ltreeBC_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v);
      
      if(p_edge.source == VConfig::null_vertex()) {
        *(eit_out++) = edge_value_type();
      } else {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(eit_out++) = std::move(get_stored_edge(p_edge).data);
#else
        *(eit_out++) = get_stored_edge(p_edge).data;
#endif
      };
      
      std::vector< vertex_descriptor > death_row;
      death_row.push_back(v);
      std::queue< vertex_descriptor > bft_queue;
      bft_queue.push(v);
      // Put all children on death-row:
      while( ! bft_queue.empty() ) {
        vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(vit_out++) = std::move(v_value.data);
#else
        *(vit_out++) = v_value.data;
#endif
        bft_queue.pop(); 
        for(OEIter ei = BC_get_begin_iter(v_value.out_edges); ei != BC_get_end_iter(v_value.out_edges); ++ei) {
          if(!BC_is_elem_valid(*ei))
            continue;
          death_row.push_back(BC_get_value(*ei).target);
          bft_queue.push(BC_get_value(*ei).target);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
          *(eit_out++) = std::move(BC_get_value(*ei).data);
#else
          *(eit_out++) = BC_get_value(*ei).data;
#endif
        };
      };
      
      // Check if we removed the root:
      if(p_edge.source == VConfig::null_vertex()) {
        // v must be the root vertex.
        clear();
        return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
      } else {
        // remove the edge.
        ltreeBC_erase_edge(get_stored_vertex(p_edge.source).out_edges, p_edge, m_vertices, p_edge.source);
        
        // Execute the death-row vertices!
        ltreeBC_erase_vertices(m_vertices, death_row);
        
        return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
      };
      
    };
    
    template <typename OutputIter>
    OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      return remove_branch_impl(v, it_out, ignore_output_iter()).first;
    };
    
    void remove_branch_impl(vertex_descriptor v) { remove_branch_impl(v, ignore_output_iter(), ignore_output_iter()); };
    
    
  };
  
  template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
            typename VertexProperties, typename EdgeProperties>
  void swap(ltreeBC_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
            ltreeBC_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& rhs) {
    lhs.swap(rhs);
  };
  
  
  
}; // namespace detail


}; // namespace graph


}; // namespace boost


#endif


















