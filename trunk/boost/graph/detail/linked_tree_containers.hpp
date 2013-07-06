// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file linked_tree_containers.hpp
 * 
 * This library implements (doubly-)linked tree data structures that underpin the linked_tree template. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_LINKED_TREE_CONTAINERS_HPP
#define BOOST_LINKED_TREE_CONTAINERS_HPP

#include <boost/config.hpp>

#include <boost/variant.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <boost/graph/detail/linked_tree_iterators.hpp>

#include <stack>
#include <queue>
#include <iterator>
#include <utility>

#include <boost/container/vector.hpp>
#include <boost/container/list.hpp>

#if !defined BOOST_NO_SLIST
#include <boost/container/slist.hpp>
#endif

#include <boost/container/set.hpp>
#include <boost/unordered_set.hpp>


namespace boost {


  template <class Selector, class ValueType>
  struct boost_container_gen { };

  template <class ValueType>
  struct boost_container_gen<listS, ValueType> {
    typedef ::boost::container::list<ValueType> type;
  };
#if !defined BOOST_NO_SLIST
  template <class ValueType>
  struct boost_container_gen<slistS, ValueType> {
    typedef ::boost::container::slist<ValueType> type;
  };
#endif
  template <class ValueType>
  struct boost_container_gen<vecS, ValueType> {
    typedef ::boost::container::vector<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<mapS, ValueType> {
    typedef ::boost::container::set<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<setS, ValueType> {
    typedef ::boost::container::set<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<multisetS, ValueType> {
    typedef ::boost::container::multiset<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<multimapS, ValueType> {
    typedef ::boost::container::multiset<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<hash_setS, ValueType> {
    typedef ::boost::unordered_set<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<hash_mapS, ValueType> {
    typedef ::boost::unordered_set<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<hash_multisetS, ValueType> {
    typedef ::boost::unordered_multiset<ValueType> type;
  };

  template <class ValueType>
  struct boost_container_gen<hash_multimapS, ValueType> {
    typedef ::boost::unordered_multiset<ValueType> type;
  };


namespace graph { 

namespace detail {
  
  
  template <typename DescType>
  struct ltree_null_desc {
    static DescType value() { return DescType(); };
  };
  
  template <typename DescType>
  struct ltree_null_desc<DescType*> {
    static DescType* value() { return NULL; };
  };
  
  template <>
  struct ltree_null_desc<std::size_t> {
    static std::size_t value() { return (std::numeric_limits<std::size_t>::max)(); };
  };
  
  
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  struct ltree_edge_desc {
    typedef typename ::boost::detail::is_random_access<OutEdgeListS>::type is_rand_access;
    typedef typename boost_container_gen<OutEdgeListS, EdgeValueType>::type edge_container;
    typedef typename edge_container::iterator edge_iter;
    typedef typename ::boost::mpl::if_< is_rand_access, std::size_t, edge_iter >::type edge_id_type;
    typedef Vertex source_descriptor;
    
    Vertex source;
    edge_id_type edge_id;
    
    ltree_edge_desc(Vertex aSrc = ltree_null_desc<Vertex>::value(), 
                    edge_id_type aEdgeId = ltree_null_desc<edge_id_type>::value()) : 
                    source(aSrc), edge_id(aEdgeId) { };
    
    static const ltree_edge_desc null_value;
  };
  
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS> ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>::null_value = ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>();
  
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator==(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return (lhs.source == rhs.source) && (lhs.edge_id == rhs.edge_id);
  };
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator!=(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return (lhs.source != rhs.source) || (lhs.edge_id != rhs.edge_id);
  };
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator >(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return (lhs.source > rhs.source) || ((lhs.source == rhs.source) && (lhs.edge_id > rhs.edge_id)); 
  };
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator >=(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return (lhs.source > rhs.source) || ((lhs.source == rhs.source) && (lhs.edge_id >= rhs.edge_id)); 
  };
  
  // NOTE: the extra set of parentheses are there to prevent the compiler for assuming that the left of the < is the start of a template.
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator <(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return ((lhs.source) < (rhs.source)) || (((lhs.source) == (rhs.source)) && ((lhs.edge_id) < (rhs.edge_id))); 
  };
  template <typename Vertex, typename EdgeValueType, typename OutEdgeListS>
  bool operator <=(const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& lhs, const ltree_edge_desc<Vertex, EdgeValueType, OutEdgeListS>& rhs) { 
    return ((lhs.source) < (rhs.source)) || (((lhs.source) == (rhs.source)) && ((lhs.edge_id) <= (rhs.edge_id))); 
  };
  
  
  
  
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_vertex_stored_type;  // forward declaration.
  
  
  
  // this is fine because boost-containers allow incomplete types:
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_vertex_config {
    typedef ltree_vertex_stored_type<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> value_type;
    typedef typename boost_container_gen<VertexListS, value_type >::type container;
    typedef typename container::iterator descriptor;
    
    static descriptor null_vertex() { return ltree_null_desc<descriptor>::value(); };
  };
  
  
  struct ltree_hole_desc {
    std::size_t value;
    explicit ltree_hole_desc(std::size_t aValue = ltree_null_desc<std::size_t>::value()) : value(aValue) { };
  };
  
  template <typename OutEdgeListS, typename VertexProperties, typename EdgeProperties>
  struct ltree_vertex_config<vecS, OutEdgeListS, VertexProperties, EdgeProperties> {
    typedef ltree_vertex_stored_type<vecS, OutEdgeListS, VertexProperties, EdgeProperties> stored_type;
    typedef variant< stored_type, ltree_hole_desc > value_type;
    typedef typename boost_container_gen<vecS, value_type >::type container;
    typedef std::size_t descriptor;
    typedef ltree_hole_desc hole_descriptor;
    
    static descriptor null_vertex() { return ltree_null_desc<descriptor>::value(); };
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_edge_stored_type {
    typedef ltree_vertex_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::descriptor vertex_descriptor;
    
    vertex_descriptor target;
    EdgeProperties data;
    
    ltree_edge_stored_type(vertex_descriptor aTarget = VConfig::null_vertex(), const EdgeProperties& aData = EdgeProperties()) : target(aTarget), data(aData) { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltree_edge_stored_type(vertex_descriptor aTarget, EdgeProperties&& aData) : target(aTarget), data(std::move(aData)) { };
#endif
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_edge_config {
    typedef typename ltree_vertex_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties>::descriptor vertex_descriptor;
    
    typedef ltree_edge_stored_type<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> value_type;
    typedef typename boost_container_gen<OutEdgeListS, value_type >::type container;
    typedef ltree_edge_desc<vertex_descriptor, value_type, OutEdgeListS> descriptor;
    
    static descriptor null_edge() { return ltree_null_desc<descriptor>::value(); };
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_vertex_stored_type {
    typedef typename ltree_edge_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties>::container edge_container;
    typedef typename ltree_edge_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties>::descriptor edge_descriptor;
    
    VertexProperties data;
    edge_container out_edges;
    edge_descriptor in_edge;
    
    ltree_vertex_stored_type(const VertexProperties& aData) : data(aData), out_edges(), in_edge() { };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ltree_vertex_stored_type(VertexProperties&& aData) : data(std::move(aData)), out_edges(), in_edge() { };
#endif
  };
  
  
  
  
  
  
  template <typename Container>
  typename Container::iterator desc_to_iterator(Container& c, std::size_t d) {
    return c.begin() + d;
  };
  
  template <typename Container>
  typename Container::const_iterator desc_to_iterator(const Container& c, std::size_t d) {
    return c.begin() + d;
  };
  
  template <typename Container, typename Iter>
  Iter desc_to_iterator(const Container&, Iter it) { return it; };
  
  template <typename IsRandAccess, typename Container>
  typename enable_if< IsRandAccess,
  std::size_t >::type iterator_to_desc(Container& c, typename Container::iterator it) {
    return it - c.begin();
  };
  
  template <typename IsRandAccess, typename Container>
  typename enable_if< IsRandAccess,
  std::size_t >::type iterator_to_desc(const Container& c, typename Container::const_iterator it) {
    return it - c.begin();
  };
  
  template <typename IsRandAccess, typename Container, typename Iter>
  typename disable_if< IsRandAccess,
  Iter >::type iterator_to_desc(const Container&, Iter it) { return it; };
  
  
  /* Dummy "ignore" output iterator that is used in remove-branch functions. */
  struct ignore_output_iter {
    struct value_type {
      template <typename T>
      value_type& operator=(const T&) { return *this; };
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      template <typename T>
      value_type& operator=(T&&) { return *this; };
#endif
    };
    typedef std::ptrdiff_t difference_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::random_access_iterator_tag iterator_category;
    
    ignore_output_iter& operator++() { return *this; };
    ignore_output_iter& operator++(int) { return *this; };
    ignore_output_iter& operator--() { return *this; };
    ignore_output_iter& operator--(int) { return *this; };
    
    value_type operator*() const { return value_type(); };
  };
  
  
  
  template <typename Container>
  struct ltree_vertex_is_not_hole {
    const Container* p_cont;
    explicit ltree_vertex_is_not_hole(const Container* aPCont = NULL) : p_cont(aPCont) { };
    template <typename Desc>
    bool operator()(Desc d) {
      return ((*p_cont)[d].which() == 0);
    };
  };
  
  
  template <typename Vertex>
  struct ltree_edge_has_a_source {
    template <typename Desc>
    bool operator()(Desc d) {
      return (d.source != ltree_null_desc<Vertex>::value());
    };
  };
  
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_pooled_vector_container {
    
    typedef typename ::boost::detail::is_random_access<VertexListS>::type vertex_rand_access;
    typedef typename ::boost::detail::is_random_access<OutEdgeListS>::type edge_rand_access;
    
    typedef ltree_vertex_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::container vertex_container;
    typedef typename vertex_container::size_type vertices_size_type;
    typedef typename VConfig::descriptor vertex_descriptor;
    typedef typename VConfig::stored_type vertex_stored_type;
    typedef typename VConfig::value_type vertex_value_type;
    typedef ltree_hole_desc hole_descriptor;
    
    typedef ltree_edge_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> EConfig;
    typedef typename EConfig::container edge_container;
    typedef typename edge_container::size_type edges_size_type;
    typedef typename EConfig::descriptor edge_descriptor;
    typedef typename EConfig::value_type edge_value_type;
    
    typedef typename vertex_container::iterator v_iter;
    typedef typename edge_container::iterator e_iter;
    
    
    typedef filter_iterator< ltree_vertex_is_not_hole<vertex_container>, ltree_viter_from_index > vertex_iterator;
    typedef ltree_child_viter_from_iter<vertex_descriptor, e_iter> child_vertex_iterator;
    typedef child_vertex_iterator adjacency_iterator;
    
    typedef typename mpl::if_< edge_rand_access,
      ltree_oeiter_from_index< edge_descriptor >,
      ltree_oeiter_from_iter< edge_descriptor > >::type out_edge_iterator;
    typedef const edge_descriptor* in_edge_iterator;
    
    typedef ltree_eiter_from_pooliter< edge_descriptor, typename vertex_container::const_iterator > edge_iter_impl;
    typedef filter_iterator< ltree_edge_has_a_source<vertex_descriptor>, edge_iter_impl > edge_iterator;
    
    
    mutable vertex_container m_vertices;
    vertex_descriptor m_root;
    hole_descriptor m_first_hole;
    vertices_size_type m_num_vertices;
    
    ltree_pooled_vector_container() : m_vertices(), m_root(VConfig::null_vertex()), m_first_hole(VConfig::null_vertex()), m_num_vertices(0) { };
    
    std::size_t size() const { return m_num_vertices; };
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    vertex_stored_type& get_stored_vertex(vertex_descriptor v) const { return get<vertex_stored_type>(m_vertices[v]); };
    
    edge_value_type& get_stored_edge(edge_descriptor e) { 
      return *desc_to_iterator(get<vertex_stored_type>(m_vertices[e.source]).out_edges, e.edge_id);
    };
    const edge_value_type& get_stored_edge(edge_descriptor e) const { 
      return *desc_to_iterator(get<vertex_stored_type>(m_vertices[e.source]).out_edges, e.edge_id);
    };
    
    
    vertex_descriptor add_new_vertex(const VertexProperties& vp) {
      v_iter result;
      if(m_first_hole.value == VConfig::null_vertex()) {
        result = m_vertices.insert(m_vertices.end(), vertex_value_type(vertex_stored_type(vp)));
      } else {
        result = desc_to_iterator(m_vertices, m_first_hole.value);
        m_first_hole = get< hole_descriptor >(*result);
        *result = vertex_value_type(vertex_stored_type(vp));
      };
      ++m_num_vertices;
      return iterator_to_desc<vertex_rand_access>(m_vertices, result);
    };
    void add_root_vertex(const VertexProperties& vp) { m_root = add_new_vertex(vp); };
    
    
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, const VertexProperties& vp, const EdgeProperties& ep) {
      vertex_descriptor new_node = add_new_vertex(vp);
      edge_container& v_oe = get_stored_vertex(v).out_edges;
      v_oe.insert(v_oe.end(), edge_value_type(new_node, ep));
      // update the edge-descriptors in the children nodes.
      for(e_iter ei = v_oe.begin(); ei != v_oe.end(); ++ei) {
        edge_descriptor& tmp_e = get_stored_vertex(ei->target).in_edge;
        tmp_e.source = v;
        tmp_e.edge_id = iterator_to_desc<edge_rand_access>(v_oe, ei);
      };
      return std::make_pair(new_node, get_stored_vertex(new_node).in_edge);
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, VertexProperties&& vp, EdgeProperties&& ep) {
      vertex_descriptor new_node = add_new_vertex(std::move(vp));
      edge_container& v_oe = get_stored_vertex(v).out_edges;
      v_oe.insert(v_oe.end(), edge_value_type(new_node, std::move(ep)));
      // update the edge-descriptors in the children nodes.
      for(e_iter ei = v_oe.begin(); ei != v_oe.end(); ++ei) {
        edge_descriptor& tmp_e = get_stored_vertex(ei->target).in_edge;
        tmp_e.source = v;
        tmp_e.edge_id = iterator_to_desc<edge_rand_access>(v_oe, ei);
      };
      return std::make_pair(new_node, get_stored_vertex(new_node).in_edge);
    };
#endif
    
    
    template <typename OutputIter>
    OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      vertex_descriptor parent = get_stored_vertex(v).in_edge.source;
      std::stack< vertex_descriptor > death_row;
      death_row.push(v);
      std::queue< vertex_descriptor > bft_queue;
      bft_queue.push(v);
      // Put all children on death-row:
      while( ! bft_queue.empty() ) {
        vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(it_out++) = std::move(v_value.data);
#else
        *(it_out++) = v_value.data;
#endif
        bft_queue.pop(); 
        for(e_iter ei = v_value.out_edges.begin(); ei != v_value.out_edges.end(); ++ei) {
          death_row.push(ei->target);
          bft_queue.push(ei->target);
        };
      };
      
      // Check if we removed the (or a) root:
      e_iter ei;
      if(parent == VConfig::null_vertex()) {
        if(v == m_root) {
          clear();
          return it_out;
        };
      } else {
        // find the edge to be removed:
        edge_container& p_oe = get_stored_vertex(parent).out_edges;
        ei = p_oe.begin();
        while( (ei != p_oe.end()) && (ei->target != v) )
          ++ei;
      };
      
      // Execute them!
      while( ! death_row.empty() ) {
        v_iter tmp_it = desc_to_iterator(m_vertices, death_row.top()); 
        death_row.pop(); 
        *tmp_it = vertex_value_type(m_first_hole);
        m_first_hole.value = iterator_to_desc<vertex_rand_access>(m_vertices, tmp_it);
        --m_num_vertices;
      };
      
      if(parent == VConfig::null_vertex())
        return it_out;
      
      // remove the edge.
      edge_container& p_oe = get_stored_vertex(parent).out_edges;
      p_oe.erase(ei);
      // update the edge-descriptors in the children nodes.
      for(ei = p_oe.begin(); ei != p_oe.end(); ++ei)
        get_stored_vertex(ei->target).in_edge.edge_id = iterator_to_desc<edge_rand_access>(p_oe, ei);
      
      return it_out;
    };
    void remove_branch_impl(vertex_descriptor v) { remove_branch_impl(v, ignore_output_iter()); };
    
    void clear() { 
      m_vertices.clear();
      m_root = VConfig::null_vertex();
      m_first_hole.value = VConfig::null_vertex();
      m_num_vertices = 0;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_new_vertex(VertexProperties&& vp) {
      v_iter result;
      if(m_first_hole.value == VConfig::null_vertex()) {
        result = m_vertices.insert(m_vertices.end(), vertex_value_type(vertex_stored_type(std::move(vp))));
      } else {
        result = desc_to_iterator(m_vertices, m_first_hole.value);
        m_first_hole = get< hole_descriptor >(*result);
        *result = vertex_value_type(vertex_stored_type(std::move(vp)));
      };
      ++m_num_vertices;
      return iterator_to_desc<vertex_rand_access>(m_vertices, result);
    };
    void add_root_vertex(VertexProperties&& vp) { m_root = add_new_vertex(std::move(vp)); };
#endif
    
    
    std::size_t get_depth(const vertex_descriptor& u) const {
      const edge_container& u_oe = get_stored_vertex(u).out_edges;
      std::size_t max_depth = 0;
      typedef typename edge_container::const_iterator EdgeIter;
      for(EdgeIter ei = u_oe.begin(); ei != u_oe.end(); ++ei) {
        std::size_t temp = get_depth( ei->target );
        if(temp > max_depth)
          max_depth = temp;
      };
      return max_depth + 1;
    };
    
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor v) const {
      edge_container& v_oe = get_stored_vertex(v).out_edges;
      return std::make_pair(
        out_edge_iterator(
          edge_descriptor(
            v, 
            iterator_to_desc<edge_rand_access>(
              v_oe, 
              v_oe.begin()
            )
          )
        ),
        out_edge_iterator(
          edge_descriptor(
            v, 
            iterator_to_desc<edge_rand_access>(
              v_oe, 
              v_oe.end()
            )
          )
        )
      );
    };
    
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) const {
      const edge_descriptor& tmp_e = get_stored_vertex(v).in_edge;
      if(tmp_e.source == VConfig::null_vertex())
        return std::make_pair(&tmp_e, &tmp_e);
      else
        return std::make_pair(&tmp_e, &tmp_e + 1);
    };
    
    std::pair< vertex_iterator, vertex_iterator > vertices() const {
      ltree_viter_from_index v_beg = ltree_viter_from_index::begin(m_vertices);
      ltree_viter_from_index v_end = ltree_viter_from_index::end(m_vertices);
      return std::pair< vertex_iterator, vertex_iterator >(
        vertex_iterator(ltree_vertex_is_not_hole<vertex_container>(&m_vertices), v_beg, v_end),
        vertex_iterator(ltree_vertex_is_not_hole<vertex_container>(&m_vertices), v_end, v_end)
      );
    };
    
    std::pair< edge_iterator, edge_iterator > edges() const {
      edge_iter_impl e_beg = edge_iter_impl(m_vertices.begin());
      edge_iter_impl e_end = edge_iter_impl(m_vertices.end());
      return std::pair< edge_iterator, edge_iterator >( edge_iterator(e_beg, e_end), edge_iterator(e_end, e_end) );
    };
    
    std::pair<edge_descriptor,bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
      edge_container& u_oe = get_stored_vertex(u).out_edges;
      for(e_iter ei = u_oe.begin(); ei != u_oe.end(); ++ei) {
        if(ei->target == v)
          return std::make_pair(edge_descriptor(u, iterator_to_desc<edge_rand_access>(u_oe, ei)),true);
      };
      return std::make_pair(edge_descriptor(),false);
    };
    
    std::pair< child_vertex_iterator, child_vertex_iterator > child_vertices(vertex_descriptor v) const {
      edge_container& v_oe = get_stored_vertex(v).out_edges;
      return std::pair< child_vertex_iterator, child_vertex_iterator >(
        child_vertex_iterator(v_oe.begin()),
        child_vertex_iterator(v_oe.end())
      );
    };
    
  };
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_linked_list_container {
    
    
    typedef typename ::boost::detail::is_random_access<VertexListS>::type vertex_rand_access;
    typedef typename ::boost::detail::is_random_access<OutEdgeListS>::type edge_rand_access;
    
    typedef ltree_vertex_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> VConfig;
    typedef typename VConfig::container vertex_container;
    typedef typename vertex_container::size_type vertices_size_type;
    typedef typename VConfig::descriptor vertex_descriptor;
    typedef typename VConfig::value_type vertex_value_type;
    typedef vertex_value_type vertex_stored_type;
    
    typedef ltree_edge_config<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> EConfig;
    typedef typename EConfig::container edge_container;
    typedef typename edge_container::size_type edges_size_type;
    typedef typename EConfig::descriptor edge_descriptor;
    typedef typename EConfig::value_type edge_value_type;
    
    typedef typename vertex_container::iterator v_iter;
    typedef typename edge_container::iterator e_iter;
    
    
    typedef ltree_viter_from_iter< vertex_descriptor > vertex_iterator;
    typedef ltree_child_viter_from_iter<vertex_descriptor, e_iter> child_vertex_iterator;
    typedef child_vertex_iterator adjacency_iterator;
    
    typedef typename mpl::if_< edge_rand_access,
      ltree_oeiter_from_index< edge_descriptor >,
      ltree_oeiter_from_iter< edge_descriptor > >::type out_edge_iterator;
    
    typedef const edge_descriptor* in_edge_iterator;
    
    typedef ltree_eiter_from_iter< edge_descriptor, typename vertex_container::const_iterator > edge_iter_impl;
    typedef filter_iterator< ltree_edge_has_a_source<vertex_descriptor>, edge_iter_impl > edge_iterator;
    
    
    mutable vertex_container m_vertices;
    vertex_descriptor m_root;
    
    ltree_linked_list_container() : m_vertices(), m_root(VConfig::null_vertex()) { };
    
    std::size_t size() const { return m_vertices.size(); };
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    
    vertex_stored_type& get_stored_vertex(vertex_descriptor v) const { return *v; };
    
    edge_value_type& get_stored_edge(edge_descriptor e) { 
      return *desc_to_iterator(e.source->out_edges, e.edge_id);
    };
    const edge_value_type& get_stored_edge(edge_descriptor e) const { 
      return *desc_to_iterator(e.source->out_edges, e.edge_id);
    };
    
    
    vertex_descriptor add_new_vertex(const VertexProperties& vp) {
      return m_vertices.insert(m_vertices.end(), vertex_value_type(VertexProperties(vp)));
    };
    void add_root_vertex(const VertexProperties& vp) { m_root = add_new_vertex(vp); };
    
    
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, const VertexProperties& vp, const EdgeProperties& ep) {
      vertex_descriptor new_node = add_new_vertex(vp);
      v_iter v_it = desc_to_iterator(m_vertices, v);
      v_it->out_edges.insert(v_it->out_edges.end(), edge_value_type(new_node, std::move(ep)));
      // update the edge-descriptors in the children nodes.
      for(e_iter ei = v_it->out_edges.begin(); ei != v_it->out_edges.end(); ++ei) {
        edge_descriptor& tmp_e = desc_to_iterator(m_vertices, ei->target)->in_edge;
        tmp_e.source = v;
        tmp_e.edge_id = iterator_to_desc<edge_rand_access>(v_it->out_edges, ei);
      };
      return std::make_pair(new_node, desc_to_iterator(m_vertices, new_node)->in_edge);
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, VertexProperties&& vp, EdgeProperties&& ep) {
      vertex_descriptor new_node = add_new_vertex(std::move(vp));
      v_iter v_it = desc_to_iterator(m_vertices, v);
      v_it->out_edges.insert(v_it->out_edges.end(), edge_value_type(new_node, std::move(ep)));
      // update the edge-descriptors in the children nodes.
      for(e_iter ei = v_it->out_edges.begin(); ei != v_it->out_edges.end(); ++ei) {
        edge_descriptor& tmp_e = desc_to_iterator(m_vertices, ei->target)->in_edge;
        tmp_e.source = v;
        tmp_e.edge_id = iterator_to_desc<edge_rand_access>(v_it->out_edges, ei);
      };
      return std::make_pair(new_node, desc_to_iterator(m_vertices, new_node)->in_edge);
    };
#endif
    
    void clear() { 
      m_vertices.clear();
      m_root = VConfig::null_vertex();
    };
    
    template <typename OutputIter>
    OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
      vertex_descriptor parent = v->in_edge.source;
      std::stack< vertex_descriptor > death_row;
      death_row.push(v);
      std::queue< vertex_descriptor > bft_queue;
      bft_queue.push(v);
      // Put all children on death-row:
      while( ! bft_queue.empty() ) {
        vertex_descriptor tmp_it = bft_queue.front(); 
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(it_out++) = std::move(tmp_it->data);
#else
        *(it_out++) = tmp_it->data;
#endif
        bft_queue.pop(); 
        for(e_iter ei2 = tmp_it->out_edges.begin(); ei2 != tmp_it->out_edges.end(); ++ei2) {
          death_row.push(ei2->target);
          bft_queue.push(ei2->target);
        };
      };
      
      // Check if we removed the (or a) root:
      e_iter ei;
      if(parent == VConfig::null_vertex()) {
        if(v == m_root) {
          clear();
          return it_out;
        };
      } else {
        // find the edge to be removed:
        ei = parent->out_edges.begin();
        while( (ei != parent->out_edges.end()) && (ei->target != v) )
          ++ei;
      };
      
      // Execute them!
      while( ! death_row.empty() ) {
        vertex_descriptor tmp_it = death_row.top(); 
        death_row.pop(); 
        m_vertices.erase(tmp_it);
      };
      
      if(parent == VConfig::null_vertex())
        return it_out;
      
      // remove the edge.
      parent->out_edges.erase(ei);
      // update the edge-descriptors in the children nodes.
      for(ei = parent->out_edges.begin(); ei != parent->out_edges.end(); ++ei)
        ei->target->in_edge.edge_id = iterator_to_desc<edge_rand_access>(parent->out_edges, ei);
      
      return it_out;
    };
    void remove_branch_impl(vertex_descriptor v) { remove_branch_impl(v, ignore_output_iter()); };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    vertex_descriptor add_new_vertex(VertexProperties&& vp) {
      return m_vertices.insert(m_vertices.end(), vertex_value_type(VertexProperties(std::move(vp))));
    };
    void add_root_vertex(VertexProperties&& vp) { m_root = add_new_vertex(std::move(vp)); };
#endif
    
    std::size_t get_depth(const vertex_descriptor& u) const {
      v_iter u_it = desc_to_iterator(m_vertices, u);
      std::size_t max_depth = 0;
      typedef typename edge_container::const_iterator EdgeIter;
      for(EdgeIter ei = u_it->out_edges.begin(); ei != u_it->out_edges.end(); ++ei) {
        std::size_t temp = get_depth( ei->target );
        if(temp > max_depth)
          max_depth = temp;
      };
      return max_depth + 1;
    };
    
    
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor v) const {
      v_iter v_it = desc_to_iterator(m_vertices, v);
      return std::make_pair(
        out_edge_iterator(
          edge_descriptor(
            v, 
            iterator_to_desc<edge_rand_access>(
              v_it->out_edges, 
              v_it->out_edges.begin()
            )
          )
        ),
        out_edge_iterator(
          edge_descriptor(
            v, 
            iterator_to_desc<edge_rand_access>(
              v_it->out_edges, 
              v_it->out_edges.end()
            )
          )
        )
      );
    };
    
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) const {
      const edge_descriptor& tmp_e = get_stored_vertex(v).in_edge;
      if(tmp_e.source == VConfig::null_vertex())
        return std::make_pair(&tmp_e, &tmp_e);
      else
        return std::make_pair(&tmp_e, &tmp_e + 1);
    };
    
    std::pair< vertex_iterator, vertex_iterator > vertices() const {
      return std::pair< vertex_iterator, vertex_iterator >(
        vertex_iterator::begin(m_vertices), 
        vertex_iterator::end(m_vertices)
      );
    };
    
    std::pair< edge_iterator, edge_iterator > edges() const {
      edge_iter_impl e_beg = edge_iter_impl(m_vertices.begin());
      edge_iter_impl e_end = edge_iter_impl(m_vertices.end());
      return std::pair< edge_iterator, edge_iterator >( edge_iterator(e_beg, e_end), edge_iterator(e_end, e_end) );
    };
    
    std::pair<edge_descriptor,bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
      v_iter u_it = desc_to_iterator(m_vertices, u);
      for(e_iter ei = u_it->out_edges.begin(); ei != u_it->out_edges.end(); ++ei) {
        if(ei->target == v)
          return std::make_pair(edge_descriptor(u, iterator_to_desc<edge_rand_access>(u_it->out_edges, ei)),true);
      };
      return std::make_pair(edge_descriptor(),false);
    };
    
    std::pair< child_vertex_iterator, child_vertex_iterator > child_vertices(vertex_descriptor v) const {
      v_iter v_it = desc_to_iterator(m_vertices, v);
      return std::pair< child_vertex_iterator, child_vertex_iterator >(
        child_vertex_iterator(v_it->out_edges.begin()),
        child_vertex_iterator(v_it->out_edges.end())
      );
    };
    
    
  };
  
  
  
  template <typename VertexListS, typename OutEdgeListS, 
            typename VertexProperties, typename EdgeProperties>
  struct ltree_container_select {
    typedef typename ::boost::detail::is_random_access<VertexListS>::type vertex_rand_access;
    typedef typename mpl::if_< vertex_rand_access,
      ltree_pooled_vector_container<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties>,
      ltree_linked_list_container<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties> >::type type;
  };
    
  
  
  
}; 


};


};


#endif


















