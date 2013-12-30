// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file tree_adaptor.hpp
 * 
 * This library provides function templates to adapt a MutableGraph from the Boost Graph Library 
 * such that it has the mutable interface of a tree structure.
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */

#ifndef BOOST_TREE_ADAPTOR_HPP
#define BOOST_TREE_ADAPTOR_HPP

#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/tree_traits.hpp>

#include <boost/graph/detail/boost_container_generators.hpp>

#include <boost/type_traits/is_convertible.hpp>

#include <queue>

namespace boost {


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor
  get_root_vertex( const Graph& g) {
  return *(vertices(g).first);
};

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template <typename Graph>
std::pair<typename tree_traits<Graph>::child_vertex_iterator, 
          typename tree_traits<Graph>::child_vertex_iterator>
  child_vertices(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g) {
  return adjacent_vertices(v, g);
};


/***********************************************************************************************
 *                             BidirectionalTreeConcept
 * ********************************************************************************************/

/**
 * Returns the parent vertex of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The parent vertex of the given vertex descriptor (will be null_vertex() if it is the root (no parent)).
 */
template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor
  parent_vertex(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g) {
  typedef typename graph_traits<Graph>::in_edge_iterator InEdgeIter;
  InEdgeIter ei, ei_end;
  tie(ei, ei_end) = in_edges(v, g);
  if(ei == ei_end)
    return graph_traits<Graph>::null_vertex();
  else
    return source(*ei, g);
};




/***********************************************************************************************
 *                             MutablePropertyTreeConcept / MutableTreeConcept
 * ********************************************************************************************/

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor create_root(Graph& g);

#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <typename Graph, typename VProp>
typename graph_traits<Graph>::vertex_descriptor create_root(const VProp& vp, Graph& g);
#else
/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <typename Graph, typename VProp>
typename graph_traits<Graph>::vertex_descriptor create_root( VProp&& vp, Graph& g);
#endif



/**
 * Adds a child vertex to the given parent vertex, and default-initializes the properties of 
 * the newly created vertex and edge.
 * \param v The parent vertex to which a child will be added.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <typename Graph >
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( typename graph_traits<Graph>::vertex_descriptor u, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v, e);
};

#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
 * vertex to the given property value.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <typename Graph, typename VProp>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( typename graph_traits<Graph>::vertex_descriptor u,
                    const VProp& vp, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(vp, g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v,e);
};

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
 * vertex and edge to the given property values.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param ep The property value to be moved into the newly created edge.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <typename Graph, typename VProp, typename EProp>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( typename graph_traits<Graph>::vertex_descriptor u,
                    const VProp& vp, const EProp& ep, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(vp, g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, ep, g).first;
  return std::make_pair(v,e);
};

#else

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
 * vertex to the given property value.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <typename Graph, typename VProp>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( typename graph_traits<Graph>::vertex_descriptor u,
                    VProp&& vp, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(std::forward<VProp>(vp), g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v,e);
};

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
 * vertex and edge to the given property values.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param ep The property value to be moved into the newly created edge.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <typename Graph, typename VProp, typename EProp>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    VProp&& vp, EProp&& ep, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(std::forward<VProp>(vp), g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, std::forward<EProp>(ep), g).first;
  return std::make_pair(v,e);
};

#endif



namespace graph { namespace detail {

  template <typename Graph>
  void record_vertex_prop_impl(typename graph_traits<Graph>::vertex_descriptor, ignore_output_iter&, Graph&) { 
    /* do nothing */
  };

  template <typename Graph, typename VertexOIter>
  void record_vertex_prop_impl(typename graph_traits<Graph>::vertex_descriptor u, VertexOIter& vit_out, Graph& g) {
  #ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    *(vit_out++) = std::move(g[u]);
  #else
    *(vit_out++) = g[u];
  #endif
  };

  template <typename Graph>
  void record_edge_prop_impl(typename graph_traits<Graph>::edge_descriptor, ignore_output_iter&, Graph&) { 
    /* do nothing */
  };

  template <typename Graph, typename EdgeOIter>
  void record_edge_prop_impl(typename graph_traits<Graph>::edge_descriptor e, EdgeOIter& eit_out, Graph& g) {
  #ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    *(eit_out++) = std::move(g[e]);
  #else
    *(eit_out++) = g[e];
  #endif
  };
  
  
  template <typename Graph>
  void record_in_edge_prop_impl(typename graph_traits<Graph>::vertex_descriptor, ignore_output_iter&, Graph&) {
    /* do nothing */
  };
  
  template <typename Graph, typename EdgeOIter>
  typename enable_if< is_convertible< typename graph_traits<Graph>::traversal_category*, bidirectional_graph_tag* >,
  void >::type record_in_edge_prop_impl(typename graph_traits<Graph>::vertex_descriptor u, EdgeOIter& eit_out, Graph& g) {
    /* depends on bidirectionality of the graph type: */
    typedef typename graph_traits<Graph>::in_edge_iterator InEdgeIter;
    typedef typename Graph::edge_bundled EdgeProp;
    InEdgeIter iei, iei_end;
    tie(iei, iei_end) = in_edges(u, g);
    if(iei == iei_end) {
      *(eit_out++) = EdgeProp();
    } else {
      record_edge_prop_impl(*iei, eit_out, g);
    };
  };

  template <typename Graph, typename EdgeOIter>
  typename disable_if< is_convertible< typename graph_traits<Graph>::traversal_category*, bidirectional_graph_tag* >,
  void >::type record_in_edge_prop_impl(typename graph_traits<Graph>::vertex_descriptor u, EdgeOIter& eit_out, Graph& g) {
    /* fall-back solution for a uni-directional graph type: */
    typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
    typedef typename Graph::edge_bundled EdgeProp;
    typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
    
    std::queue<Vertex> v_queue;
    Vertex v = get_root_vertex(g);
    if( v == u ) {
      *(eit_out++) = EdgeProp();  // fill in a dummy edge-property.
      return;
    };
    v_queue.push(v);
    
    while( ! v_queue.empty() ) {
      v = v_queue.front(); v_queue.pop();
      EdgeIter ei, ei_end;
      for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
        if( target(*ei, g) == u ) {
          record_edge_prop_impl(*ei, eit_out, g);
          return;
        };
        v_queue.push(target(*ei, g));
      };
    };
  };
  
  
}; };



template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter,
          typename DummyTag>
std::pair<VertexOIter, EdgeOIter> clear_children_impl(
    typename graph_traits<Graph>::vertex_descriptor u, 
    VertexOIter vit_out, EdgeOIter eit_out, Graph& g, const DummyTag&) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  while( ! v_queue.empty() ) {
    Vertex v = v_queue.front(); v_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      v_queue.push(target(*ei, g));
      graph::detail::record_vertex_prop_impl(target(*ei, g), vit_out, g);
      graph::detail::record_edge_prop_impl(*ei, eit_out, g);
    };
    if( v != u ) {
      clear_vertex(v, g);
      remove_vertex(v, g);
    };
  };
  return std::pair<VertexOIter, EdgeOIter>(vit_out,eit_out);
};

template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children_impl(
    typename graph_traits<Graph>::vertex_descriptor u, 
    VertexOIter vit_out, EdgeOIter eit_out, Graph& g, const std::random_access_iterator_tag&) {
  /* TODO: Figure out a way to implement this better! */
  
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  using std::swap;
  
  Graph g_tmp;
  std::queue<Vertex> v_tmp_queue;
  std::queue<Vertex> v_ori_queue;
  Vertex v = get_root_vertex(g);
  v_ori_queue.push(v);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  v_tmp_queue.push(create_root(std::move(g[v]), g_tmp));
#else
  v_tmp_queue.push(create_root(g[v], g_tmp));
#endif
  
  while( ! v_ori_queue.empty() ) {
    v = v_ori_queue.front(); v_ori_queue.pop();
    Vertex v_tmp = v_tmp_queue.front(); v_tmp_queue.pop();
    if( v == u )
      continue;
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      v_ori_queue.push(target(*ei, g));
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      v_tmp_queue.push(add_child_vertex(v_tmp, std::move(g[target(*ei, g)]), std::move(g[*ei]), g_tmp).first);
#else
      v_tmp_queue.push(add_child_vertex(v_tmp, g[target(*ei, g)], g[*ei], g_tmp).first);
#endif
    };
  };
  
  v_ori_queue.push(u);
  while( ! v_ori_queue.empty() ) {
    v = v_ori_queue.front(); v_ori_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      v_ori_queue.push(target(*ei, g));
      graph::detail::record_vertex_prop_impl(target(*ei, g), vit_out, g);
      graph::detail::record_edge_prop_impl(*ei, eit_out, g);
    };
  };
  
  swap(g, g_tmp);
  
  return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
};



/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template <typename Graph>
void clear_children(typename graph_traits<Graph>::vertex_descriptor u, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  clear_children_impl(u, graph::detail::ignore_output_iter(), graph::detail::ignore_output_iter(), g, category());
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <typename Graph,
          typename OutputIter>
OutputIter clear_children(typename graph_traits<Graph>::vertex_descriptor u, OutputIter it_out, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  return clear_children_impl(u, it_out, graph::detail::ignore_output_iter(), g, category()).first;
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children(typename graph_traits<Graph>::vertex_descriptor u, VertexOIter vit_out, EdgeOIter eit_out, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  return clear_children_impl(u, vit_out, eit_out, g, category());
};






template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter,
          typename DummyTag>
std::pair<VertexOIter, EdgeOIter> remove_branch_impl(
    typename graph_traits<Graph>::vertex_descriptor u, 
    VertexOIter vit_out, EdgeOIter eit_out, Graph& g, const DummyTag&) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  graph::detail::record_vertex_prop_impl(u, vit_out, g);
  graph::detail::record_in_edge_prop_impl(u, eit_out, g);
  
  while( ! v_queue.empty() ) {
    Vertex v = v_queue.front(); v_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      v_queue.push(target(*ei, g));
      graph::detail::record_vertex_prop_impl(target(*ei, g), vit_out, g);
      graph::detail::record_edge_prop_impl(*ei, eit_out, g);
    };
    clear_vertex(v, g);
    remove_vertex(v, g);
  };
  
  return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
};

template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch_impl(
    typename graph_traits<Graph>::vertex_descriptor u, 
    VertexOIter vit_out, EdgeOIter eit_out, Graph& g, const std::random_access_iterator_tag&) {
  /* TODO: Figure out a way to implement this better! */
  
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename Graph::edge_bundled EdgeProp;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  using std::swap;
  
  Graph g_tmp;
  std::queue<Vertex> v_tmp_queue;
  std::queue<Vertex> v_ori_queue;
  Vertex v = get_root_vertex(g);
  if( v != u ) {
    v_ori_queue.push(v);
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    v_tmp_queue.push(create_root(std::move(g[v]), g_tmp));
#else
    v_tmp_queue.push(create_root(g[v], g_tmp));
#endif
  } else {
    graph::detail::record_vertex_prop_impl(v, vit_out, g);
    *(eit_out++) = EdgeProp();  // fill in a dummy edge-property.
  };
  
  while( ! v_ori_queue.empty() ) {
    v = v_ori_queue.front(); v_ori_queue.pop();
    Vertex v_tmp = v_tmp_queue.front(); v_tmp_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      if( target(*ei, g) == u ) {
        graph::detail::record_vertex_prop_impl(u, vit_out, g);
        graph::detail::record_edge_prop_impl(*ei, eit_out, g);
        continue;
      };
      v_ori_queue.push(target(*ei, g));
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      v_tmp_queue.push(add_child_vertex(v_tmp, std::move(g[target(*ei, g)]), std::move(g[*ei]), g_tmp).first);
#else
      v_tmp_queue.push(add_child_vertex(v_tmp, g[target(*ei, g)], g[*ei], g_tmp).first);
#endif
    };
  };
  
  v_ori_queue.push(u);
  while( ! v_ori_queue.empty() ) {
    v = v_ori_queue.front(); v_ori_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
      v_ori_queue.push(target(*ei, g));
      graph::detail::record_vertex_prop_impl(target(*ei, g), vit_out, g);
      graph::detail::record_edge_prop_impl(*ei, eit_out, g);
    };
  };
  
  swap(g, g_tmp);
  
  return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
};



/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template <typename Graph>
void remove_branch(typename graph_traits<Graph>::vertex_descriptor u, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  remove_branch_impl(u, graph::detail::ignore_output_iter(), graph::detail::ignore_output_iter(), g, category());
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template <typename Graph,
          typename OutputIter>
OutputIter remove_branch(typename graph_traits<Graph>::vertex_descriptor u, OutputIter it_out, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  return remove_branch_impl(u, it_out, graph::detail::ignore_output_iter(), g, category()).first;
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while 
 * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template <typename Graph,
          typename VertexOIter,
          typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(typename graph_traits<Graph>::vertex_descriptor u, VertexOIter vit_out, EdgeOIter eit_out, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  return remove_branch_impl(u, vit_out, eit_out, g, category());
};



template <typename Graph >
typename graph_traits<Graph>::vertex_descriptor
  create_root(Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(g);
};

#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename Graph, typename VProp>
typename graph_traits<Graph>::vertex_descriptor create_root(const VProp& vp, Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(vp,g);
};
#else
template <typename Graph, typename VProp>
typename graph_traits<Graph>::vertex_descriptor
  create_root( VProp&& vp, Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(std::forward<VProp>(vp),g);
};
#endif




/***********************************************************************************************
 *                             PropertyGraphConcept
 * ********************************************************************************************/

template <typename Graph>
typename vertex_bundle_type<Graph>::type& get_property(typename graph_traits<Graph>::vertex_descriptor v_i, Graph& g) {
  return g[v_i];
};

template <typename Graph>
const typename vertex_bundle_type<Graph>::type& get_property(typename graph_traits<Graph>::vertex_descriptor v_i, const Graph& g) {
  return g[v_i];
};

template <typename Graph>
typename edge_bundle_type<Graph>::type& get_property(typename graph_traits<Graph>::edge_descriptor e_i, Graph& g) {
  return g[e_i];
};

template <typename Graph>
const typename edge_bundle_type<Graph>::type& get_property(typename graph_traits<Graph>::edge_descriptor e_i, const Graph& g) {
  return g[e_i];
};



};


#endif


















