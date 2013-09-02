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

#include <queue>

namespace boost {


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor
  get_root_vertex( const Graph& g) {
  return *(vertices(g).first);
};

template <typename Graph>
std::pair<typename tree_traits<Graph>::child_vertex_iterator, 
          typename tree_traits<Graph>::child_vertex_iterator>
  child_vertices(typename graph_traits<Graph>::vertex_descriptor v, const Graph& g) {
  return adjacent_vertices(v, g);
};


/***********************************************************************************************
 *                             BidirectionalTreeConcept
 * ********************************************************************************************/

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
 *                             MutableTreeConcept
 * ********************************************************************************************/

template <typename Graph >
typename graph_traits<Graph>::vertex_descriptor
  create_root(Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(g);
};

template <typename Graph >
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex(const typename graph_traits<Graph>::vertex_descriptor& u, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v, e);
};



template <typename Graph, typename DummyTag >
void remove_branch_impl(const typename graph_traits<Graph>::vertex_descriptor& u, 
                        Graph& g, const DummyTag&) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  while( ! v_queue.empty() ) {
    Vertex v = v_queue.front(); v_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei)
      v_queue.push(target(*ei, g));
    clear_vertex(v, g);
    remove_vertex(v, g);
  };
};


template <typename Graph >
void remove_branch_impl(const typename graph_traits<Graph>::vertex_descriptor& u, 
                        Graph& g, const std::random_access_iterator_tag&) {
  /* TODO: Figure out a way to implement this! */
};

template <typename Graph >
void remove_branch(const typename graph_traits<Graph>::vertex_descriptor& u, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  remove_branch_impl(u,g,category());
};




/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor create_root(const typename Graph::vertex_bundled& vp, Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(vp,g);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename Graph>
typename graph_traits<Graph>::vertex_descriptor
  create_root( typename Graph::vertex_bundled&& vp, 
               Graph& g) {
  if(num_vertices(g) > 0)
    remove_branch(get_root_vertex(g), g);
  return add_vertex(std::move(vp),g);
};
#endif

template <typename Graph>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    const typename Graph::vertex_bundled& vp, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(vp, g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v,e);
};

template <typename Graph>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    const typename Graph::vertex_bundled& vp,
                    const typename Graph::edge_bundled& ep, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(vp, g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, ep, g).first;
  return std::make_pair(v,e);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template <typename Graph>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    typename Graph::vertex_bundled&& vp, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(std::move(vp), g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, g).first;
  return std::make_pair(v,e);
};

template <typename Graph>
std::pair< 
typename graph_traits<Graph>::vertex_descriptor,
typename graph_traits<Graph>::edge_descriptor >
  add_child_vertex( const typename graph_traits<Graph>::vertex_descriptor& u,
                    typename Graph::vertex_bundled&& vp,
                    typename Graph::edge_bundled&& ep, Graph& g) {
  typename graph_traits<Graph>::vertex_descriptor v = add_vertex(std::move(vp), g);
  typename graph_traits<Graph>::edge_descriptor e = add_edge(u, v, std::move(ep), g).first;
  return std::make_pair(v,e);
};
#endif



template <typename Graph,
          typename OutputIter,
          typename DummyTag>
OutputIter remove_branch_impl(const typename graph_traits<Graph>::vertex_descriptor& u, 
                              OutputIter it_out, Graph& g, const DummyTag&) {
  typedef typename graph_traits<Graph>::out_edge_iterator EdgeIter;
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  while( ! v_queue.empty() ) {
    Vertex v = v_queue.front(); v_queue.pop();
    EdgeIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei)
      v_queue.push(target(*ei, g));
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    *(it_out++) = std::move(g[v]);
#else
    *(it_out++) = g[v];
#endif
    clear_vertex(v, g);
    remove_vertex(v, g);
  };
  return it_out;
};

template <typename Graph,
          typename OutputIter>
OutputIter remove_branch_impl(const typename graph_traits<Graph>::vertex_descriptor& u, 
                              OutputIter it_out, Graph& g, const std::random_access_iterator_tag&) {
  /* TODO: Figure out a way to implement this! */
};



template <typename Graph,
          typename OutputIter>
OutputIter remove_branch(const typename graph_traits<Graph>::vertex_descriptor& u, OutputIter it_out, Graph& g) {
  typedef typename graph_traits<Graph>::vertex_iterator VIter;
  typedef typename std::iterator_traits<VIter>::iterator_category category;
  return remove_branch_impl(u, it_out, g, category());
};







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


















