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


namespace detail {
  
  template <typename Graph>
  class child_vertex_iter_impl {
    public:
      typedef Graph graph_type;
      typedef typename graph_traits< graph_type >::out_edge_iterator out_edge_iter;
    private:
      out_edge_iter ei;
      graph_type const * g;
    public:
      
      child_vertex_iter_impl() : ei(), g(NULL) { };
      child_vertex_iter_impl(const out_edge_iter& aEi, graph_type const * aG) : ei(aEi), g(aG) { };
      
      typedef child_vertex_iter_impl<Graph> self;
      
      typedef std::ptrdiff_t difference_type;
      typedef typename graph_traits< graph_type >::vertex_descriptor value_type;
      typedef value_type* pointer;
      typedef value_type& reference;
      typedef typename std::iterator_traits<out_edge_iter>::iterator_category iterator_category;
      
      bool operator==(const self& rhs) { return ei == rhs.ei; };
      bool operator!=(const self& rhs) { return ei != rhs.ei; };
      bool operator >(const self& rhs) { return ei > rhs.ei; };
      bool operator >=(const self& rhs) { return ei >= rhs.ei; };
      bool operator <(const self& rhs) { return ei < rhs.ei; };
      bool operator <=(const self& rhs) { return ei <= rhs.ei; };
      
      self& operator++() { ++ei; return *this; };
      self operator++(int) { self result(*this); ++ei; return result; };
      self& operator--() { --ei; return *this; };
      self operator--(int) { self result(*this); --ei; return result; };
      
      self operator+(difference_type i) const {
        return self(ei + i, g);
      };
      self operator-(difference_type i) const {
        return self(ei - i, g);
      };
      difference_type operator-(const self& rhs) const {
        return difference_type(ei - rhs.ei);
      };
      
      self& operator +=(difference_type i) { ei += i; return *this; };
      self& operator -=(difference_type i) { ei -= i; return *this; };
      
      value_type operator[](difference_type i) const { return target(*(ei + i), *g); };
      value_type operator*() { return target(*ei, *g); };
      pointer operator->() { return &target(*ei, *g); };
      
  };
  
  template <typename Graph>
  child_vertex_iter_impl<Graph> operator+(std::ptrdiff_t i, const child_vertex_iter_impl<Graph>& rhs) {
    return rhs + i;
  };
  
  
};


/***********************************************************************************************
 *                             tree traits for adjacency_list
 * ********************************************************************************************/


// forward-declaration:
template <class OutEdgeListS,
          class VertexListS,
          class DirectedS,
          class VertexProperty,
          class EdgeProperty,
          class GraphProperty,
          class EdgeListS>
class adjacency_list;

struct adj_list_tree_storage { };

template <typename VertexProperty = no_property, 
          typename EdgeProperty = no_property, 
          typename GraphProperty = no_property,
          typename TreeStorage = adj_list_tree_storage>
struct tree_storage {
  typedef adjacency_list< vecS, listS, bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS > type;
};

template <typename TreeStorage = adj_list_tree_storage>
struct tree_storage_traits :
  adjacency_list_traits< vecS, listS, bidirectionalS, listS > { };

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty>
struct tree_traits< adjacency_list< vecS, listS, bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS> > {
  typedef detail::child_vertex_iter_impl< adjacency_list< vecS, listS, bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS> > child_vertex_iterator;
};


/***********************************************************************************************
 *                             tree traits for pooled_adjacency_list
 * ********************************************************************************************/

// forward-declaration:
template <typename DirectedS,
          typename VertexProperty,
          typename EdgeProperty,
          typename GraphProperty,
          typename EdgeListS>
class pooled_adjacency_list;

struct pooled_adj_list_tree_storage { };

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty>
struct tree_storage< VertexProperty, EdgeProperty, GraphProperty, pooled_adj_list_tree_storage > {
  typedef pooled_adjacency_list< bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS> type;
};

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty>
struct tree_traits< pooled_adjacency_list< bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS> > {
  typedef detail::child_vertex_iter_impl< pooled_adjacency_list< bidirectionalS, VertexProperty, EdgeProperty, GraphProperty, listS> > child_vertex_iterator;
};




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
  typedef typename tree_traits<Graph>::child_vertex_iterator CVIter;
  typename graph_traits<Graph>::out_edge_iterator ei,ei_end;
  tie(ei,ei_end) = out_edges(v, g);
  return std::make_pair(CVIter(ei,&g),CVIter(ei_end,&g));
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

template <typename Graph >
void remove_branch(const typename graph_traits<Graph>::vertex_descriptor& u, Graph& g) {
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
          typename OutputIter>
OutputIter remove_branch(const typename graph_traits<Graph>::vertex_descriptor& u, OutputIter it_out, Graph& g) {
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




/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/


template <typename Graph>
bool is_vertex_valid( typename graph_traits<Graph>::vertex_descriptor, const Graph&) {
  return true;
};

template <typename Graph>
bool is_edge_valid( typename graph_traits<Graph>::edge_descriptor, const Graph&) {
  return true;
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


















