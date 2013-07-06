// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


/**
 * \file pooled_adjacency_list.hpp
 * 
 * This library provides a class template that implements a pooled-memory adjacency-list, as pooled_adjacency_list. 
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */


#ifndef BOOST_POOLED_ADJACENCY_LIST_HPP
#define BOOST_POOLED_ADJACENCY_LIST_HPP

#include <boost/config.hpp>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/variant.hpp>
#include <boost/iterator/filter_iterator.hpp>

#include <boost/graph/more_property_maps.hpp>

#include <utility>


namespace boost {



namespace graph { namespace detail {
  
  template <typename Graph>
  struct pooled_adj_list_validity_pred {
    const Graph* p_graph;
    explicit pooled_adj_list_validity_pred(const Graph* aPGraph = NULL) : p_graph(aPGraph) { };
    template <typename Desc>
    bool operator()(Desc d) {
      return ((*p_graph)[d].which() == 0);
    };
  };
  
}; };


// NOTE: The class and functions below will not work when BOOST_GRAPH_NO_BUNDLED_PROPERTIES is defined.
//       This is because operator[] of adjacency_list are disabled in that case.
//       However, there should be a way to enable those operators regardless.
//       This is because BOOST_GRAPH_NO_BUNDLED_PROPERTIES is linked to not being able to use 
//       template specializations based on a pointer-to-member type, but this should only lead 
//       to calls like get(&Bundle::Member, g) to be impossible to realize (due to having to 
//       use property_map< Graph, T Bundle::* > specialization).


/**
 * This class template is a simple wrapper of the boost::adjacency_list class. This pooled-memory 
 * adjacency-list is uses a vector storage of the vertex-list without invalidating vertex-descriptors 
 * and vertex-iterators when removing nodes from the graph. In other words, instead of removing nodes 
 * from the graph, it invalidates them and revives them later when nodes are added to the graph, meaning 
 * that as long as the rate of insertion is greater than or equal to the rate of removal of nodes, 
 * there isn't much waste of memory. The use of a vector storage allows for more efficient traversals
 * through the vertices due to locality of reference, but having to internally "skip" invalid vertices 
 * implies a certain overhead (test validity) and means that random-access with vertex-iterators is not possible.
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */
template <typename DirectedS = directedS,
          typename VertexProperty = no_property,
          typename EdgeProperty = no_property,
          typename GraphProperty = no_property,
          typename EdgeListS = listS>
class pooled_adjacency_list {
    
  public:
    
    typedef pooled_adjacency_list self;
    
    struct hole_descriptor {
      std::size_t value;
      explicit hole_descriptor(std::size_t aValue = 0) : value(aValue) { };
    };
    
    // PropertyGraph traits:
    typedef EdgeProperty edge_property_type;
    typedef VertexProperty vertex_property_type;
    
    typedef variant<VertexProperty, hole_descriptor> raw_vertex_property_type;
    typedef edge_property_type raw_edge_property_type;
    
    typedef edge_property_type edge_bundled;
    typedef vertex_property_type vertex_bundled;
    
    typedef adjacency_list< 
      vecS, vecS, DirectedS, 
      raw_vertex_property_type, 
      edge_property_type, 
      GraphProperty, EdgeListS > graph_type;
    typedef graph_traits< graph_type > Traits;
    
    // Graph traits:
    typedef typename Traits::vertex_descriptor vertex_descriptor;
    typedef typename Traits::edge_descriptor edge_descriptor;
    typedef typename Traits::directed_category directed_category;
    typedef typename Traits::edge_parallel_category edge_parallel_category;
    typedef typename Traits::traversal_category traversal_category;
    
    static vertex_descriptor null_vertex() { return Traits::null_vertex(); };
    
    // IncidenceGraph traits:
    typedef typename Traits::out_edge_iterator out_edge_iterator;
    typedef typename Traits::degree_size_type degree_size_type;
    
    // BidirectionalGraph traits:
    typedef typename Traits::in_edge_iterator in_edge_iterator;
    
    // VertexListGraph traits:
    typedef filter_iterator<
      graph::detail::pooled_adj_list_validity_pred<graph_type>, 
      typename Traits::vertex_iterator > vertex_iterator;
    
    typedef typename Traits::vertices_size_type vertices_size_type;
    
    // EdgeListGraph traits:
    typedef typename Traits::edge_iterator edge_iterator;
    typedef typename Traits::edges_size_type edges_size_type;
    
    // AdjacencyGraph traits:
    typedef filter_iterator<
      graph::detail::pooled_adj_list_validity_pred<graph_type>, 
      typename Traits::adjacency_iterator > adjacency_iterator;
    
    typedef typename graph_type::graph_bundled graph_bundled;
    
  public:  // private:  would be private is friends were more portable.
    
    graph_type m_graph;
    
    hole_descriptor m_first_hole;
    vertices_size_type m_num_vertices;
    
  public:
    
    pooled_adjacency_list() : m_graph(), m_first_hole(null_vertex()), m_num_vertices(0) { };
    
    
    // Bundled Property-map functions (used by the property_map< self, T Bundle::* > classes).
    
    vertex_bundled& operator[]( const vertex_descriptor& v_i) {
      return get<vertex_bundled>(m_graph[v_i]);
    };
    const vertex_bundled& operator[]( const vertex_descriptor& v_i) const {
      return get<vertex_bundled>(m_graph[v_i]);
    };
    edge_bundled& operator[]( const edge_descriptor& e_i) {
      return m_graph[e_i];
    };
    const edge_bundled& operator[]( const edge_descriptor& e_i) const {
      return m_graph[e_i];
    };
    
    graph_bundled& operator[](graph_bundle_t gb) { return m_graph[gb]; };
    const graph_bundled& operator[](graph_bundle_t gb) const { return m_graph[gb]; };
    
    
    /*
#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
    // Directly access a vertex or edge bundle
    vertex_bundled& operator[](vertex_descriptor v)
    { return get(vertex_bundle, *this)[v]; }

    const vertex_bundled& operator[](vertex_descriptor v) const
    { return get(vertex_bundle, *this)[v]; }

    edge_bundled& operator[](edge_descriptor e)
    { return get(edge_bundle, *this)[e]; }

    const edge_bundled& operator[](edge_descriptor e) const
    { return get(edge_bundle, *this)[e]; }

    graph_bundled& operator[](graph_bundle_t)
    { return get_property(*this); }

    graph_bundled const& operator[](graph_bundle_t) const
    { return get_property(*this); }
#endif
    */
    
    
  private:
    
    template <typename DirectedCategory>
    void rewire_in_edges_impl(vertex_descriptor, vertex_descriptor, DirectedCategory) { };
    
    void rewire_in_edges_impl(vertex_descriptor u, vertex_descriptor v, directed_tag) {
      in_edge_iterator in_ei, in_ei_end;
      for(tie(in_ei, in_ei_end) = in_edges(v, m_graph); in_ei != in_ei_end; ++in_ei) {
        std::pair<edge_descriptor, bool> e = add_edge(source(*in_ei, m_graph), u, m_graph);
        swap(m_graph[e], m_graph[*in_ei]);
      };
    };
  
  public:
    
    void repack() {
      using std::swap;
      
      if((m_first_hole == null_vertex()) || (num_vertices(m_graph) == 0)) // nothing to re-pack.
        return;
      
      // a simple remove loop won't cut it.
      // the idea here is to achieve a behavior similar to the std::remove function, but by taking
      // elements from the back of the vertex list and swapping them with the first hole.
      // including the re-wiring of edges.
      vertex_descriptor v_end = reinterpret_cast<vertex_descriptor>(num_vertices(m_graph));
      vertex_descriptor v = m_first_hole.value;
      while(v != null_vertex()) {
        do {
          --v_end;
          while((v != null_vertex()) && (v > v_end))
            v = get< hole_descriptor >(m_graph[v]).value;
        } while( ( v_end > 0 ) && ( m_graph[v_end].which() == 1 ) ); 
        // at this point, v is a hole, and v_end is the highest-index valid vertex.
        if((v_end == 0) || (v == null_vertex()))
          break;
        
        swap(m_graph[v], m_graph[v_end]);
        out_edge_iterator ei, ei_end;
        for(tie(ei,ei_end) = out_edges(v_end, m_graph); ei != ei_end; ++ei) {
          std::pair<edge_descriptor, bool> e = add_edge(v, target(*ei, m_graph), m_graph);
          swap(m_graph[e],m_graph[*ei]);
        };
        rewire_in_edges_impl(v, v_end, directed_category());
        clear_vertex(v_end, m_graph);
        
        v = get< hole_descriptor >(m_graph[v_end]).value;
      };
      m_first_hole.value = null_vertex();
      // at this point v_end is the first invalid vertex remaining.
      vertex_descriptor v_rm = reinterpret_cast<vertex_descriptor>(num_vertices(m_graph));
      while(v_rm > v_end)
        remove_vertex(--v_rm, m_graph);
    };

};




/******************************************************************************************
 * *************************************************************************************
 *                             Graph functions
 * *************************************************************************************
 * ***************************************************************************************/

/*******************************************************************************************
 *                  IncidenceGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::out_edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::out_edge_iterator >
  out_edges(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
            const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return out_edges(v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  source(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return source(e, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  target(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return target(e, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  out_degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
             const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return out_degree(v, g.m_graph);
};


/*******************************************************************************************
 *                  BidirectionalGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::in_edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::in_edge_iterator >
  in_edges(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return in_edges(v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  in_degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
            const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return in_degree(v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::degree_size_type
  degree(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
         const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return degree(e, g.m_graph);
};


/*******************************************************************************************
 *                  AdjacencyGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::adjacency_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::adjacency_iterator >
  adjacent_vertices(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                    const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> Graph;
  typedef graph::detail::pooled_adj_list_validity_pred<typename Graph::graph_type> Predicate;
  typedef typename graph_traits< Graph >::adjacency_iterator AdjIter;
  typedef typename Graph::graph_type BaseGraph;
  typedef typename graph_traits< BaseGraph >::adjacency_iterator BaseAdjIter;
  
  std::pair< BaseAdjIter, BaseAdjIter > tmp_res = adjacent_vertices(u, g.m_graph);
  return std::pair< AdjIter, AdjIter >(
    AdjIter(Predicate(&g), tmp_res.first,  tmp_res.second),
    AdjIter(Predicate(&g), tmp_res.second, tmp_res.second)
  );
};


/*******************************************************************************************
 *                  VertexListGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_iterator >
  vertices(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> Graph;
  typedef graph::detail::pooled_adj_list_validity_pred<typename Graph::graph_type> Predicate;
  typedef typename graph_traits< Graph >::vertex_iterator VIter;
  typedef typename Graph::graph_type BaseGraph;
  typedef typename graph_traits< BaseGraph >::vertex_iterator BaseVIter;
  
  std::pair< BaseVIter, BaseVIter > tmp_res = vertices(g.m_graph);
  return std::pair< VIter, VIter >(
    VIter(Predicate(&g.m_graph), tmp_res.first,  tmp_res.second),
    VIter(Predicate(&g.m_graph), tmp_res.second, tmp_res.second)
  );
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertices_size_type
  num_vertices(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return g.m_num_vertices;
};


/*******************************************************************************************
 *                  EdgeListGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair< typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator >
  edges(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return edges(g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edges_size_type
  num_edges(const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return num_edges(g.m_graph);
};


/*******************************************************************************************
 *                  AdjacencyMatrix concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
       typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
       const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return edge(u, v, g.m_graph);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> Graph;
  typedef typename graph_traits< Graph >::vertex_descriptor Vertex;
  typedef typename Graph::vertex_property_type VProp;
  typedef typename Graph::hole_descriptor Hole;
  
  Vertex v = Vertex();
  if( g.m_first_hole.value == Graph::null_vertex() ) {
    // add a new node in the graph (this is safe, it will not invalidate vertex descriptors).
    v = add_vertex(g.m_graph);
  } else {
    // resurrect a node from the graveyard.
    v = g.m_first_hole.value;
    g.m_first_hole = get< Hole >(g.m_graph[v]);
  };
  
  g.m_graph[v] = VProp();
  ++(g.m_num_vertices);
  
  return v;
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void clear_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                  pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  clear_vertex(v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                   pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  clear_vertex(v, g.m_graph);
  g.m_graph[v] = g.m_first_hole;
  g.m_first_hole.value = v;
  --(g.m_num_vertices);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return add_edge(u, v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                 typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(u, v, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(e, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator e_iter,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(e_iter, g.m_graph);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type& vp,
             pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor Vertex;
  Vertex v = add_vertex(g);
  g.m_graph[v] = vp;
  return v;
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_vertex(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                   typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type& vp,
                   pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  vp = g.m_graph[v];
  remove_vertex(v, g);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return add_edge(u, v, ep, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                 typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(u, v, ep, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(e, ep, g.m_graph);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void remove_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_iterator e_iter,
                 typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type& ep,
                 pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  remove_edge(e_iter, ep, g.m_graph);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor
  add_vertex(typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_property_type&& vp,
             pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor Vertex;
  Vertex v = add_vertex(g);
  g.m_graph[v] = std::move(vp);
  return v;
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
std::pair<typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor, bool>
  add_edge(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
           typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor v,
           typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_property_type&& ep,
           pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return add_edge(u, v, std::move(ep), g.m_graph);
};

#endif

/***********************************************************************************************
 *                             NonCompactGraphConcept
 * ********************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
bool is_vertex_valid(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::vertex_descriptor u,
                     const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  typedef pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> Graph;
  if( ( u != Graph::null_vertex() ) && ( g.m_graph[u].which() == 0 ) )
    return true;
  else
    return false;
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
bool is_edge_valid(typename graph_traits< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS> >::edge_descriptor e,
                   const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return true;
};



/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/


template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_bundled& 
get( const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, 
     typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_descriptor v_i) 
{
  return g[v_i];
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_bundled& 
get( const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, 
     typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_descriptor e_i) 
{
  return g[e_i];
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void put( pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, 
          typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_descriptor v_i, 
          const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_bundled& value) {
  g[v_i] = value;
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS>
void put( pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, 
          typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_descriptor e_i, 
          const typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::edge_bundled& value) {
  g[e_i] = value;
};



#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, typename T, typename Bundle>
struct property_map< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_same< non_const_Bundle, typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::vertex_bundled > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS,
          typename T, typename Bundle>
typename property_map< 
  pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, 
  T Bundle::* 
>::type
get( T Bundle::* p, pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return typename property_map< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, T Bundle::* >::type(&g, p);
};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS,
          typename T, typename Bundle>
typename property_map< 
  pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, 
  T Bundle::* 
>::const_type
get( T Bundle::* p, const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g) {
  return typename property_map< pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, T Bundle::* >::const_type(&g, p);
};

#endif


/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, 
          typename Property>
struct property_map<pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>, Property>
  : property_map<
      typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::graph_type, 
      Property> {};

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, 
          typename Property>
typename property_map< typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::graph_type, Property>::type
get(Property p, pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g)
{
  typedef typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::graph_type BaseGraph;
  return get(p, g.m_graph);
}

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, 
          typename Property>
typename property_map< typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::graph_type, Property>::const_type
get(Property p, const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g)
{
  return get(p, g.m_graph);
}

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, 
          typename Property, typename Key>
typename property_map_value< typename pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>::graph_type, Property>::type
get(Property p, const pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, const Key& k)
{
  return get(p, g.m_graph, k);
}

template <typename DirectedS, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename EdgeListS, 
          typename Property, typename Key, typename Value>
void
put(Property p, pooled_adjacency_list<DirectedS, VertexProperty, EdgeProperty, GraphProperty, EdgeListS>& g, 
    const Key& k, const Value& val)
{
  put(p, g.m_graph, k, val);
}




};


#endif


















