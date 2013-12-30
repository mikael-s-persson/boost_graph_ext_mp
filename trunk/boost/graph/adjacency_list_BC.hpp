// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file adjacency_list_BC.hpp
 * 
 * This library provides a class that implements an adjacency-list structure based on Boost.Containers. 
 * This is a classic tree implementation in which each node contain a list of edges to its children, 
 * and a link to its parent.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date July 2013
 */

#ifndef BOOST_ADJACENCY_LIST_BC_HPP
#define BOOST_ADJACENCY_LIST_BC_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/detail/adjlistBC_iterators.hpp>
#include <boost/graph/detail/adjlistBC_containers.hpp>
#include <boost/graph/more_property_maps.hpp>
#include <boost/unordered_map.hpp>

#include <boost/graph/tree_traits.hpp>

#include <utility>


namespace boost {


template <class OutEdgeListS = vecBC,
          class VertexListS = vecBC,
          class DirectedS = directedS>
struct adjacency_list_BC_traits
{
  typedef typename detail::is_random_access_BC<VertexListS>::type is_rand_access;
  typedef typename DirectedS::is_bidir_t is_bidir;
  typedef typename DirectedS::is_directed_t is_directed;
  
  typedef typename mpl::if_<is_bidir,
    bidirectional_tag,
    typename mpl::if_<is_directed,
      directed_tag, undirected_tag
    >::type
  >::type directed_category;
  
  typedef typename detail::parallel_edge_BC_traits<OutEdgeListS>::type edge_parallel_category;
  
  typedef graph::detail::adjlistBC_traversal_tag<DirectedS> traversal_category;
  
  typedef std::size_t vertices_size_type;
  typedef std::size_t edges_size_type;
  
};


struct adj_list_BC_tag { };

template <typename VertexListS>
struct adjacency_list_BC_disallowed_vertex_list {
  typedef void type;
};

template <> struct adjacency_list_BC_disallowed_vertex_list<setBC> { };
template <> struct adjacency_list_BC_disallowed_vertex_list<multisetBC> { };
template <> struct adjacency_list_BC_disallowed_vertex_list<unordered_setBC> { };
template <> struct adjacency_list_BC_disallowed_vertex_list<unordered_multisetBC> { };



// forward-declare:
template <typename OutEdgeListS = vecBC,
          typename VertexListS = vecBC,
          typename DirectedS = directedS,
          typename VertexProperties = no_property,
          typename EdgeProperties = no_property>
class adjacency_list_BC;


template <typename OutEdgeListS, typename VertexListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
void do_graph_deep_copy(adjacency_list_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>& lhs, 
                        const adjacency_list_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>& rhs);


/**
 * This class implements an adjacency-list based on Boost.Containers that is tailored 
 * to store elements of a graph.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam DirectedS A type tag to choose the directionality of the edges.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS,
          typename VertexListS,
          typename DirectedS,
          typename VertexProperties,
          typename EdgeProperties>
class adjacency_list_BC
{
  public:
    typedef typename adjacency_list_BC_disallowed_vertex_list<VertexListS>::type check_allowed_vertex_list;
    
    typedef adjacency_list_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties> self;
    
    typedef OutEdgeListS out_edge_list_selector;
    typedef VertexListS vertex_list_selector;
    typedef DirectedS directed_selector;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    typedef typename graph::detail::adjlistBC_vertex_container<
      VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties> storage_type;
    
    
    typedef typename storage_type::vertex_descriptor vertex_descriptor;
    typedef typename storage_type::vertices_size_type vertices_size_type;
    
    typedef typename storage_type::edge_descriptor edge_descriptor;
    typedef typename storage_type::edges_size_type edges_size_type;
    typedef edges_size_type degree_size_type;
    
    
    typedef typename storage_type::vertex_iterator vertex_iterator;
    
    typedef typename storage_type::out_edge_iterator out_edge_iterator;
    typedef typename storage_type::in_edge_iterator in_edge_iterator;
    typedef typename storage_type::edge_iterator edge_iterator;
    
    typedef typename adjacency_iterator_generator<self, vertex_descriptor, out_edge_iterator>::type adjacency_iterator;
    typedef typename inv_adjacency_iterator_generator<self, vertex_descriptor, in_edge_iterator>::type inv_adjacency_iterator;
    
    typedef typename storage_type::vertex_stored_type vertex_stored_impl;
    typedef typename storage_type::vertex_value_type vertex_value_impl;
    
    typedef typename storage_type::edge_stored_type edge_stored_impl;
    typedef typename storage_type::edge_value_type edge_value_impl;
    
    
    typedef adjacency_list_BC_traits<OutEdgeListS, VertexListS, DirectedS> Traits;
    
    typedef typename Traits::directed_category directed_category;
    typedef typename Traits::edge_parallel_category edge_parallel_category;
    typedef typename Traits::traversal_category traversal_category;
    
    
    typedef adj_list_BC_tag graph_tag;
    
    
    /**
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() { return graph::detail::BC_null_desc<vertex_descriptor>::value(); };
    
    /**
     * This static member function outputs the null-edge (invalid edge descriptor).
     * \return A null-edge descriptor (invalid edge descriptor).
     */
    static edge_descriptor null_edge() { return graph::detail::BC_null_desc<edge_descriptor>::value(); };
    
    
  public:  // private:
    
    storage_type m_pack;
    
  public:
    
    
    /**
     * Constructs an empty adjacency-list.
     */
    adjacency_list_BC() : m_pack() { };
    
    
    adjacency_list_BC(const self& rhs) : m_pack() { do_graph_deep_copy(*this,rhs); };
    self& operator=(const self& rhs) {
      do_graph_deep_copy(*this,rhs);
      return *this;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    adjacency_list_BC(self&& rhs) : m_pack(std::move(rhs.m_pack)) { };
    self& operator=(self&& rhs) {
      m_pack = std::move(rhs.m_pack);
      return *this;
    };
#endif
    
    /**
     * Swaps the adjacency-list with another.
     */
    void swap(self& rhs) { 
      m_pack.swap(rhs.m_pack);
    };
    
    /**
     * Clears the adjacency-list of all vertices and edges.
     */
    void clear() { 
      m_pack.clear();
    };
    
    /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[](vertex_descriptor v) {
      return m_pack.get_stored_vertex(v).data;
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[](vertex_descriptor v) const {
      return m_pack.get_stored_vertex(v).data;
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[](const edge_descriptor& e) {
      return m_pack.get_stored_edge(e).data;
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[](const edge_descriptor& e) const {
      return m_pack.get_stored_edge(e).data;
    };
    
};




// the safest adj-list tree-storage style is <boost::vecBC, boost::listBC, boost::bidirectionalS>

/**
 * This is the tree-storage specifier for an adjacency-list being used as a tree storage.
 */
template <typename OutEdgeListS = vecBC, typename VertexListS = listBC, typename DirectedS = directedS>
struct adj_list_BC_as_tree_storage { };

template <typename VertexProperty, typename EdgeProperty, typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage< VertexProperty, EdgeProperty, adj_list_BC_as_tree_storage<OutEdgeListS, VertexListS, DirectedS> > {
  typedef adjacency_list_BC< OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty> type;
};

template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage_traits< adj_list_BC_as_tree_storage<OutEdgeListS, VertexListS, DirectedS> > :
  adjacency_list_BC_traits<OutEdgeListS, VertexListS, DirectedS> { };

template <typename VertexProperty, typename EdgeProperty, typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_traits< adjacency_list_BC< OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty > > {
  typedef typename adjacency_list_BC< OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty>::adjacency_iterator child_vertex_iterator;
};





#define BGL_ADJACENCY_LIST_BC_ARGS typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties, typename EdgeProperties
#define BGL_ADJACENCY_LIST_BC adjacency_list_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>


template < BGL_ADJACENCY_LIST_BC_ARGS >
void swap(BGL_ADJACENCY_LIST_BC & lhs, BGL_ADJACENCY_LIST_BC & rhs) { 
  lhs.swap(rhs);
};



/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_descriptor
  source( const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, const BGL_ADJACENCY_LIST_BC &) {
  return e.source;
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_descriptor
  target( const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.get_stored_edge(e).target;
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC::out_edge_iterator,
 typename BGL_ADJACENCY_LIST_BC::out_edge_iterator >
  out_edges( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.out_edges(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::size_t
  out_degree( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC::in_edge_iterator,
 typename BGL_ADJACENCY_LIST_BC::in_edge_iterator >
  in_edges( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.in_edges(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::size_t in_degree( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.get_in_degree(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::size_t degree( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
};



/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC::vertex_iterator,
 typename BGL_ADJACENCY_LIST_BC::vertex_iterator > vertices( const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.vertices();
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertices_size_type num_vertices( const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.size();
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC::edge_iterator,
 typename BGL_ADJACENCY_LIST_BC::edge_iterator > edges( const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.edges();
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::edges_size_type num_edges( const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.num_edges();
};



/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair< 
typename BGL_ADJACENCY_LIST_BC::adjacency_iterator,
typename BGL_ADJACENCY_LIST_BC::adjacency_iterator >
  adjacent_vertices( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  typedef typename BGL_ADJACENCY_LIST_BC::adjacency_iterator AdjIter;
  typedef typename BGL_ADJACENCY_LIST_BC::out_edge_iterator OEIter;
  
  std::pair<OEIter, OEIter> oe_pair = out_edges(v, g);
  return std::pair<AdjIter, AdjIter>(AdjIter(oe_pair.first, &g), AdjIter(oe_pair.second, &g));
};



/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair< 
typename BGL_ADJACENCY_LIST_BC::inv_adjacency_iterator,
typename BGL_ADJACENCY_LIST_BC::inv_adjacency_iterator >
  inv_adjacent_vertices( typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  typedef typename BGL_ADJACENCY_LIST_BC::inv_adjacency_iterator InvAdjIter;
  typedef typename BGL_ADJACENCY_LIST_BC::in_edge_iterator IEIter;
  
  std::pair<IEIter, IEIter> ie_pair = in_edges(v, g);
  return std::pair<InvAdjIter, InvAdjIter>(InvAdjIter(ie_pair.first, &g), InvAdjIter(ie_pair.second, &g));
};



/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair< typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool > 
  edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u, 
       typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
       const BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.get_edge(u,v);
};





/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_descriptor add_vertex( BGL_ADJACENCY_LIST_BC & g) {
  typedef typename BGL_ADJACENCY_LIST_BC::vertex_bundled VertexBundled;
  return g.m_pack.add_vertex( VertexBundled() );
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void clear_vertex(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, BGL_ADJACENCY_LIST_BC & g) {
  g.m_pack.clear_vertex(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_vertex(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, BGL_ADJACENCY_LIST_BC & g) {
  g.m_pack.remove_vertex(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
           BGL_ADJACENCY_LIST_BC & g) {
  typedef typename BGL_ADJACENCY_LIST_BC::edge_bundled EdgeBundled;
  return g.m_pack.add_edge(u, v, EdgeBundled());
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u,
                 typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
                 BGL_ADJACENCY_LIST_BC & g) {
  std::pair< typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool > e_p = edge(u,v,g);
  if(e_p.second)
    g.m_pack.remove_edge(e_p.first);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, BGL_ADJACENCY_LIST_BC & g) {
  g.m_pack.remove_edge(e);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(typename BGL_ADJACENCY_LIST_BC::edge_iterator e_iter, BGL_ADJACENCY_LIST_BC & g) {
  g.m_pack.remove_edge(*e_iter);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_descriptor
  add_vertex(const typename BGL_ADJACENCY_LIST_BC::vertex_property_type& vp, BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.add_vertex( vp );
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_vertex(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
                   typename BGL_ADJACENCY_LIST_BC::vertex_property_type& vp,
                   BGL_ADJACENCY_LIST_BC & g) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  vp = std::move(g[v]);
#else
  vp = g[v];
#endif
  g.m_pack.remove_vertex(v);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
           const typename BGL_ADJACENCY_LIST_BC::edge_property_type& ep,
           BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.add_edge(u, v, ep);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u,
                 typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
                 typename BGL_ADJACENCY_LIST_BC::edge_property_type& ep,
                 BGL_ADJACENCY_LIST_BC & g) {
  std::pair< typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool > e_p = edge(u,v,g);
  if(e_p.second) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    ep = std::move(g[e_p.first]);
#else
    ep = g[e_p.first];
#endif
    g.m_pack.remove_edge(e_p.first);
  };
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e,
                 typename BGL_ADJACENCY_LIST_BC::edge_property_type& ep,
                 BGL_ADJACENCY_LIST_BC & g) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  ep = std::move(g[e]);
#else
  ep = g[e];
#endif
  g.m_pack.remove_edge(e);
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
void remove_edge(typename BGL_ADJACENCY_LIST_BC::edge_iterator e_iter,
                 typename BGL_ADJACENCY_LIST_BC::edge_property_type& ep,
                 BGL_ADJACENCY_LIST_BC& g) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  ep = std::move(g[*e_iter]);
#else
  ep = g[*e_iter];
#endif
  g.m_pack.remove_edge(*e_iter);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_descriptor
  add_vertex(typename BGL_ADJACENCY_LIST_BC::vertex_property_type&& vp,
             BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.add_vertex( std::move(vp) );
};

template < BGL_ADJACENCY_LIST_BC_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v,
           typename BGL_ADJACENCY_LIST_BC::edge_property_type&& ep,
           BGL_ADJACENCY_LIST_BC & g) {
  return g.m_pack.add_edge(u, v, std::move(ep));
};

#endif






/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which to draw the vertex.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
const typename BGL_ADJACENCY_LIST_BC::vertex_property_type& 
  get(const BGL_ADJACENCY_LIST_BC & g, typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v) {
  return g[v];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param g The tree from which to draw the edge.
  * \param e The edge descriptor of the sought-after edge-property.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
const typename BGL_ADJACENCY_LIST_BC::edge_property_type& 
  get( const BGL_ADJACENCY_LIST_BC & g, const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e) {
  return g[e];
};

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by const-reference, to be associated to the given vertex.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
void put( BGL_ADJACENCY_LIST_BC & g, 
          typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, 
          const typename BGL_ADJACENCY_LIST_BC::vertex_property_type& value) {
  g[v] = value;
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by const-reference, to be associated to the given edge.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
void put( BGL_ADJACENCY_LIST_BC & g, 
          const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, 
          const typename BGL_ADJACENCY_LIST_BC::edge_property_type& value) {
  g[e] = value;
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
void put( BGL_ADJACENCY_LIST_BC & g, 
          typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, 
          typename BGL_ADJACENCY_LIST_BC::vertex_property_type&& value) {
  g[v] = std::move(value);
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
void put( BGL_ADJACENCY_LIST_BC & g, 
          const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, 
          typename BGL_ADJACENCY_LIST_BC::edge_property_type&& value) {
  g[e] = std::move(value);
};

#endif


/**
  * Returns a reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by reference, associated to the given vertex descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::vertex_property_type& 
  get_property(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, BGL_ADJACENCY_LIST_BC & g) {
  return g[v];
};

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
const typename BGL_ADJACENCY_LIST_BC::vertex_property_type& 
  get_property(typename BGL_ADJACENCY_LIST_BC::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC & g) {
  return g[v];
};

/**
  * Returns a reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by reference, associated to the given edge descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
typename BGL_ADJACENCY_LIST_BC::edge_property_type& 
  get_property(const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, BGL_ADJACENCY_LIST_BC & g) {
  return g[e];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_ADJACENCY_LIST_BC_ARGS >
const typename BGL_ADJACENCY_LIST_BC::edge_property_type& 
  get_property(const typename BGL_ADJACENCY_LIST_BC::edge_descriptor& e, const BGL_ADJACENCY_LIST_BC & g) {
  return g[e];
};


#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle>
struct property_map< BGL_ADJACENCY_LIST_BC, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename BGL_ADJACENCY_LIST_BC::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, BGL_ADJACENCY_LIST_BC,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const BGL_ADJACENCY_LIST_BC,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle>
typename property_map< BGL_ADJACENCY_LIST_BC, T Bundle::* >::type
get( T Bundle::* p, BGL_ADJACENCY_LIST_BC& g) {
  return typename property_map< BGL_ADJACENCY_LIST_BC, T Bundle::* >::type(&g, p);
};

template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle>
typename property_map< BGL_ADJACENCY_LIST_BC, T Bundle::* >::const_type
get( T Bundle::* p, const BGL_ADJACENCY_LIST_BC& g) {
  return typename property_map< BGL_ADJACENCY_LIST_BC, T Bundle::* >::const_type(&g, p);
};



template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle, typename Key>
const typename remove_const< T >::type& get( T Bundle::* p, const BGL_ADJACENCY_LIST_BC& g, const Key& k) {
  return (g[k]).*p;
};

template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle, typename Key>
void put( T Bundle::* p, BGL_ADJACENCY_LIST_BC& g, const Key& k, const T& val) {
  (g[k]).*p = val;
};


#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_ADJACENCY_LIST_BC_ARGS, typename T, typename Bundle, typename Key>
void put( T Bundle::* p, BGL_ADJACENCY_LIST_BC& g, const Key& k, T&& val) {
  (g[k]).*p = std::move(val);
};

#endif



#endif


/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct adj_list_BC_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    typedef typename property_value<Property, Tag>::type value_type;
    
    typedef tagged_from_bundle_property_map<value_type, Graph, Tag> type;
    typedef tagged_from_bundle_property_map<const value_type, const Graph, Tag> const_type;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<adj_list_BC_tag> {
  typedef adj_list_BC_property_selector type;
};

template <>
struct edge_property_selector<adj_list_BC_tag> {
  typedef adj_list_BC_property_selector type;
};


template < BGL_ADJACENCY_LIST_BC_ARGS, typename Property>
typename property_map< BGL_ADJACENCY_LIST_BC, Property>::type
get(Property p, BGL_ADJACENCY_LIST_BC& g) {
  typedef typename property_map< BGL_ADJACENCY_LIST_BC, Property>::type Map;
  return Map(&g, p);
};

template < BGL_ADJACENCY_LIST_BC_ARGS, typename Property>
typename property_map< BGL_ADJACENCY_LIST_BC, Property>::const_type
get(Property p, const BGL_ADJACENCY_LIST_BC& g) {
  typedef typename property_map< BGL_ADJACENCY_LIST_BC, Property>::const_type Map;
  return Map(&g, p);
};


template < BGL_ADJACENCY_LIST_BC_ARGS, typename Property, typename Key>
typename property_map_value< BGL_ADJACENCY_LIST_BC, Property>::type
get(Property p, const BGL_ADJACENCY_LIST_BC& g, const Key& k) {
  return get_property_value(g[k], p);
};

template < BGL_ADJACENCY_LIST_BC_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_ADJACENCY_LIST_BC& g, const Key& k, const Value& val) {
  get_property_value(g[k], p) = val;
};


#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_ADJACENCY_LIST_BC_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_ADJACENCY_LIST_BC& g, const Key& k, Value&& val) {
  get_property_value(g[k], p) = std::move(val);
};

#endif



template < BGL_ADJACENCY_LIST_BC_ARGS >
void do_graph_deep_copy(BGL_ADJACENCY_LIST_BC& lhs, const BGL_ADJACENCY_LIST_BC& rhs) {
  typedef typename BGL_ADJACENCY_LIST_BC::vertex_descriptor Vertex;
  typedef typename BGL_ADJACENCY_LIST_BC::vertex_iterator VIter;
  typedef typename BGL_ADJACENCY_LIST_BC::edge_iterator EIter;
  
  lhs.m_pack.clear();
  
  // first, add all the vertices (keep unordered_map of correspondance).
  ::boost::unordered_map< Vertex, Vertex, graph::detail::BC_desc_hasher > v_map;
  VIter vi, vi_end;
  for(tie(vi, vi_end) = vertices(rhs); vi != vi_end; ++vi) {
    Vertex v = add_vertex(rhs[*vi], lhs);
    v_map[*vi] = v;
  };
  
  // then, go through all the edges and add them to the lhs:
  EIter ei, ei_end;
  for(tie(ei, ei_end) = edges(rhs); ei != ei_end; ++ei)
    add_edge(v_map[source(*ei, rhs)], v_map[target(*ei, rhs)], rhs[*ei], lhs);
};


#undef BGL_ADJACENCY_LIST_BC_ARGS
#undef BGL_ADJACENCY_LIST_BC




/**
 * This class implements an adjacency-list based on Boost.Containers that is tailored 
 * to store elements of an undirected graph.
 * \note This is a partial specialization for undirected edges.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS, typename VertexListS,
          typename VertexProperties, typename EdgeProperties>
class adjacency_list_BC<OutEdgeListS, VertexListS, undirectedS, VertexProperties, EdgeProperties>
{
  public:
    typedef typename adjacency_list_BC_disallowed_vertex_list<VertexListS>::type check_allowed_vertex_list;
    
    typedef adjacency_list_BC<OutEdgeListS, VertexListS, undirectedS, VertexProperties, EdgeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    typedef typename graph::detail::adjlistBC_vertex_container<
      VertexListS, OutEdgeListS, bidirectionalS, VertexProperties, EdgeProperties> bidir_storage_type;
    
    
    typedef typename bidir_storage_type::vertex_descriptor vertex_descriptor;
    typedef typename bidir_storage_type::vertices_size_type vertices_size_type;
    
    typedef typename bidir_storage_type::edge_descriptor bidir_edge_descriptor;
    typedef graph::detail::BC_undir_edge_desc<bidir_edge_descriptor> edge_descriptor;
    typedef typename bidir_storage_type::edges_size_type edges_size_type;
    typedef edges_size_type degree_size_type;
    
    
    typedef typename bidir_storage_type::vertex_iterator vertex_iterator;
    
    typedef typename bidir_storage_type::out_edge_iterator bidir_out_edge_iterator;
    typedef typename bidir_storage_type::in_edge_iterator bidir_in_edge_iterator;
    typedef typename bidir_storage_type::edge_iterator bidir_edge_iterator;
    
    typedef graph::detail::adjlistBC_undir_ioeiter<edge_descriptor, bidir_in_edge_iterator, bidir_out_edge_iterator> out_edge_iterator;
    typedef out_edge_iterator in_edge_iterator;
    
    typedef graph::detail::adjlistBC_undir_eiter<edge_descriptor, bidir_edge_iterator> edge_iterator;
    
    typedef typename adjacency_iterator_generator<self, vertex_descriptor, out_edge_iterator>::type adjacency_iterator;
    typedef typename inv_adjacency_iterator_generator<self, vertex_descriptor, in_edge_iterator>::type inv_adjacency_iterator;
    
    typedef typename bidir_storage_type::vertex_stored_type vertex_stored_impl;
    typedef typename bidir_storage_type::vertex_value_type vertex_value_impl;
    
    typedef typename bidir_storage_type::edge_stored_type edge_stored_impl;
    typedef typename bidir_storage_type::edge_value_type edge_value_impl;
    
    
    typedef adjacency_list_BC_traits<OutEdgeListS, VertexListS, undirectedS> Traits;
    
    typedef typename Traits::directed_category directed_category;
    typedef typename Traits::edge_parallel_category edge_parallel_category;
    typedef typename Traits::traversal_category traversal_category;
    
    
    typedef adj_list_BC_tag graph_tag;
    
    
    /**
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() { return graph::detail::BC_null_desc<vertex_descriptor>::value(); };
    
    /**
     * This static member function outputs the null-edge (invalid edge descriptor).
     * \return A null-edge descriptor (invalid edge descriptor).
     */
    static edge_descriptor null_edge() { return edge_descriptor(graph::detail::BC_null_desc<bidir_edge_descriptor>::value()); };
    
    
  public:  // private:
    
    bidir_storage_type m_pack;
    
  public:
    
    
    /**
     * Constructs an empty linked-tree.
     */
    adjacency_list_BC() : m_pack() { };
    
    
    adjacency_list_BC(const self& rhs) : m_pack() { do_graph_deep_copy(*this, rhs); };
    self& operator=(const self& rhs) {
      do_graph_deep_copy(*this, rhs);
      return *this;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    adjacency_list_BC(self&& rhs) : m_pack(std::move(rhs.m_pack)) { };
    self& operator=(self&& rhs) {
      m_pack = std::move(rhs.m_pack);
      return *this;
    };
#endif
    
    /**
     * Swaps the adjacency-list with another.
     */
    void swap(self& rhs) { 
      m_pack.swap(rhs.m_pack);
    };
    
    /**
     * Clears the adjacency-list of all vertices and edges.
     */
    void clear() { 
      m_pack.clear();
    };
    
    /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[](vertex_descriptor v) {
      return m_pack.get_stored_vertex(v).data;
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[](vertex_descriptor v) const {
      return m_pack.get_stored_vertex(v).data;
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[](const edge_descriptor& e) {
      return m_pack.get_stored_edge(bidir_edge_descriptor(e)).data;
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[](const edge_descriptor& e) const {
      return m_pack.get_stored_edge(bidir_edge_descriptor(e)).data;
    };
    
};



#define BGL_ADJACENCY_LIST_BC_UNDIR_ARGS typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties
#define BGL_ADJACENCY_LIST_BC_UNDIR adjacency_list_BC<OutEdgeListS, VertexListS, undirectedS, VertexProperties, EdgeProperties>




/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor
  source( const typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor& e, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  if( e.is_reversed )
    return g.m_pack.get_stored_edge(e).target;
  else
    return e.source;
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor
  target( const typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor& e, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  if( e.is_reversed )
    return e.source;
  else
    return g.m_pack.get_stored_edge(e).target;
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC_UNDIR::out_edge_iterator,
 typename BGL_ADJACENCY_LIST_BC_UNDIR::out_edge_iterator >
  out_edges( typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::out_edge_iterator EIter;
  
  std::pair<
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_out_edge_iterator,
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_out_edge_iterator > oe_pair = g.m_pack.out_edges(v);
  
  std::pair<
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_in_edge_iterator,
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_in_edge_iterator > ie_pair = g.m_pack.in_edges(v);
  
  return std::pair<EIter,EIter>(
    EIter(true, ie_pair.first, ie_pair.second, oe_pair.first, oe_pair.second),
    EIter(true, ie_pair.second, oe_pair.first, oe_pair.second)
  );
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::size_t
  out_degree( typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  return g.m_pack.get_out_degree(v) + g.m_pack.get_in_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC_UNDIR::in_edge_iterator,
 typename BGL_ADJACENCY_LIST_BC_UNDIR::in_edge_iterator >
  in_edges( typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::in_edge_iterator EIter;
  
  std::pair<
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_out_edge_iterator,
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_out_edge_iterator > oe_pair = g.m_pack.out_edges(v);
  
  std::pair<
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_in_edge_iterator,
    typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_in_edge_iterator > ie_pair = g.m_pack.in_edges(v);
  
  return std::pair<EIter,EIter>(
    EIter(false, ie_pair.first, ie_pair.second, oe_pair.first, oe_pair.second),
    EIter(false, ie_pair.second, oe_pair.first, oe_pair.second)
  );
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::size_t in_degree( typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::size_t degree( typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<
 typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_iterator,
 typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_iterator > edges( const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_iterator EIter;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_edge_iterator BiEIter;
  std::pair<BiEIter,BiEIter> be_pair = g.m_pack.edges();
  return std::pair<EIter,EIter>(EIter(be_pair.first), EIter(be_pair.second));
};

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
typename BGL_ADJACENCY_LIST_BC_UNDIR::edges_size_type num_edges( const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  return g.m_pack.num_edges();
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair< typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor, bool > 
  edge(typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor u, 
       typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v,
       const BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_edge_descriptor BiEdge;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor Edge;
  std::pair< BiEdge, bool > be_result = g.m_pack.get_edge(u,v);
  if( be_result.second) 
    return std::pair< Edge, bool >(Edge(be_result.first), true);
  else {
    be_result = g.m_pack.get_edge(v,u);
    if( be_result.second )
      return std::pair< Edge, bool >(Edge(be_result.first, true), true);
  };
  return std::pair< Edge, bool >(Edge(be_result.first), false);
};


/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v,
           BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_bundled EdgeBundled;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_edge_descriptor BiEdge;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor Edge;
  std::pair<BiEdge, bool> be_pair = g.m_pack.add_edge(u, v, EdgeBundled());
  return std::pair<Edge, bool>(Edge(be_pair.first), be_pair.second);
};

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v,
           const typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_property_type& ep,
           BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_edge_descriptor BiEdge;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor Edge;
  std::pair<BiEdge, bool> be_pair = g.m_pack.add_edge(u, v, ep);
  return std::pair<Edge, bool>(Edge(be_pair.first), be_pair.second);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_ADJACENCY_LIST_BC_UNDIR_ARGS >
std::pair<typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor, bool>
  add_edge(typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor u,
           typename BGL_ADJACENCY_LIST_BC_UNDIR::vertex_descriptor v,
           typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_property_type&& ep,
           BGL_ADJACENCY_LIST_BC_UNDIR & g) {
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::bidir_edge_descriptor BiEdge;
  typedef typename BGL_ADJACENCY_LIST_BC_UNDIR::edge_descriptor Edge;
  std::pair<BiEdge, bool> be_pair = g.m_pack.add_edge(u, v, std::move(ep));
  return std::pair<Edge, bool>(Edge(be_pair.first), be_pair.second);
};

#endif



#undef BGL_ADJACENCY_LIST_BC_UNDIR_ARGS
#undef BGL_ADJACENCY_LIST_BC_UNDIR





};


#endif


















