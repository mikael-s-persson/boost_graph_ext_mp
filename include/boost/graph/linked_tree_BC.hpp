// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file linked_tree_BC.hpp
 * 
 * This library provides a class that implements a (doubly-)linked tree structure. This is a 
 * classic tree implementation in which each node contain a list of edges to its children, 
 * and a link to its parent.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_LINKED_TREE_BC_HPP
#define BOOST_LINKED_TREE_BC_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/adjacency_iterator.hpp>
#include <boost/graph/detail/ltreeBC_containers.hpp>
#include <boost/graph/more_property_maps.hpp>

#include <boost/graph/tree_traits.hpp>

#include <utility>


namespace boost {


/**
 * This traits class template is used to obtain the types (and meta-values) that describe 
 * the basic types used in a linked-tree with the given out-edge-list and vertex-list storage 
 * policies. This traits class is useful to obtain type-erased (or type-agnostic) vertex and edge
 * descriptors. Note, this traits class is essentially the linked-tree equivalent of the BGL adjacency_list_traits class.
 */
template <typename OutEdgeListS = vecBC, 
          typename VertexListS = vecBC,
          typename DirectedS = directedS>
struct linked_tree_BC_traits {
  /** This meta-value tells if the edges are bidirectional, or not. */
  typedef typename DirectedS::is_bidir_t is_bidir;
  /** This meta-value tells if the edges are directional, or not. */
  typedef typename DirectedS::is_directed_t is_directed;
  
  /** This tag gives the edges' directional categorization. */
  typedef typename mpl::if_<is_bidir,
    bidirectional_tag,
    typename mpl::if_<is_directed,
      directed_tag, undirected_tag
    >::type
  >::type directed_category;
  
  /** This meta-value tells if the parallel edges are allowed, or not. */
  typedef typename detail::parallel_edge_BC_traits<OutEdgeListS>::type edge_parallel_category;
  
  typedef graph::detail::adjlistBC_traversal_tag<DirectedS> traversal_category;
  
  /** This meta-value tells if the vertex storage is random-access, or not. */
  typedef typename detail::is_random_access_BC<VertexListS>::type vertex_rand_access;
  typedef vertex_rand_access is_rand_access;
  
  /** This meta-value tells if the vertex storage is random-access, or not. */
  typedef typename detail::is_random_access_BC<OutEdgeListS>::type edge_rand_access;
  
  typedef std::size_t vertices_size_type;
  typedef std::size_t edges_size_type;
  
};





struct linked_tree_BC_tag { };

template <typename VertexListS>
struct linked_tree_BC_disallowed_vertex_list {
  typedef void type;
};

template <> struct linked_tree_BC_disallowed_vertex_list<setBC> { };
template <> struct linked_tree_BC_disallowed_vertex_list<multisetBC> { };
template <> struct linked_tree_BC_disallowed_vertex_list<unordered_setBC> { };
template <> struct linked_tree_BC_disallowed_vertex_list<unordered_multisetBC> { };



template <typename DirectedS>
struct linked_tree_BC_disallowed_undirected {
  typedef void type;
};

template <> struct linked_tree_BC_disallowed_undirected<undirectedS> { };



/**
 * This class implements a D-Ary Breadth-first tree that is tailored 
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
 * \tparam OutEdgeList A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexList A type tag to choose the storage policy for the vertices.
 * \tparam DirectedS A type tag to choose the directionality of the edges.
 * \tparam VertexProperties A POD type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A POD type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS = vecBC,
          typename VertexListS = vecBC,
          typename DirectedS = directedS,
          typename VertexProperties = no_property,
          typename EdgeProperties = no_property >
class linked_tree_BC
{
  public:
    typedef typename linked_tree_BC_disallowed_vertex_list<VertexListS>::type check_allowed_vertex_list;
    typedef typename linked_tree_BC_disallowed_undirected<DirectedS>::type check_allowed_directionality;
    
    typedef linked_tree_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    
    typedef typename graph::detail::ltreeBC_vertex_container<
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
    typedef adjacency_iterator child_vertex_iterator;
    
    
    typedef typename storage_type::vertex_stored_type vertex_stored_impl;
    typedef typename storage_type::vertex_value_type vertex_value_impl;
    
    typedef typename storage_type::edge_stored_type edge_stored_impl;
    typedef typename storage_type::edge_value_type edge_value_impl;
    
    
    typedef linked_tree_BC_traits<OutEdgeListS, VertexListS, DirectedS> Traits;
    
    typedef typename Traits::vertex_rand_access vertex_rand_access;
    typedef typename Traits::edge_rand_access edge_rand_access;
    typedef typename Traits::directed_category directed_category;
    typedef typename Traits::edge_parallel_category edge_parallel_category;
    typedef typename Traits::traversal_category traversal_category;
    
    typedef linked_tree_BC_tag graph_tag;
    
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
     * Constructs an empty linked-tree.
     */
    linked_tree_BC() : m_pack() { };
    
  private:
    void do_deep_copy_from(const self& rhs);
    
  public:
    
    /**
     * Constructs a linked-tree as a copy of the given tree.
     */
    linked_tree_BC(const self& rhs) : m_pack() { do_deep_copy_from(rhs); };
    /**
     * Assigns the linked-tree as a copy of the given tree.
     */
    self& operator=(const self& rhs) {
      do_deep_copy_from(rhs);
      return *this;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Constructs a linked-tree by moving the given tree into it.
     */
    linked_tree_BC(self&& rhs) : m_pack(std::move(rhs.m_pack)) { };
    /**
     * Assigns a linked-tree by moving the given tree into it.
     */
    linked_tree_BC& operator=(self&& rhs) {
      m_pack = std::move(rhs.m_pack);
      return *this;
    };
#endif
    
    /**
     * Standard swap function.
     */
    void swap(self& rhs) {
      m_pack.swap(rhs.m_pack);
    };
    
    
    /**
     * Returns the size of the tree (the number of vertices it contains).
     * \return The size of the tree (the number of vertices it contains).
     */
    std::size_t size() const { return m_pack.size(); };
    
    /**
     * Checks if the tree is empty.
     * \return True if the tree is empty.
     */
    bool empty() const { return size() == 0; };
    
    /**
     * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
     * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
     */
    std::size_t capacity() const { return m_pack.m_vertices.capacity(); };
    
    /**
     * Returns the depth of the tree.
     * \note This operation must recurse through all the branches of the tree (depth-first), and is 
     * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack) 
     * w.r.t. the depth of tree).
     * \return The depth of the tree.
     */
    std::size_t depth() const { 
      if(m_pack.m_root != null_vertex())
        return m_pack.get_depth( m_pack.m_root );
      return 0;
    };
    
    /**
     * Clears the tree of all vertices and edges.
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



/**
 * This is the tree-storage specifier for a linked-tree of the given edge and vertex storage policies.
 */
template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct linked_tree_BC_storage { };


template <typename VertexProperties, typename EdgeProperties, typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage<VertexProperties, EdgeProperties, linked_tree_BC_storage<OutEdgeListS, VertexListS, DirectedS> > {
  typedef linked_tree_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties> type;
};


template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage_traits< linked_tree_BC_storage<OutEdgeListS, VertexListS, DirectedS> > :
  linked_tree_BC_traits<OutEdgeListS, VertexListS, DirectedS> { };



#define BGL_LINKED_TREE_BC_ARGS typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties, typename EdgeProperties
#define BGL_LINKED_TREE_BC linked_tree_BC<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>

/**
 * Standard swap function. Swaps the contents of two objects.
 * \param lhs The left-hand-side of the swap.
 * \param rhs The right-hand-side of the swap.
 */
template < BGL_LINKED_TREE_BC_ARGS >
void swap(BGL_LINKED_TREE_BC& lhs, BGL_LINKED_TREE_BC& rhs) { lhs.swap(rhs); };


/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

/**
 * Returns the source vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The source vertex of the given edge descriptor.
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_descriptor
  source( typename BGL_LINKED_TREE_BC::edge_descriptor e, const BGL_LINKED_TREE_BC &) {
  return e.source;
};

/**
 * Returns the target vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The target vertex of the given edge descriptor.
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_descriptor
  target( typename BGL_LINKED_TREE_BC::edge_descriptor e, const BGL_LINKED_TREE_BC& g) {
  return g.m_pack.get_stored_edge(e).target;
};

/**
 * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the out-edges of a given vertex descriptor.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair<
 typename BGL_LINKED_TREE_BC::out_edge_iterator,
 typename BGL_LINKED_TREE_BC::out_edge_iterator >
  out_edges( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.out_edges(v);
};

/**
 * Returns the out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The out-degree of the given vertex descriptor.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::size_t
  out_degree( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for the in-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the in-edges of a given vertex descriptor.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair<
 typename BGL_LINKED_TREE_BC::in_edge_iterator,
 typename BGL_LINKED_TREE_BC::in_edge_iterator >
  in_edges( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.in_edges(v);
};

/**
 * Returns the in-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::size_t in_degree( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.get_in_degree(v);
};

/**
 * Returns the in-degree plus out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree plus out-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::size_t degree( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
};



/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the vertices of the tree.
 * \param g The graph.
 * \return The vertex iterator range for all the vertices of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair<
 typename BGL_LINKED_TREE_BC::vertex_iterator,
 typename BGL_LINKED_TREE_BC::vertex_iterator > vertices( const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.vertices();
};

/**
 * Returns the size of the tree (the number of vertices it contains).
 * \param g The graph.
 * \return The size of the tree (the number of vertices it contains).
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertices_size_type num_vertices( const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.size();
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for all the edges of the tree.
 * \param g The graph.
 * \return The edge iterator range for all the edges of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair<
 typename BGL_LINKED_TREE_BC::edge_iterator,
 typename BGL_LINKED_TREE_BC::edge_iterator > edges( const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.edges();
};

/**
 * Returns the number of edges in the tree.
 * \param g The graph.
 * \return The number of edges in the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::edges_size_type num_edges( const BGL_LINKED_TREE_BC & g) {
  std::size_t tmp = g.m_pack.size();
  if(tmp > 0)
    return tmp - 1;
  else
    return 0;
};



/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_TREE_BC_ARGS >
std::pair< 
typename BGL_LINKED_TREE_BC::adjacency_iterator,
typename BGL_LINKED_TREE_BC::adjacency_iterator >
  adjacent_vertices( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  typedef typename BGL_LINKED_TREE_BC::adjacency_iterator AdjIter;
  typedef typename BGL_LINKED_TREE_BC::out_edge_iterator OEIter;
  
  std::pair<OEIter, OEIter> oe_pair = out_edges(v, g);
  return std::pair<AdjIter, AdjIter>(AdjIter(oe_pair.first, &g), AdjIter(oe_pair.second, &g));
};



/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_TREE_BC_ARGS >
std::pair< 
typename BGL_LINKED_TREE_BC::inv_adjacency_iterator,
typename BGL_LINKED_TREE_BC::inv_adjacency_iterator >
  inv_adjacent_vertices( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  typedef typename BGL_LINKED_TREE_BC::inv_adjacency_iterator InvAdjIter;
  typedef typename BGL_LINKED_TREE_BC::in_edge_iterator IEIter;
  
  std::pair<IEIter, IEIter> ie_pair = in_edges(v, g);
  return std::pair<InvAdjIter, InvAdjIter>(InvAdjIter(ie_pair.first, &g), InvAdjIter(ie_pair.second, &g));
};






/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

/**
 * Returns the edge descriptor for the edge between two given vertex descriptors.
 * \param u The vertex descriptor of the source vertex.
 * \param v The vertex descriptor of the target vertex.
 * \param g The graph.
 * \return The edge descriptor for the given vertex descriptor pair.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair<
  typename BGL_LINKED_TREE_BC::edge_descriptor,
  bool > edge( typename BGL_LINKED_TREE_BC::vertex_descriptor u, 
               typename BGL_LINKED_TREE_BC::vertex_descriptor v,
               const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.get_edge(u,v);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_descriptor get_root_vertex( const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.m_root;
};

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair< 
typename BGL_LINKED_TREE_BC::child_vertex_iterator,
typename BGL_LINKED_TREE_BC::child_vertex_iterator >
  child_vertices( typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return adjacent_vertices(v,g);
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
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_descriptor
  parent_vertex(typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g.m_pack.get_parent(v);
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_descriptor create_root( BGL_LINKED_TREE_BC& g) {
  typedef typename BGL_LINKED_TREE_BC::vertex_property_type VProp;
  if(g.m_pack.size())
    g.m_pack.clear();
  g.m_pack.add_root_vertex(VProp());
  return g.m_pack.m_root;
};

/**
 * Adds a child vertex to the given parent vertex, and default-initializes the properties of 
 * the newly created vertex and edge.
 * \param v The parent vertex to which a child will be added.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template < BGL_LINKED_TREE_BC_ARGS >
std::pair< 
typename BGL_LINKED_TREE_BC::vertex_descriptor,
typename BGL_LINKED_TREE_BC::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_TREE_BC::vertex_descriptor v, BGL_LINKED_TREE_BC & g) {
  typedef typename BGL_LINKED_TREE_BC::vertex_property_type VProp;
  typedef typename BGL_LINKED_TREE_BC::edge_property_type EProp;
  return g.m_pack.add_child(v, VProp(), EProp());
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template < BGL_LINKED_TREE_BC_ARGS >
void clear_children( typename BGL_LINKED_TREE_BC::vertex_descriptor v, BGL_LINKED_TREE_BC & g) {
  g.m_pack.clear_children_impl(v);
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template < BGL_LINKED_TREE_BC_ARGS >
void remove_branch( typename BGL_LINKED_TREE_BC::vertex_descriptor v, BGL_LINKED_TREE_BC & g) {
  g.m_pack.remove_branch_impl(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS, typename VProp >
typename BGL_LINKED_TREE_BC::vertex_descriptor create_root( const VProp& vp,  BGL_LINKED_TREE_BC & g) {
  if(g.m_pack.size())
    g.m_pack.clear();
  g.m_pack.add_root_vertex(vp);
  return g.m_pack.m_root;
};
#else
/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_LINKED_TREE_BC_ARGS, typename VProp >
typename BGL_LINKED_TREE_BC::vertex_descriptor create_root( VProp&& vp, BGL_LINKED_TREE_BC & g) {
  if(g.m_pack.size())
    g.m_pack.clear();
  g.m_pack.add_root_vertex(std::forward<VProp>(vp));
  return g.m_pack.m_root;
};
#endif


#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
 * vertex to the given property value.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template < BGL_LINKED_TREE_BC_ARGS, typename VProp >
std::pair< 
typename BGL_LINKED_TREE_BC::vertex_descriptor,
typename BGL_LINKED_TREE_BC::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_TREE_BC::vertex_descriptor v,
                    const VProp& vp, BGL_LINKED_TREE_BC & g) {
  typedef typename BGL_LINKED_TREE_BC::edge_property_type EProp;
  return g.m_pack.add_child(v, vp, EProp());
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
template < BGL_LINKED_TREE_BC_ARGS, typename VProp, typename EProp >
std::pair< 
typename BGL_LINKED_TREE_BC::vertex_descriptor,
typename BGL_LINKED_TREE_BC::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_TREE_BC::vertex_descriptor v,
                    const VProp& vp, const EProp& ep, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.add_child(v, vp, ep);
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
template < BGL_LINKED_TREE_BC_ARGS, typename VProp >
std::pair< 
typename BGL_LINKED_TREE_BC::vertex_descriptor,
typename BGL_LINKED_TREE_BC::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_TREE_BC::vertex_descriptor v,
                    VProp&& vp, BGL_LINKED_TREE_BC & g) {
  typedef typename BGL_LINKED_TREE_BC::edge_property_type EProp;
  return g.m_pack.add_child(v, std::forward<VProp>(vp), EProp());
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
template < BGL_LINKED_TREE_BC_ARGS, typename VProp, typename EProp >
std::pair< 
typename BGL_LINKED_TREE_BC::vertex_descriptor,
typename BGL_LINKED_TREE_BC::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_TREE_BC::vertex_descriptor v,
                    VProp&& vp, EProp&& ep, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.add_child(v, std::forward<VProp>(vp), std::forward<EProp>(ep));
};
#endif


/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template < BGL_LINKED_TREE_BC_ARGS , typename OutputIter>
OutputIter clear_children( typename BGL_LINKED_TREE_BC::vertex_descriptor v, OutputIter it_out, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.clear_children_impl(v, it_out);
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
template < BGL_LINKED_TREE_BC_ARGS , typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children( typename BGL_LINKED_TREE_BC::vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.clear_children_impl(v, vit_out, eit_out);
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
template < BGL_LINKED_TREE_BC_ARGS , typename OutputIter>
OutputIter remove_branch( typename BGL_LINKED_TREE_BC::vertex_descriptor v, OutputIter it_out, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.remove_branch_impl(v, it_out);
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
template < BGL_LINKED_TREE_BC_ARGS , typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch( typename BGL_LINKED_TREE_BC::vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out, BGL_LINKED_TREE_BC & g) {
  return g.m_pack.remove_branch_impl(v, vit_out, eit_out);
};




/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which to draw the vertex.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
const typename BGL_LINKED_TREE_BC::vertex_property_type& get(const BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::vertex_descriptor v) {
  return g[v];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param g The tree from which to draw the edge.
  * \param e The edge descriptor of the sought-after edge-property.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
const typename BGL_LINKED_TREE_BC::edge_property_type& get( const BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::edge_descriptor e) {
  return g[e];
};

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by const-reference, to be associated to the given vertex.
  */
template < BGL_LINKED_TREE_BC_ARGS >
void put( BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::vertex_descriptor v, const typename BGL_LINKED_TREE_BC::vertex_property_type& value) {
  g[v] = value;
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by const-reference, to be associated to the given edge.
  */
template < BGL_LINKED_TREE_BC_ARGS >
void put( BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::edge_descriptor e, const typename BGL_LINKED_TREE_BC::edge_property_type& value) {
  g[e] = value;
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
  */
template < BGL_LINKED_TREE_BC_ARGS >
void put( BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::vertex_descriptor v, typename BGL_LINKED_TREE_BC::vertex_property_type&& value) {
  g[v] = std::move(value);
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
  */
template < BGL_LINKED_TREE_BC_ARGS >
void put( BGL_LINKED_TREE_BC & g, typename BGL_LINKED_TREE_BC::edge_descriptor e, typename BGL_LINKED_TREE_BC::edge_property_type&& value) {
  g[e] = std::move(value);
};

#endif


/**
  * Returns a reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by reference, associated to the given vertex descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::vertex_property_type& get_property(typename BGL_LINKED_TREE_BC::vertex_descriptor v, BGL_LINKED_TREE_BC & g) {
  return g[v];
};

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
const typename BGL_LINKED_TREE_BC::vertex_property_type& get_property(typename BGL_LINKED_TREE_BC::vertex_descriptor v, const BGL_LINKED_TREE_BC & g) {
  return g[v];
};

/**
  * Returns a reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
typename BGL_LINKED_TREE_BC::edge_property_type& get_property(typename BGL_LINKED_TREE_BC::edge_descriptor e, BGL_LINKED_TREE_BC & g) {
  return g[e];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_TREE_BC_ARGS >
const typename BGL_LINKED_TREE_BC::edge_property_type& get_property(typename BGL_LINKED_TREE_BC::edge_descriptor e, const BGL_LINKED_TREE_BC & g) {
  return g[e];
};


#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template < BGL_LINKED_TREE_BC_ARGS, typename T, typename Bundle>
struct property_map< BGL_LINKED_TREE_BC, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename BGL_LINKED_TREE_BC::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, BGL_LINKED_TREE_BC,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const BGL_LINKED_TREE_BC,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


template < BGL_LINKED_TREE_BC_ARGS, typename T, typename Bundle>
typename property_map< BGL_LINKED_TREE_BC, T Bundle::* >::type
get( T Bundle::* p, BGL_LINKED_TREE_BC& g) {
  return typename property_map< BGL_LINKED_TREE_BC, T Bundle::* >::type(&g, p);
};

template < BGL_LINKED_TREE_BC_ARGS, typename T, typename Bundle>
typename property_map< BGL_LINKED_TREE_BC, T Bundle::* >::const_type
get( T Bundle::* p, const BGL_LINKED_TREE_BC& g) {
  return typename property_map< BGL_LINKED_TREE_BC, T Bundle::* >::const_type(&g, p);
};

#endif


/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/


struct linked_tree_BC_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    typedef typename property_value<Property, Tag>::type value_type;
    
    typedef tagged_from_bundle_property_map<value_type, Graph, Tag> type;
    typedef tagged_from_bundle_property_map<const value_type, const Graph, Tag> const_type;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<linked_tree_BC_tag> {
  typedef linked_tree_BC_property_selector type;
};

template <>
struct edge_property_selector<linked_tree_BC_tag> {
  typedef linked_tree_BC_property_selector type;
};



template < BGL_LINKED_TREE_BC_ARGS, typename Property>
typename property_map< BGL_LINKED_TREE_BC, Property>::type
get(Property p, BGL_LINKED_TREE_BC& g) {
  typedef typename property_map< BGL_LINKED_TREE_BC, Property>::type Map;
  return Map(&g, p);
};

template < BGL_LINKED_TREE_BC_ARGS, typename Property>
typename property_map< BGL_LINKED_TREE_BC, Property>::const_type
get(Property p, const BGL_LINKED_TREE_BC& g) {
  typedef typename property_map< BGL_LINKED_TREE_BC, Property>::const_type Map;
  return Map(&g, p);
};


template < BGL_LINKED_TREE_BC_ARGS, typename Property, typename Key>
typename property_map_value< BGL_LINKED_TREE_BC, Property>::type
get(Property p, const BGL_LINKED_TREE_BC& g, const Key& k) {
  return get_property_value(g[k], p);
};

template < BGL_LINKED_TREE_BC_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_LINKED_TREE_BC& g, const Key& k, const Value& val) {
  get_property_value(g[k], p) = val;
};


#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_LINKED_TREE_BC_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_LINKED_TREE_BC& g, const Key& k, Value&& val) {
  get_property_value(g[k], p) = std::move(val);
};

#endif








template < BGL_LINKED_TREE_BC_ARGS >
void BGL_LINKED_TREE_BC::do_deep_copy_from(const BGL_LINKED_TREE_BC& rhs) {
  typedef typename BGL_LINKED_TREE_BC::vertex_descriptor Vertex;
  typedef typename BGL_LINKED_TREE_BC::edge_descriptor Edge;
  typedef typename BGL_LINKED_TREE_BC::out_edge_iterator OEIter;
  typedef std::pair< Vertex, Vertex > TaskType;
  
  this->m_pack.clear();
  
  // NOTE: Do the copy in breadth-first order because it leads to the best memory traversal patterns later on.
  std::queue< TaskType > bft_queue;
  TaskType cur;
  cur.first = get_root_vertex(rhs);
  cur.second = create_root(rhs[cur.first], *this);
  bft_queue.push( cur );
  
  while( !bft_queue.empty() ) {
    cur = bft_queue.front(); bft_queue.pop();
    OEIter ei, ei_end;
    for(tie(ei, ei_end) = out_edges(cur.first, rhs); ei != ei_end; ++ei) {
      std::pair< Vertex, Edge > tmp = add_child_vertex(cur.second, rhs[target(*ei, rhs)], rhs[*ei], *this);
      bft_queue.push( TaskType(target(*ei, rhs), tmp.first) );
    };
  };
  
};






#undef BGL_LINKED_TREE_BC_ARGS
#undef BGL_LINKED_TREE_BC



};


#endif


















