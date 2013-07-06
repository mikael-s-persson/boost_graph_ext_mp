// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file linked_tree.hpp
 * 
 * This library provides a class that implements a (doubly-)linked tree structure. This is a 
 * classic tree implementation in which each node contain a list of edges to its children, 
 * and a link to its parent.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_LINKED_TREE_HPP
#define BOOST_LINKED_TREE_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/detail/linked_tree_iterators.hpp>
#include <boost/graph/detail/linked_tree_containers.hpp>
#include <boost/graph/more_property_maps.hpp>

#include <utility>


namespace boost {


/**
 * This traits class template is used to obtain the types (and meta-values) that describe 
 * the basic types used in a linked-tree with the given out-edge-list and vertex-list storage 
 * policies. This traits class is useful to obtain type-erased (or type-agnostic) vertex and edge
 * descriptors. Note, this traits class is essentially the linked-tree equivalent of the BGL adjacency_list_traits class.
 */
template <typename OutEdgeListS = vecS, 
          typename VertexListS = vecS>
struct linked_tree_traits {
  /** This meta-value tells if the edges are bidirectional, or not. */
  typedef mpl::true_ is_bidir;
  /** This meta-value tells if the edges are directional, or not. */
  typedef mpl::true_ is_directed;
  
  /** This tag gives the edges' directional categorization. */
  typedef bidirectional_tag directed_category;
  
  /** This meta-value tells if the parallel edges are allowed, or not. */
  typedef disallow_parallel_edge_tag edge_parallel_category;
  
  /** This meta-value tells if the vertex storage is random-access, or not. */
  typedef typename detail::is_random_access<VertexListS>::type vertex_rand_access;
  
  /** This meta-value tells if the vertex storage is random-access, or not. */
  typedef typename detail::is_random_access<OutEdgeListS>::type edge_rand_access;
  
};



/**
 * This class implements a D-Ary Breadth-first tree that is tailored 
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
 * \tparam OutEdgeList A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexList A type tag to choose the storage policy for the vertices.
 * \tparam VertexProperties A POD type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A POD type to be attached to each edge in the tree.
 * \tparam GraphProperties A POD type to be attached to each tree object.
 */
template <typename OutEdgeListS = vecS,
          typename VertexListS = vecS,
          typename VertexProperties = no_property,
          typename EdgeProperties = no_property,
          typename TreeProperties = no_property >
class linked_tree
{
  public:
    typedef linked_tree<OutEdgeListS, VertexListS, VertexProperties, EdgeProperties, TreeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    
    typedef typename graph::detail::ltree_container_select<VertexListS, OutEdgeListS, VertexProperties, EdgeProperties>::type storage_type;
    
    typedef typename storage_type::vertex_descriptor vertex_descriptor;
    typedef typename storage_type::vertices_size_type vertices_size_type;
    
    typedef typename storage_type::edge_descriptor edge_descriptor;
    typedef typename storage_type::edges_size_type edges_size_type;
    typedef edges_size_type degree_size_type;
    
    typedef linked_tree_traits<OutEdgeListS, VertexListS> ltree_traits;
    
    typedef typename ltree_traits::vertex_rand_access vertex_rand_access;
    typedef typename ltree_traits::edge_rand_access edge_rand_access;
    
    typedef typename ltree_traits::directed_category directed_category;
    typedef typename ltree_traits::edge_parallel_category edge_parallel_category;
    struct traversal_category : 
      bidirectional_graph_tag, 
      adjacency_graph_tag,
      vertex_list_graph_tag,
      edge_list_graph_tag,
      adjacency_matrix_tag { };
    
    /**
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() { return graph::detail::ltree_null_desc<vertex_descriptor>::value(); };
    
    /**
     * This static member function outputs the null-edge (invalid edge descriptor).
     * \return A null-edge descriptor (invalid edge descriptor).
     */
    static edge_descriptor null_edge() { return graph::detail::ltree_null_desc<edge_descriptor>::value(); };
    
    typedef typename storage_type::edge_value_type edge_value_type;
    typedef typename storage_type::vertex_value_type vertex_value_type;
    
    
    typedef typename storage_type::vertex_iterator vertex_iterator;
    typedef typename storage_type::child_vertex_iterator child_vertex_iterator;
    typedef typename storage_type::adjacency_iterator adjacency_iterator;
    typedef typename storage_type::out_edge_iterator out_edge_iterator;
    typedef typename storage_type::in_edge_iterator in_edge_iterator;
    typedef typename storage_type::edge_iterator edge_iterator;
    
  private:
    
    storage_type m_pack;
    
  public:
    
    
    /**
     * Constructs an empty linked-tree.
     */
    linked_tree() : m_pack() { };
    
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
    
    /**
     * Returns the out-degree of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The out-degree of the given vertex descriptor.
     */
    edges_size_type get_out_degree(vertex_descriptor v) const {
      return m_pack.get_stored_vertex(v).out_edges.size();
    };
    
    /**
     * Returns the in-degree of a given vertex descriptor in the tree.
     * \param v_i The vertex descriptor.
     * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
     */
    edges_size_type get_in_degree(vertex_descriptor v_i) const {
      if((v_i != null_vertex()) && (v_i != m_pack.m_root))
        return 1;
      return 0;
    };
    
    /**
     * Returns the parent vertex of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The parent vertex of the given vertex descriptor (will be null_vertex() if it is the root (no parent)).
     */
    vertex_descriptor get_parent_vertex(vertex_descriptor v) const {
      m_pack.get_stored_vertex(v).in_edge.source;
    };
    
    /**
     * Returns the target vertex of a given edge descriptor in the tree.
     * \param e The edge descriptor.
     * \return The target vertex of the given edge descriptor.
     */
    vertex_descriptor get_edge_target(const edge_descriptor& e) const {
      return m_pack.get_stored_edge(e).target;
    };
    
    /**
     * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The edge iterator range for the out-edges of a given vertex descriptor.
     */
    std::pair< out_edge_iterator, out_edge_iterator > out_edges(vertex_descriptor v) {
      return m_pack.out_edges(v);
    };
    
    /**
     * Returns the edge iterator range for the in-edges of a given vertex descriptor in the tree.
     * \param v The vertex descriptor.
     * \return The edge iterator range for the in-edges of a given vertex descriptor.
     */
    std::pair< in_edge_iterator, in_edge_iterator > in_edges( vertex_descriptor v) {
      return m_pack.in_edges(v);
    };
    
    /**
     * Returns the vertex iterator range for all the vertices of the tree.
     * \return The vertex iterator range for all the vertices of the tree.
     */
    std::pair< vertex_iterator, vertex_iterator > vertices() {
      return m_pack.vertices();
    };
    
    /**
     * Returns the edge iterator range for all the edges of the tree.
     * \return The edge iterator range for all the edges of the tree.
     */
    std::pair< edge_iterator, edge_iterator > edges() const {
      return m_pack.edges();
    };
    
    /**
     * Returns the edge descriptor for the edge between two given vertex descriptors.
     * \param u The vertex descriptor of the source vertex.
     * \param v The vertex descriptor of the target vertex.
     * \return The edge descriptor for the given vertex descriptor pair.
     */
    std::pair<edge_descriptor,bool> get_edge( vertex_descriptor u, vertex_descriptor v) {
      return m_pack.get_edge(u,v);
    };
    
    /**
     * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
     * \param v The vertex descriptor whose children are sought.
     * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
     */
    std::pair< child_vertex_iterator, child_vertex_iterator > child_vertices(vertex_descriptor v) {
      return m_pack.child_vertices(v);
    };
    
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values.
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value for the newly created vertex.
     * \param ep The property value for the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
    std::pair< vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, 
                                                             const vertex_property_type& vp = vertex_property_type(), 
                                                             const edge_property_type& ep = edge_property_type()) {
      return m_pack.add_child(v, vp, ep);
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values, by move-semantics (C++11).
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value to be moved into the newly created vertex.
     * \param ep The property value to be moved into the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
                                                             vertex_property_type&& vp, 
                                                             edge_property_type&& ep = edge_property_type()) {
      return m_pack.add_child(v, std::move(vp), std::move(ep));
    };
#endif
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
    void remove_branch(vertex_descriptor v) {
      m_pack.remove_branch_impl(v);
    };
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex, while 
     * recording the vertex-properties of all the removed vertices into an output-iterator.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \return The output-iterator after the collection of all the removed vertices.
     * \note The first vertex-property to figure in the output range is that of the vertex v.
     */
    template <typename OutputIter>
    OutputIter remove_branch(vertex_descriptor v, OutputIter it_out) {
      return m_pack.remove_branch_impl(v, it_out);
    };
    
    /**
     * Returns the vertex-descriptor of the root of the tree.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor get_root_vertex() const {
      return m_pack.m_root; 
    };
    
    /**
     * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
     * \param vp The vertex-property to assign to the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(m_pack.size())
        m_pack.clear();
      m_pack.add_root_vertex(vp);
      return m_pack.m_root;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
     * \param vp The vertex-property to move into the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(vertex_property_type&& vp) {
      if(m_pack.size())
        m_pack.clear();
      m_pack.add_root_vertex(std::move(vp));
      return m_pack.m_root;
    };
#endif
    
    
    
};



#define BGL_LINKED_LIST_ARGS typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties, typename TreeProperties
#define BGL_LINKED_LIST linked_tree<OutEdgeListS, VertexListS,VertexProperties, EdgeProperties, TreeProperties>




/**
 * This is the tree-storage specifier for a linked-tree of the given edge and vertex storage policies.
 */
template <typename OutEdgeListS, typename VertexListS>
struct linked_tree_storage { };


template <typename VertexProperties, typename EdgeProperties, typename OutEdgeListS, typename VertexListS>
struct tree_storage<VertexProperties, EdgeProperties, linked_tree_storage<OutEdgeListS, VertexListS> > {
  typedef linked_tree<OutEdgeListS, VertexListS, VertexProperties, EdgeProperties> type;
};


template <typename OutEdgeListS, typename VertexListS>
struct tree_storage_traits< linked_tree_storage<OutEdgeListS, VertexListS> > :
  linked_tree_traits<OutEdgeListS, VertexListS> { };





/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor
  source( typename BGL_LINKED_LIST::edge_descriptor e, const BGL_LINKED_LIST &) {
  return e.source;
};

template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor
  target( typename BGL_LINKED_LIST::edge_descriptor e, const BGL_LINKED_LIST& g) {
  return g.get_edge_target(e);
};

template < BGL_LINKED_LIST_ARGS >
std::pair<
 typename BGL_LINKED_LIST::out_edge_iterator,
 typename BGL_LINKED_LIST::out_edge_iterator >
  out_edges( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.out_edges(v);
};

template < BGL_LINKED_LIST_ARGS >
std::size_t
  out_degree( typename BGL_LINKED_LIST::vertex_descriptor v, const BGL_LINKED_LIST & g) {
  return g.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
std::pair<
 typename BGL_LINKED_LIST::in_edge_iterator,
 typename BGL_LINKED_LIST::in_edge_iterator >
  in_edges( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.in_edges(v);
};

template < BGL_LINKED_LIST_ARGS >
std::size_t in_degree( typename BGL_LINKED_LIST::vertex_descriptor v, const BGL_LINKED_LIST & g) {
  return g.get_in_degree(v);
};

template < BGL_LINKED_LIST_ARGS >
std::size_t degree( typename BGL_LINKED_LIST::vertex_descriptor v, const BGL_LINKED_LIST & g) {
  return g.get_in_degree(v) + g.get_out_degree(v);
};



/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
std::pair<
 typename BGL_LINKED_LIST::vertex_iterator,
 typename BGL_LINKED_LIST::vertex_iterator > vertices( BGL_LINKED_LIST & g) {
  return g.vertices();
};

template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertices_size_type num_vertices( const BGL_LINKED_LIST & g) {
  return g.size();
};


/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
std::pair<
 typename BGL_LINKED_LIST::edge_iterator,
 typename BGL_LINKED_LIST::edge_iterator > edges( BGL_LINKED_LIST & g) {
  return g.edges();
};

template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::edges_size_type num_edges( const BGL_LINKED_LIST & g) {
  std::size_t tmp = g.size();
  if(tmp > 0)
    return tmp - 1;
  else
    return 0;
};



/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::adjacency_iterator,
typename BGL_LINKED_LIST::adjacency_iterator >
  adjacent_vertices( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.child_vertices(v);
};


/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
std::pair<
  typename BGL_LINKED_LIST::edge_descriptor,
  bool > edge( typename BGL_LINKED_LIST::vertex_descriptor u, 
               typename BGL_LINKED_LIST::vertex_descriptor v,
               BGL_LINKED_LIST & g) {
  return g.get_edge(u,v);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor get_root_vertex( const BGL_LINKED_LIST & g) {
  return g.get_root_vertex();
};

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::child_vertex_iterator,
typename BGL_LINKED_LIST::child_vertex_iterator >
  child_vertices( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.child_vertices(v);
};


/***********************************************************************************************
 *                             BidirectionalTreeConcept
 * ********************************************************************************************/

template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor
  parent_vertex(typename BGL_LINKED_LIST::vertex_descriptor v, const BGL_LINKED_LIST & g) {
  return g.get_parent_vertex(v);
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor create_root( BGL_LINKED_LIST& g) {
  return g.create_root_vertex();
};

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::vertex_descriptor,
typename BGL_LINKED_LIST::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.add_child(v);
};

template < BGL_LINKED_LIST_ARGS >
void remove_branch( typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g.remove_branch(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor
  create_root( const typename BGL_LINKED_LIST::vertex_property_type& vp,  BGL_LINKED_LIST & g) {
  return g.create_root_vertex(vp);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_descriptor
  create_root( typename BGL_LINKED_LIST::vertex_property_type&& vp, BGL_LINKED_LIST & g) {
  return g.create_root_vertex(std::move(vp));
};
#endif

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::vertex_descriptor,
typename BGL_LINKED_LIST::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_LIST::vertex_descriptor v,
                    const typename BGL_LINKED_LIST::vertex_property_type& vp,
                    BGL_LINKED_LIST & g) {
  return g.add_child(v,vp);
};

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::vertex_descriptor,
typename BGL_LINKED_LIST::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_LIST::vertex_descriptor v,
                    const typename BGL_LINKED_LIST::vertex_property_type& vp,
                    const typename BGL_LINKED_LIST::edge_property_type& ep,
                    BGL_LINKED_LIST & g) {
  return g.add_child(v,vp,ep);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::vertex_descriptor,
typename BGL_LINKED_LIST::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_LIST::vertex_descriptor v,
                    typename BGL_LINKED_LIST::vertex_property_type&& vp,
                    BGL_LINKED_LIST & g) {
  return g.add_child(v,std::move(vp));
};

template < BGL_LINKED_LIST_ARGS >
std::pair< 
typename BGL_LINKED_LIST::vertex_descriptor,
typename BGL_LINKED_LIST::edge_descriptor >
  add_child_vertex( typename BGL_LINKED_LIST::vertex_descriptor v,
                    typename BGL_LINKED_LIST::vertex_property_type&& vp,
                    typename BGL_LINKED_LIST::edge_property_type&& ep,
                    BGL_LINKED_LIST & g) {
  return g.add_child(v,std::move(vp),std::move(ep));
};
#endif

template < BGL_LINKED_LIST_ARGS , typename OutputIter>
OutputIter remove_branch( typename BGL_LINKED_LIST::vertex_descriptor v, OutputIter it_out, BGL_LINKED_LIST & g) {
  return g.remove_branch(v,it_out);
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
template < BGL_LINKED_LIST_ARGS >
const typename BGL_LINKED_LIST::vertex_property_type& get(const BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::vertex_descriptor v) {
  return g[v];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param g The tree from which to draw the edge.
  * \param e The edge descriptor of the sought-after edge-property.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_LIST_ARGS >
const typename BGL_LINKED_LIST::edge_property_type& get( const BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::edge_descriptor e) {
  return g[e];
};

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by const-reference, to be associated to the given vertex.
  */
template < BGL_LINKED_LIST_ARGS >
void put( BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::vertex_descriptor v, const typename BGL_LINKED_LIST::vertex_property_type& value) {
  g[v] = value;
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by const-reference, to be associated to the given edge.
  */
template < BGL_LINKED_LIST_ARGS >
void put( BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::edge_descriptor e, const typename BGL_LINKED_LIST::edge_property_type& value) {
  g[e] = value;
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
  */
template < BGL_LINKED_LIST_ARGS >
void put( BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::vertex_descriptor v, typename BGL_LINKED_LIST::vertex_property_type&& value) {
  g[v] = std::move(value);
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
  */
template < BGL_LINKED_LIST_ARGS >
void put( BGL_LINKED_LIST & g, typename BGL_LINKED_LIST::edge_descriptor e, typename BGL_LINKED_LIST::edge_property_type&& value) {
  g[e] = std::move(value);
};

#endif


/**
  * Returns a reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by reference, associated to the given vertex descriptor.
  */
template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::vertex_property_type& get_property(typename BGL_LINKED_LIST::vertex_descriptor v, BGL_LINKED_LIST & g) {
  return g[v];
};

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_LINKED_LIST_ARGS >
const typename BGL_LINKED_LIST::vertex_property_type& get_property(typename BGL_LINKED_LIST::vertex_descriptor v, const BGL_LINKED_LIST & g) {
  return g[v];
};

/**
  * Returns a reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_LIST_ARGS >
typename BGL_LINKED_LIST::edge_property_type& get_property(typename BGL_LINKED_LIST::edge_descriptor e, BGL_LINKED_LIST & g) {
  return g[e];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_LINKED_LIST_ARGS >
const typename BGL_LINKED_LIST::edge_property_type& get_property(typename BGL_LINKED_LIST::edge_descriptor e, const BGL_LINKED_LIST & g) {
  return g[e];
};


#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template < BGL_LINKED_LIST_ARGS, typename T, typename Bundle>
struct property_map< BGL_LINKED_LIST, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_same< non_const_Bundle, typename BGL_LINKED_LIST::vertex_bundled > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, BGL_LINKED_LIST,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const BGL_LINKED_LIST,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


template < BGL_LINKED_LIST_ARGS, typename T, typename Bundle>
typename property_map< BGL_LINKED_LIST, T Bundle::* >::type
get( T Bundle::* p, BGL_LINKED_LIST& g) {
  return typename property_map< BGL_LINKED_LIST, T Bundle::* >::type(&g, p);
};

template < BGL_LINKED_LIST_ARGS, typename T, typename Bundle>
typename property_map< BGL_LINKED_LIST, T Bundle::* >::const_type
get( T Bundle::* p, const BGL_LINKED_LIST& g) {
  return typename property_map< BGL_LINKED_LIST, T Bundle::* >::const_type(&g, p);
};

#endif


/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

/* I'm not sure if the default version will work.

template < BGL_LINKED_LIST_ARGS, typename Property>
struct property_map<BGL_LINKED_LIST, Property>
  : property_map<
      typename BGL_LINKED_LIST::graph_type, 
      Property> {};

template < BGL_LINKED_LIST_ARGS, typename Property>
typename property_map< typename BGL_LINKED_LIST::graph_type, Property>::type
get(Property p, BGL_LINKED_LIST& g) {
  typedef typename BGL_LINKED_LIST::graph_type BaseGraph;
  return get(p, g.m_graph);
};

template < BGL_LINKED_LIST_ARGS, typename Property>
typename property_map< typename BGL_LINKED_LIST::graph_type, Property>::const_type
get(Property p, const BGL_LINKED_LIST& g) {
  return get(p, g.m_graph);
};

template < BGL_LINKED_LIST_ARGS, typename Property, typename Key>
typename property_map_value< typename BGL_LINKED_LIST::graph_type, Property>::type
get(Property p, const BGL_LINKED_LIST& g, const Key& k) {
  return get(p, g.m_graph, k);
};

template < BGL_LINKED_LIST_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_LINKED_LIST& g, const Key& k, const Value& val) {
  put(p, g.m_graph, k, val);
};
*/

#undef BGL_LINKED_LIST_ARGS
#undef BGL_LINKED_LIST



};


#endif


















