// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file vebl_d_ary_tree.hpp
 * 
 * This library provides a class that implements a D-Ary von Emde Boas Layout tree that stores
 * elements in a contiguous array in a cache-friendly recursive layout. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date September 2013
 */


#ifndef BOOST_VEBL_D_ARY_TREE_HPP
#define BOOST_VEBL_D_ARY_TREE_HPP

#include <boost/graph/properties.hpp>
#include <boost/config.hpp>

#include <vector>
#include <stdexcept>
#include <utility>
#include <iterator>

#include <boost/graph/more_property_maps.hpp>
#include <boost/graph/tree_concepts.hpp>
#include <boost/graph/tree_traits.hpp>
#include <boost/graph/detail/vebl_tree_iterators.hpp>

namespace boost {



struct vebl_d_ary_tree_tag { };
 


/**
 * This class implements a D-Ary von Emde Boas Layout tree that stores
 * elements in a contiguous array in a cache-friendly recursive layout. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. 
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <std::size_t Arity = 2,
          typename VertexProperties = boost::no_property,
          typename EdgeProperties = boost::no_property >
class vebl_d_ary_tree
{
  public:
    typedef vebl_d_ary_tree<Arity, VertexProperties, EdgeProperties> self;
    
    typedef VertexProperties vertex_property_type;
    typedef EdgeProperties edge_property_type;
    
    typedef vertex_property_type vertex_bundled;
    typedef edge_property_type edge_bundled;
    typedef void graph_bundled;
    
    typedef graph::detail::bfltree_value_type<vertex_property_type, edge_property_type> value_type;
    
    typedef std::vector< value_type > container_type;
    
    typedef std::size_t vertex_descriptor;
    typedef graph::detail::bfltree_edge_desc edge_descriptor;
    
    typedef std::size_t vertices_size_type;
    typedef vertices_size_type edges_size_type;
    typedef vertices_size_type degree_size_type;
    
    /**
     * This static member function outputs the null-vertex (invalid vertex descriptor).
     * \return A null-vertex descriptor (invalid vertex descriptor).
     */
    static vertex_descriptor null_vertex() { 
      return std::numeric_limits<std::size_t>::max();
    };
    
    typedef graph::detail::vebltree_edge_validity<container_type, Arity> edge_validity;
    typedef filter_iterator< edge_validity, graph::detail::bfltree_eiter > edge_iterator;
    typedef edge_iterator out_edge_iterator;
    typedef graph::detail::bfltree_eiter in_edge_iterator;
    
    typedef graph::detail::vebltree_vertex_validity<container_type, Arity> vertex_validity;
    typedef filter_iterator< vertex_validity, graph::detail::bfltree_viter > vertex_iterator;
    typedef vertex_iterator adjacency_iterator;
    typedef vertex_iterator child_vertex_iterator;
    typedef graph::detail::bfltree_viter inv_adjacency_iterator;
    
    typedef boost::directed_tag directed_category;
    typedef boost::disallow_parallel_edge_tag edge_parallel_category;
    
    struct traversal_category : 
      virtual public boost::incidence_graph_tag,
      virtual public boost::adjacency_graph_tag,
      virtual public boost::bidirectional_graph_tag,
      virtual public boost::vertex_list_graph_tag,
      virtual public boost::edge_list_graph_tag { };
    
    typedef vebl_d_ary_tree_tag graph_tag;
    
    
  public: //private:
    
    container_type m_vertices;
    vertices_size_type m_vertex_count;
    
    graph::detail::vebl_depth_records m_depth_recs;
    
    value_type& get_vertex_value_impl(vertex_descriptor v) {
      return m_vertices[graph::detail::convert_bfl_to_vebl<Arity>(v, m_depth_recs)];
    };
    const value_type& get_vertex_value_impl(vertex_descriptor v) const {
      return m_vertices[graph::detail::convert_bfl_to_vebl<Arity>(v, m_depth_recs)];
    };
    
    static vertex_descriptor get_child_impl(vertex_descriptor v, std::size_t edge_id) {
      return Arity * v + 1 + edge_id;
    };
    
    
  public:
    
    
    /**
     * Construct the D-ary BF-tree with a given reserved depth.
     * \param aDepth The depth of the graph to reserve space for.
     */
    vebl_d_ary_tree(vertices_size_type aDepth = 0) : m_vertex_count(0) {
      m_depth_recs.T.push_back(0);
      m_depth_recs.B.push_back(1);
      m_depth_recs.D.push_back(0);
      for(vertices_size_type i = 0; i < aDepth; ++i)
        graph::detail::extend_vebl_depth_records<Arity>(m_depth_recs);
      m_vertices.resize(graph::detail::s_treesize<Arity>(aDepth));
    };
    
    /**
     * Checks if the tree is empty.
     * \return True if the tree is empty.
     */
    bool empty() const { return m_vertex_count == 0; };
    
    /**
     * Returns the size of the tree (the number of vertices it contains).
     * \return The size of the tree (the number of vertices it contains).
     */
    std::size_t size() const { return m_vertex_count; };
    
    /**
     * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
     * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
     */
    std::size_t capacity() const { return m_vertices.capacity(); };
    
    /**
     * Returns the depth of the tree.
     * \return The depth of the tree.
     */
    std::size_t depth() const { return m_depth_recs.T.size(); };
    
    /**
     * Standard swap function.
     */
    void swap(self& rhs) {
      using std::swap;
      m_vertices.swap(rhs.m_vertices);
      swap(m_vertex_count, rhs.m_vertex_count);
      m_depth_recs.T.swap(rhs.m_depth_recs.T);
      m_depth_recs.B.swap(rhs.m_depth_recs.B);
      m_depth_recs.D.swap(rhs.m_depth_recs.D);
    };
    
    /**
     * Clears the tree of all vertices and edges.
     */
    void clear() { 
      m_vertices.clear();
      m_vertices.resize(1);
      m_depth_recs.T.clear();
      m_depth_recs.T.push_back(0);
      m_depth_recs.B.clear();
      m_depth_recs.B.push_back(1);
      m_depth_recs.D.clear();
      m_depth_recs.D.push_back(0);
      m_vertex_count = 0;
    };
    
    /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[]( const vertex_descriptor& v) {
      return get_vertex_value_impl(v).vertex();
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[]( const vertex_descriptor& v) const {
      return get_vertex_value_impl(v).vertex();
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[]( const edge_descriptor& e) {
      return get_vertex_value_impl(e.target_vertex).edge();
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[]( const edge_descriptor& e) const {
      return get_vertex_value_impl(e.target_vertex).edge();
    };
    
    
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values.
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value for the newly created vertex.
     * \param ep The property value for the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
    std::pair< vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, const vertex_property_type& vp, const edge_property_type& ep) {
#else
    template <typename VP, typename EP>
    std::pair< vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v, VP&& vp, EP&& ep) {
#endif
      if( ( v >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(v)) )
        return std::pair<vertex_descriptor, edge_descriptor>(null_vertex(), edge_descriptor(null_vertex()));
      vertex_descriptor result = get_child_impl(v, 0);
      for(; result < get_child_impl(v, Arity); ++result)
        if( ( result >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(result)) )
          break;
      if( result == get_child_impl(v, Arity) ) 
        return std::pair<vertex_descriptor, edge_descriptor>(null_vertex(), edge_descriptor(null_vertex()));
      
      if( result >= m_vertices.size() )
        graph::detail::extend_vebl_storage<Arity>(m_vertices, m_depth_recs);
      
      value_type& result_prop = get_vertex_value_impl(result);
      result_prop.out_degree = 0;
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
      result_prop.vertex() = vp;
      result_prop.edge()   = ep;
#else
      result_prop.vertex() = std::forward<VP>(vp);
      result_prop.edge()   = std::forward<EP>(ep);
#endif
      ++(get_vertex_value_impl(v).out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(result));
    };
    
    
    template <typename VertexOIter, typename EdgeOIter>
    void clear_children_impl(vertex_descriptor v, VertexOIter& vit_out, EdgeOIter& eit_out) {
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(std::size_t i = 0; i < Arity; ++i) {
        vertex_descriptor next_v = get_child_impl(v, i);
        if(next_v >= m_vertices.size())
          break;
        value_type& next_v_prop = get_vertex_value_impl(next_v);
        if( !graph::detail::bfltree_is_vertex_valid(next_v_prop) )
          continue;
        --m_vertex_count;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(vit_out++) = std::move(next_v_prop.vertex());
        *(eit_out++) = std::move(next_v_prop.edge());
#else
        *(vit_out++) = next_v_prop.vertex();
        *(eit_out++) = next_v_prop.edge();
#endif
        clear_children_impl(next_v, vit_out, eit_out);
      };
      get_vertex_value_impl(v).out_degree = std::numeric_limits<std::size_t>::max();
      if( v != 0 )  // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
    };
    
    /**
     * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
     * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
     * \param v The root of the sub-tree to be removed.
     * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
     * \return The output-iterator after the collection of all the removed vertices.
     */
    template <typename VertexOIter, typename EdgeOIter>
    std::pair<VertexOIter, EdgeOIter> clear_children(vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out) {
      if(v >= m_vertices.size())
        return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out); // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);  // vertex is already deleted.
      clear_children_impl(v, vit_out, eit_out);
      return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
    };
    
    /**
     * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
     * recording the vertex-properties of all the removed vertices into an output-iterator.
     * \param v The root of the sub-tree to be removed.
     * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \return The output-iterator after the collection of all the removed vertices.
     */
    template <typename VertexOIter>
    VertexOIter clear_children(vertex_descriptor v, VertexOIter vit_out) {
      if(v >= m_vertices.size())
        return vit_out; // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return vit_out;  // vertex is already deleted.
      graph::detail::ignore_output_iter eit_out;
      clear_children_impl(v, vit_out, eit_out);
      return vit_out;
    };
    
    /**
     * Removes a branch (sub-tree) starting from but excluding the given vertex.
     * \param v The root of the sub-tree to be removed.
     */
    void clear_children(vertex_descriptor v) {
      if(v >= m_vertices.size())
        return; // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return;  // vertex is already deleted.
      graph::detail::ignore_output_iter vit_out, eit_out;
      clear_children_impl(v, vit_out, eit_out);
    };
    
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex, while 
     * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
     * \return The output-iterator after the collection of all the removed vertices.
     * \note The first vertex-property to figure in the output range is that of the vertex v.
     */
    template <typename VertexOIter, typename EdgeOIter>
    std::pair<VertexOIter, EdgeOIter> remove_branch(vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out) {
      if(v >= m_vertices.size())
        return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out); // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);  // vertex is already deleted.
      if( v == 0 ) {
        *(eit_out++) = edge_property_type();
      } else {
        // in-edge:  u = (v - 1) / Arity;  e_id = (v - 1) % Arity;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
        *(eit_out++) = std::move(v_prop.edge());
#else
        *(eit_out++) = v_prop.edge();
#endif
        // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
      };
      --m_vertex_count;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      *(vit_out++) = std::move(v_prop.vertex());
#else
      *(vit_out++) = v_prop.vertex();
#endif
      clear_children_impl(v, vit_out, eit_out);
      return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
    };
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex, while 
     * recording the vertex-properties of all the removed vertices into an output-iterator.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
     * \return The output-iterator after the collection of all the removed vertices.
     * \note The first vertex-property to figure in the output range is that of the vertex v.
     */
    template <typename VertexOIter>
    VertexOIter remove_branch(vertex_descriptor v, VertexOIter vit_out) {
      if(v >= m_vertices.size())
        return vit_out; // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return vit_out;  // vertex is already deleted.
      if( v != 0 ) {
        // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
      };
      --m_vertex_count;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      *(vit_out++) = std::move(v_prop.vertex());
#else
      *(vit_out++) = v_prop.vertex();
#endif
      graph::detail::ignore_output_iter eit_out;
      clear_children_impl(v, vit_out, eit_out);
      return vit_out;
    };
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
    void remove_branch(vertex_descriptor v) {
      if(v >= m_vertices.size())
        return; // vertex is already deleted.
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return;  // vertex is already deleted.
      if( v != 0 ) {
        // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
      };
      --m_vertex_count;
      graph::detail::ignore_output_iter vit_out, eit_out;
      clear_children_impl(v, vit_out, eit_out);
    };
    
    
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
     * \param vp The vertex-property to assign to the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(graph::detail::bfltree_is_vertex_valid(m_vertices[0]))
        remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].vertex() = vp;
      ++m_vertex_count;
      return 0;
    };
#else
    /**
     * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
     * \param vp The vertex-property to move into the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    template <typename VProp>
    vertex_descriptor create_root_vertex(VProp&& vp) {
      if(graph::detail::bfltree_is_vertex_valid(m_vertices[0]))
        remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].vertex() = std::forward<VProp>(vp);
      ++m_vertex_count;
      return 0;
    };
#endif
    
    
};





/**
 * This is the tree-storage specifier for a D-Ary von Emde Boas Layout tree of a given Arity.
 */
template <std::size_t Arity = 2>
struct vebl_d_ary_tree_storage { };


template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity>
struct tree_storage<VertexDescriptor, EdgeDescriptor, vebl_d_ary_tree_storage<Arity> > {
  typedef vebl_d_ary_tree<Arity, VertexDescriptor, EdgeDescriptor> type;
};


template <std::size_t Arity>
struct tree_storage_traits< vebl_d_ary_tree_storage<Arity> > {
  typedef boost::mpl::true_ is_rand_access;
  typedef boost::mpl::true_ is_bidir;
  typedef boost::mpl::true_ is_directed;
  
  typedef typename boost::mpl::if_< is_bidir,
    boost::bidirectional_tag,
    typename boost::mpl::if_< is_directed,
      boost::directed_tag, boost::undirected_tag
    >::type
  >::type directed_category;
  
  typedef boost::disallow_parallel_edge_tag edge_parallel_category;
  
  typedef std::size_t vertices_size_type;
  typedef std::size_t vertex_descriptor;
  typedef std::size_t edges_size_type;
  typedef graph::detail::bfltree_edge_desc edge_descriptor;
  
};




#define BGL_VEBL_D_ARY_TREE_ARGS std::size_t Arity, typename VertexProperties, typename EdgeProperties
#define BGL_VEBL_D_ARY_TREE vebl_d_ary_tree<Arity, VertexProperties, EdgeProperties>



/**
 * Standard swap function. Swaps the contents of two objects.
 * \param lhs The left-hand-side of the swap.
 * \param rhs The right-hand-side of the swap.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void swap(BGL_VEBL_D_ARY_TREE& lhs, BGL_VEBL_D_ARY_TREE& rhs) { lhs.swap(rhs); };



/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/


/**
 * Returns the source vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The source vertex of the given edge descriptor.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  source( typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, const BGL_VEBL_D_ARY_TREE&) {
  return (e.target_vertex - 1) / Arity;
};

/**
 * Returns the target vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The target vertex of the given edge descriptor.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  target( typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, const BGL_VEBL_D_ARY_TREE&) {
  return e.target_vertex;
};

/**
 * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the out-edges of a given vertex descriptor.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::out_edge_iterator, typename BGL_VEBL_D_ARY_TREE::out_edge_iterator >
  out_edges( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::out_edge_iterator OutIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  // Arity * v + 1 + edge_index;
  graph::detail::bfltree_eiter v_beg(graph::detail::bfltree_edge_desc(Arity * v + 1));
  graph::detail::bfltree_eiter v_end(graph::detail::bfltree_edge_desc(Arity * (v + 1) + 1));
  return std::pair< OutIter, OutIter >(
    OutIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    OutIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};

/**
 * Returns the out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The out-degree of the given vertex descriptor.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t 
  out_degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return g.get_vertex_value_impl(v).out_degree;
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
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::in_edge_iterator, typename BGL_VEBL_D_ARY_TREE::in_edge_iterator >
  in_edges( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::in_edge_iterator InIter;
  if(v == 0)
    return std::make_pair(InIter(graph::detail::bfltree_edge_desc(0)), InIter(graph::detail::bfltree_edge_desc(0)));
  else
    return std::make_pair(InIter(graph::detail::bfltree_edge_desc(v)), InIter(graph::detail::bfltree_edge_desc(v+1)));
};

/**
 * Returns the in-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t in_degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  if(v == 0)
    return 0;
  else
    return 1;
};

/**
 * Returns the in-degree plus out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree plus out-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return in_degree(v, g) + out_degree(v, g);
};


/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the vertices of the tree.
 * \param g The graph.
 * \return The vertex iterator range for all the vertices of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_iterator, typename BGL_VEBL_D_ARY_TREE::vertex_iterator >
  vertices( const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::vertex_iterator VIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  graph::detail::bfltree_viter v_beg(0);
  graph::detail::bfltree_viter v_end(g.m_vertices.size());
  return std::pair< VIter, VIter >(
    VIter( graph::detail::vebltree_vertex_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    VIter( graph::detail::vebltree_vertex_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};

/**
 * Returns the size of the tree (the number of vertices it contains).
 * \param g The graph.
 * \return The size of the tree (the number of vertices it contains).
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertices_size_type
  num_vertices( const BGL_VEBL_D_ARY_TREE& g) {
  return g.m_vertex_count;
};



/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/


/**
 * Returns the edge iterator range for all the edges of the tree.
 * \param g The graph.
 * \return The edge iterator range for all the edges of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair<
 typename BGL_VEBL_D_ARY_TREE::edge_iterator,
 typename BGL_VEBL_D_ARY_TREE::edge_iterator >
  edges( const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_iterator EIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  
  graph::detail::bfltree_eiter v_beg(graph::detail::bfltree_edge_desc(1));
  graph::detail::bfltree_eiter v_end(graph::detail::bfltree_edge_desc(g.m_vertices.size()));
  return std::pair< EIter, EIter >(
    EIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    EIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};

/**
 * Returns the number of edges in the tree.
 * \param g The graph.
 * \return The number of edges in the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertices_size_type
  num_edges( const BGL_VEBL_D_ARY_TREE& g) {
  return num_vertices(g) - 1;
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
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::edge_descriptor, bool >
  edge( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor u,
        typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
        const BGL_VEBL_D_ARY_TREE&) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_descriptor Edge;
  if( u == (v - 1) / Arity )
    return std::make_pair(Edge(v), true);
  else
    return std::make_pair(Edge(0), false);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  get_root_vertex( const BGL_VEBL_D_ARY_TREE& g) {
  return 0;
};

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_iterator, typename BGL_VEBL_D_ARY_TREE::vertex_iterator >
  child_vertices( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::vertex_iterator VIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  graph::detail::bfltree_viter v_beg(Arity * v + 1);
  graph::detail::bfltree_viter v_end(Arity * (v + 1) + 1);
  return std::pair< VIter, VIter >(
    VIter( graph::detail::vebltree_vertex_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    VIter( graph::detail::vebltree_vertex_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};



/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::adjacency_iterator, typename BGL_VEBL_D_ARY_TREE::adjacency_iterator >
  adjacent_vertices( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return child_vertices(v, g);
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
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  parent_vertex(typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BGL_VEBL_D_ARY_TREE&) {
  if( v == 0 )
    return BGL_VEBL_D_ARY_TREE::null_vertex();
  else
    return (v - 1) / Arity;
};


/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for the parent-vertex of a given vertex of the tree.
 * \param v The vertex descriptor whose parent is sought.
 * \param g The graph.
 * \return The vertex iterator range for the parent-vertex of a given vertex of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair<
 typename BGL_VEBL_D_ARY_TREE::inv_adjacency_iterator,
 typename BGL_VEBL_D_ARY_TREE::inv_adjacency_iterator >
  inv_adjacent_vertices( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::inv_adjacency_iterator InvAdjIter;
  if( v == 0 )
    return std::pair<InvAdjIter, InvAdjIter>(graph::detail::bfltree_viter(0), graph::detail::bfltree_viter(0));
  else
    return std::pair<InvAdjIter, InvAdjIter>(graph::detail::bfltree_viter((v - 1) / Arity), graph::detail::bfltree_viter((v - 1) / Arity + 1));
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void remove_branch( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, BGL_VEBL_D_ARY_TREE& g) {
  return g.remove_branch(v);
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void clear_children( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, BGL_VEBL_D_ARY_TREE& g) {
  return g.clear_children(v);
};

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( BGL_VEBL_D_ARY_TREE& g) {
  if(graph::detail::bfltree_is_vertex_valid(g.m_vertices[0]))
    remove_branch(0, g);
  g.m_vertices[0].out_degree = 0;
  ++g.m_vertex_count;
  return 0;
};

/**
 * Adds a child vertex to the given parent vertex, and default-initializes the properties of 
 * the newly created vertex and edge.
 * \param v The parent vertex to which a child will be added.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< 
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor,
typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::vertex_property_type VProp;
  typedef typename BGL_VEBL_D_ARY_TREE::edge_property_type EProp;
  return g.add_child(v, VProp(), EProp());
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( const VProp& vp, BGL_VEBL_D_ARY_TREE& g) {
  if(graph::detail::bfltree_is_vertex_valid(g.m_vertices[0]))
    remove_branch(0, g);
  g.m_vertices[0].out_degree = 0;
  g.m_vertices[0].vertex() = vp;
  ++g.m_vertex_count;
  return 0;
};

#else

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( VProp&& vp, BGL_VEBL_D_ARY_TREE& g) {
  if(graph::detail::bfltree_is_vertex_valid(g.m_vertices[0]))
    remove_branch(0, g);
  g.m_vertices[0].out_degree = 0;
  g.m_vertices[0].vertex() = std::forward<VProp>(vp);
  ++g.m_vertex_count;
  return 0;
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp >
std::pair< 
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor,
typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                    const VProp& vp, BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_property_type EProp;
  return g.add_child(v, vp, EProp());
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp, typename EProp >
std::pair< 
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor,
typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                    const VProp& vp, const EProp& ep, BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v, vp, ep);
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp >
std::pair< 
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor,
typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                    VProp&& vp,
                    BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_property_type EProp;
  return g.add_child(v, std::forward<VProp>(vp), EProp());
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VProp, typename EProp >
std::pair< 
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor,
typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                    VProp&& vp, EProp&& ep,
                    BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v, std::forward<VProp>(vp), std::forward<EProp>(ep));
};

#endif


/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS, typename OutputIter >
OutputIter remove_branch( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                          OutputIter it_out, BGL_VEBL_D_ARY_TREE& g) {
  return g.remove_branch(v, it_out);
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter >
std::pair<VertexOIter, EdgeOIter> remove_branch( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                 VertexOIter vit_out, EdgeOIter eit_out, BGL_VEBL_D_ARY_TREE& g) {
  return g.remove_branch(v, vit_out, eit_out);
};


/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while 
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template < BGL_VEBL_D_ARY_TREE_ARGS, typename OutputIter >
OutputIter clear_children( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                           OutputIter it_out, BGL_VEBL_D_ARY_TREE& g) {
  return g.clear_children(v, it_out);
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
template < BGL_VEBL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter >
std::pair<VertexOIter, EdgeOIter> clear_children( typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                  VertexOIter vit_out, EdgeOIter eit_out, BGL_VEBL_D_ARY_TREE& g) {
  return g.clear_children(v, vit_out, eit_out);
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
template < BGL_VEBL_D_ARY_TREE_ARGS >
const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& 
  get(const BGL_VEBL_D_ARY_TREE & g, typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v) {
  return g[v];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param g The tree from which to draw the edge.
  * \param e The edge descriptor of the sought-after edge-property.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
const typename BGL_VEBL_D_ARY_TREE::edge_property_type& 
  get( const BGL_VEBL_D_ARY_TREE & g, typename BGL_VEBL_D_ARY_TREE::edge_descriptor e) {
  return g[e];
};

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by const-reference, to be associated to the given vertex.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void put( BGL_VEBL_D_ARY_TREE & g, 
          typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, 
          const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& value) {
  g[v] = value;
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by const-reference, to be associated to the given edge.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void put( BGL_VEBL_D_ARY_TREE & g, 
          typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, 
          const typename BGL_VEBL_D_ARY_TREE::edge_property_type& value) {
  g[e] = value;
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void put( BGL_VEBL_D_ARY_TREE & g, 
          typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, 
          typename BGL_VEBL_D_ARY_TREE::vertex_property_type&& value) {
  g[v] = std::move(value);
};

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
void put( BGL_VEBL_D_ARY_TREE & g, 
          typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, 
          typename BGL_VEBL_D_ARY_TREE::edge_property_type&& value) {
  g[e] = std::move(value);
};

#endif


/**
  * Returns a reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by reference, associated to the given vertex descriptor.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_property_type& 
  get_property(typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, BGL_VEBL_D_ARY_TREE & g) {
  return g[v];
};

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& 
  get_property(typename BGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BGL_VEBL_D_ARY_TREE & g) {
  return g[v];
};

/**
  * Returns a reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by reference, associated to the given edge descriptor.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::edge_property_type& 
  get_property(typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, BGL_VEBL_D_ARY_TREE & g) {
  return g[e];
};

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template < BGL_VEBL_D_ARY_TREE_ARGS >
const typename BGL_VEBL_D_ARY_TREE::edge_property_type& 
  get_property(typename BGL_VEBL_D_ARY_TREE::edge_descriptor e, const BGL_VEBL_D_ARY_TREE & g) {
  return g[e];
};




#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES

template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
struct property_map< BGL_VEBL_D_ARY_TREE, T Bundle::* > {
  typedef typename remove_const< Bundle >::type non_const_Bundle;
  typedef typename remove_const< T >::type non_const_T;
  typedef is_convertible< typename BGL_VEBL_D_ARY_TREE::vertex_bundled*, non_const_Bundle* > is_vertex_bundle;
  typedef bundle_member_property_map< non_const_T, BGL_VEBL_D_ARY_TREE,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > type;
  typedef bundle_member_property_map< const non_const_T, const BGL_VEBL_D_ARY_TREE,
    typename mpl::if_< is_vertex_bundle, vertex_bundle_t, edge_bundle_t >::type > const_type;
};


template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
typename property_map< BGL_VEBL_D_ARY_TREE, T Bundle::* >::type
get( T Bundle::* p, BGL_VEBL_D_ARY_TREE& g) {
  return typename property_map< BGL_VEBL_D_ARY_TREE, T Bundle::* >::type(&g, p);
};

template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
typename property_map< BGL_VEBL_D_ARY_TREE, T Bundle::* >::const_type
get( T Bundle::* p, const BGL_VEBL_D_ARY_TREE& g) {
  return typename property_map< BGL_VEBL_D_ARY_TREE, T Bundle::* >::const_type(&g, p);
};



template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
const typename remove_const< T >::type& get( T Bundle::* p, const BGL_VEBL_D_ARY_TREE& g, const Key& k) {
  return (g[k]).*p;
};

template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
void put( T Bundle::* p, BGL_VEBL_D_ARY_TREE& g, const Key& k, const T& val) {
  (g[k]).*p = val;
};


#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
void put( T Bundle::* p, BGL_VEBL_D_ARY_TREE& g, const Key& k, T&& val) {
  (g[k]).*p = std::move(val);
};

#endif



#endif


/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct vebl_d_ary_tree_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    typedef typename property_value<Property, Tag>::type value_type;
    
    typedef tagged_from_bundle_property_map<value_type, Graph, Tag> type;
    typedef tagged_from_bundle_property_map<const value_type, const Graph, Tag> const_type;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<vebl_d_ary_tree_tag> {
  typedef vebl_d_ary_tree_property_selector type;
};

template <>
struct edge_property_selector<vebl_d_ary_tree_tag> {
  typedef vebl_d_ary_tree_property_selector type;
};


template < BGL_VEBL_D_ARY_TREE_ARGS, typename Property>
typename property_map< BGL_VEBL_D_ARY_TREE, Property>::type
get(Property p, BGL_VEBL_D_ARY_TREE& g) {
  typedef typename property_map< BGL_VEBL_D_ARY_TREE, Property>::type Map;
  return Map(&g, p);
};

template < BGL_VEBL_D_ARY_TREE_ARGS, typename Property>
typename property_map< BGL_VEBL_D_ARY_TREE, Property>::const_type
get(Property p, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename property_map< BGL_VEBL_D_ARY_TREE, Property>::const_type Map;
  return Map(&g, p);
};


template < BGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key>
typename property_map_value< BGL_VEBL_D_ARY_TREE, Property>::type
get(Property p, const BGL_VEBL_D_ARY_TREE& g, const Key& k) {
  return get_property_value(g[k], p);
};

template < BGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_VEBL_D_ARY_TREE& g, const Key& k, const Value& val) {
  get_property_value(g[k], p) = val;
};


#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES

template < BGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_VEBL_D_ARY_TREE& g, const Key& k, Value&& val) {
  get_property_value(g[k], p) = std::move(val);
};

#endif






#undef BGL_VEBL_D_ARY_TREE_ARGS
#undef BGL_VEBL_D_ARY_TREE



};



#endif


















