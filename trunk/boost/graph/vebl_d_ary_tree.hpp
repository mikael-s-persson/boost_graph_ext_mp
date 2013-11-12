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
    
    typedef VertexProperties vertex_bundled;
    typedef EdgeProperties edge_bundled;
    typedef void graph_bundled;
    
    struct value_type {
      std::size_t out_degree;
      vertex_property_type v;
      edge_property_type e[Arity];
      
      value_type() : out_degree(std::numeric_limits<std::size_t>::max()), v() { };
    };
    
    typedef std::vector< value_type > container_type;
    
    
    typedef std::size_t vertex_descriptor;
    
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
    
    
    typedef graph::detail::bfltree_edge_desc edge_descriptor;
    
    typedef graph::detail::vebltree_edge_validity<container_type, Arity> edge_validity;
    typedef filter_iterator< edge_validity, graph::detail::bfltree_oeiter > out_edge_iterator;
    typedef graph::detail::bfltree_ieiter in_edge_iterator;
    typedef filter_iterator< edge_validity, graph::detail::bfltree_eiter<Arity> > edge_iterator;
    
    typedef graph::detail::vebltree_vertex_validity<container_type, Arity> vertex_validity;
    typedef filter_iterator< vertex_validity, graph::detail::bfltree_viter > vertex_iterator;
    typedef vertex_iterator adjacency_iterator;
    typedef vertex_iterator child_vertex_iterator;
    
    
    
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
    
    value_type& get_vertex_value_impl(vertex_descriptor v_i) {
      return m_vertices[graph::detail::convert_bfl_to_vebl<Arity>(v_i, m_depth_recs)];
    };
    const value_type& get_vertex_value_impl(vertex_descriptor v_i) const {
      return m_vertices[graph::detail::convert_bfl_to_vebl<Arity>(v_i, m_depth_recs)];
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
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
    vertex_property_type& operator[]( const vertex_descriptor& v_i) {
      return get_vertex_value_impl(v_i).v;
    };
    /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
    const vertex_property_type& operator[]( const vertex_descriptor& v_i) const {
      return get_vertex_value_impl(v_i).v;
    };
    /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
    edge_property_type& operator[]( const edge_descriptor& e_i) {
      return get_vertex_value_impl(e_i.source_vertex).e[e_i.edge_index];
    };
    /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
    const edge_property_type& operator[]( const edge_descriptor& e_i) const {
      return get_vertex_value_impl(e_i.source_vertex).e[e_i.edge_index];
    };
    
    
    /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values.
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value for the newly created vertex.
     * \param ep The property value for the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
    std::pair< vertex_descriptor, edge_descriptor> add_child(const vertex_descriptor& v, 
                                                             const vertex_property_type& vp = vertex_property_type(), 
                                                             const edge_property_type& ep = edge_property_type()) {
      if( ( v >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(v)) )
        throw std::range_error("Cannot add child-node to an empty node!");
      std::size_t new_edge = 0;
      for(; new_edge < Arity; ++new_edge) {
        vertex_descriptor v_c = get_child_impl(v,new_edge);
        if( ( v_c >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(v_c)) )
          break;
      };
      if( new_edge == Arity ) 
        throw std::range_error("Cannot add child-node to a full node!");
      
      vertex_descriptor result = get_child_impl(v,new_edge);
      
      if( result >= m_vertices.size() )
        graph::detail::extend_vebl_storage<Arity>(m_vertices, m_depth_recs);
      
      value_type& v_prop      = get_vertex_value_impl(v);
      value_type& result_prop = get_vertex_value_impl(result);
      result_prop.out_degree = 0;
      result_prop.v = vp;
      v_prop.e[new_edge] = ep;
      ++(v_prop.out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
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
      if( ( v >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(v)) )
        throw std::range_error("Cannot add child-node to an empty node!");
      std::size_t new_edge = 0;
      for(; new_edge < Arity; ++new_edge) {
        vertex_descriptor v_c = get_child_impl(v,new_edge);
        if( ( v_c >= m_vertices.size() ) || !graph::detail::bfltree_is_vertex_valid(get_vertex_value_impl(v_c)) )
          break;
      };
      if( new_edge == Arity ) 
        throw std::range_error("Cannot add child-node to a full node!");
      
      vertex_descriptor result = get_child_impl(v,new_edge);
      
      if( result >= m_vertices.size() )
        graph::detail::extend_vebl_storage<Arity>(m_vertices, m_depth_recs);
      
      value_type& v_prop      = get_vertex_value_impl(v);
      value_type& result_prop = get_vertex_value_impl(result);
      result_prop.out_degree = 0;
      result_prop.v      = std::move(vp);
      v_prop.e[new_edge] = std::move(ep);
      ++(v_prop.out_degree);
      ++m_vertex_count;
      return std::make_pair(result, edge_descriptor(v, new_edge));
    };
#endif
    
    /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
    void remove_branch(const vertex_descriptor& v) {
      if( v >= m_vertices.size() )
        return;
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return;  // vertex is already deleted.
      --m_vertex_count;
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(std::size_t i = 0; i < Arity; ++i)
        remove_branch(get_child_impl(v,i));
      v_prop.out_degree = std::numeric_limits<std::size_t>::max();
      if( v != 0 )  // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
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
      if( v >= m_vertices.size() )
        return it_out;
      value_type& v_prop = get_vertex_value_impl(v);
      if( !graph::detail::bfltree_is_vertex_valid(v_prop) )
        return it_out;  // vertex is already deleted.
      --m_vertex_count; 
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
      *(it_out++) = std::move(v_prop.v);
#else
      *(it_out++) = v_prop.v;
#endif
      // this traversal order is intentional (traverse pre-order depth-first, and 
      // delay removal of empty tail elements as much as possible, such that it is only required once).
      for(std::size_t i = 0; i < Arity; ++i)
        it_out = remove_branch(get_child_impl(v,i),it_out);
      v_prop.out_degree = std::numeric_limits<std::size_t>::max();
      if( v != 0 )  // if the node is not the root one, then update the out-degree of the parent node:
        get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
      return it_out;
    };
    
    
    /**
     * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
     * \param vp The vertex-property to assign to the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(const vertex_property_type& vp = vertex_property_type()) {
      if(graph::detail::bfltree_is_vertex_valid(m_vertices[0]))
        remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].v = vp;
      ++m_vertex_count;
      return 0;
    };
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    /**
     * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
     * \param vp The vertex-property to move into the newly created root vertex.
     * \return The vertex-descriptor of the root of the tree.
     */
    vertex_descriptor create_root_vertex(vertex_property_type&& vp) {
      if(graph::detail::bfltree_is_vertex_valid(m_vertices[0]))
        remove_branch(0);
      m_vertices[0].out_degree = 0;
      m_vertices[0].v = std::move(vp);
      ++m_vertex_count;
      return 0;
    };
#endif
    
    
};



#if 0
/**
 * This is the tree-storage specifier for a D-ary COB-tree of a given Arity and CuttingDepth.
 */
template <std::size_t Arity = 2, std::size_t CuttingDepth = 8>
struct d_ary_cob_tree_storage { };


template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity, std::size_t CuttingDepth>
struct tree_storage<VertexDescriptor, EdgeDescriptor, d_ary_cob_tree_storage<Arity, CuttingDepth> > {
  typedef d_ary_cob_tree<VertexDescriptor, Arity, EdgeDescriptor, CuttingDepth> type;
};




template <std::size_t Arity, std::size_t CuttingDepth>
struct tree_storage_traits< d_ary_cob_tree_storage<Arity, CuttingDepth> > {
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
  typedef void* vertex_ptr;
  typedef typename d_ary_cob_tree<int,Arity,boost::no_property,CuttingDepth>::vertex_descriptor vertex_descriptor;  // the value-type doesn't affect the vertex_descriptor type (int is a dummy type here).
  typedef typename d_ary_cob_tree<int,Arity,boost::no_property,CuttingDepth>::edge_descriptor edge_descriptor;  // the value-type doesn't affect the edge_descriptor type (int is a dummy type here).
  typedef std::size_t edges_size_type;
  
};

#endif 



#define BGL_VEBL_D_ARY_TREE_ARGS std::size_t Arity, typename VertexProperties, typename EdgeProperties
#define BGL_VEBL_D_ARY_TREE vebl_d_ary_tree<Arity, VertexProperties, EdgeProperties>



template < BGL_VEBL_D_ARY_TREE_ARGS >
void swap(BGL_VEBL_D_ARY_TREE& lhs, BGL_VEBL_D_ARY_TREE& rhs) { lhs.swap(rhs); };




/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  source( const typename BGL_VEBL_D_ARY_TREE::edge_descriptor& e, const BGL_VEBL_D_ARY_TREE&) {
  return e.source_vertex;
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  target( const typename BGL_VEBL_D_ARY_TREE::edge_descriptor& e, const BGL_VEBL_D_ARY_TREE&) {
  return Arity * e.source_vertex + 1 + e.edge_index;
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::out_edge_iterator, typename BGL_VEBL_D_ARY_TREE::out_edge_iterator >
  out_edges( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::out_edge_iterator OutIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  
  graph::detail::bfltree_oeiter v_beg(graph::detail::bfltree_edge_desc(v, 0));
  graph::detail::bfltree_oeiter v_end(graph::detail::bfltree_edge_desc(v, Arity));
  return std::pair< OutIter, OutIter >(
    OutIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    OutIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t 
  out_degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return g.get_vertex_value_impl(v).out_degree;
};



/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/


template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::in_edge_iterator, typename BGL_VEBL_D_ARY_TREE::in_edge_iterator >
  in_edges( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::in_edge_iterator InIter;
  
  if(v == 0)
    return std::make_pair(InIter(graph::detail::bfltree_edge_desc(0,0)),
                          InIter(graph::detail::bfltree_edge_desc(0,0)));
  else {
    std::size_t u = (v - 1) / Arity;
    std::size_t e_id = (v - 1) % Arity;
    return std::make_pair(InIter(graph::detail::bfltree_edge_desc(u, e_id)),
                          InIter(graph::detail::bfltree_edge_desc(u, e_id+1)));
  };
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t
  in_degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  if(v == 0)
    return 0;
  else
    return 1;
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::size_t
  degree( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return in_degree(v, g) + out_degree(v, g);
};


/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

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

template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertices_size_type
  num_vertices( const BGL_VEBL_D_ARY_TREE& g) {
  return g.m_vertex_count;
};



/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::edge_iterator, typename BGL_VEBL_D_ARY_TREE::edge_iterator >
  edges( const BGL_VEBL_D_ARY_TREE& g) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_iterator EIter;
  typedef typename BGL_VEBL_D_ARY_TREE::container_type RawContainer;
  
  graph::detail::bfltree_eiter<Arity> v_beg(graph::detail::bfltree_edge_desc(0, 0));
  graph::detail::bfltree_eiter<Arity> v_end(graph::detail::bfltree_edge_desc(g.m_vertices.size(), Arity));
  return std::pair< EIter, EIter >(
    EIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_beg, v_end ), 
    EIter( graph::detail::vebltree_edge_validity<RawContainer,Arity>(&g.m_vertices, &g.m_depth_recs), v_end, v_end ));
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertices_size_type
  num_edges( const BGL_VEBL_D_ARY_TREE& g) {
  return num_vertices(g) - 1;
};



/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::edge_descriptor, bool >
  edge( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor&,
        const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v,
        const BGL_VEBL_D_ARY_TREE&) {
  typedef typename BGL_VEBL_D_ARY_TREE::edge_descriptor Edge;
  return std::make_pair(Edge((v - 1) / Arity, (v - 1) % Arity),true);
};


/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/


template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  get_root_vertex( const BGL_VEBL_D_ARY_TREE& g) {
  return 0;
};

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

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::adjacency_iterator, typename BGL_VEBL_D_ARY_TREE::adjacency_iterator >
  adjacent_vertices( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BGL_VEBL_D_ARY_TREE& g) {
  return child_vertices(v, g);
};


/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/


template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( BGL_VEBL_D_ARY_TREE& g) {
  return g.create_root_vertex();
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_descriptor, typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v);
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
void remove_branch( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, BGL_VEBL_D_ARY_TREE& g) {
  return g.remove_branch(v);
};



/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/


template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& vp, BGL_VEBL_D_ARY_TREE& g) {
  return g.create_root_vertex(vp);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template < BGL_VEBL_D_ARY_TREE_ARGS >
typename BGL_VEBL_D_ARY_TREE::vertex_descriptor
  create_root( typename BGL_VEBL_D_ARY_TREE::vertex_property_type&& vp, BGL_VEBL_D_ARY_TREE& g) {
  return g.create_root_vertex(std::move(vp));
};
#endif

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_descriptor, typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v, 
                    const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& vp,
                    BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v,vp);
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_descriptor, typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v,
                    const typename BGL_VEBL_D_ARY_TREE::vertex_property_type& vp,
                    const typename BGL_VEBL_D_ARY_TREE::edge_property_type& ep,
                    BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v,vp,ep);
};

#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_descriptor, typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v,
                    typename BGL_VEBL_D_ARY_TREE::vertex_property_type&& vp,
                    BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v,std::move(vp));
};

template < BGL_VEBL_D_ARY_TREE_ARGS >
std::pair< typename BGL_VEBL_D_ARY_TREE::vertex_descriptor, typename BGL_VEBL_D_ARY_TREE::edge_descriptor >
  add_child_vertex( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v,
                    typename BGL_VEBL_D_ARY_TREE::vertex_property_type&& vp,
                    typename BGL_VEBL_D_ARY_TREE::edge_property_type&& ep,
                    BGL_VEBL_D_ARY_TREE& g) {
  return g.add_child(v,std::move(vp),std::move(ep));
};
#endif

template < BGL_VEBL_D_ARY_TREE_ARGS, typename OutputIter>
OutputIter remove_branch( const typename BGL_VEBL_D_ARY_TREE::vertex_descriptor& v,
                          OutputIter it_out,
                          BGL_VEBL_D_ARY_TREE& g) {
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


















