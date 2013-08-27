
/*
 *    Copyright 2013 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/tree_adaptor.hpp>
#include <boost/graph/pooled_adjacency_list.hpp>
#include <boost/graph/adjacency_list_BC.hpp>


#if 0
#define TEST_PRINT_REACHED_MARKER std::cout << __LINE__ << " reached!" << std::endl;
#else
#define TEST_PRINT_REACHED_MARKER 
#endif



#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE bgl_graphs
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/front_inserter.hpp>

using namespace boost;


typedef mpl::list< 
// These must be disabled because the old adj-list implementation does not pass the move-semantics tests.
  // adjacency_list< vecS,  vecS,  bidirectionalS, int, int>,
  // adjacency_list< listS, vecS,  bidirectionalS, int, int>,
  // adjacency_list< vecS,  listS, bidirectionalS, int, int>,
  // adjacency_list< listS, listS, bidirectionalS, int, int>,
  // adjacency_list< vecS,  vecS,  directedS, int, int>,
  // adjacency_list< listS, vecS,  directedS, int, int>,
  // adjacency_list< vecS,  listS, directedS, int, int>,
  // adjacency_list< listS, listS, directedS, int, int>,
  // pooled_adjacency_list<bidirectionalS, int, int>,
  // pooled_adjacency_list<directedS, int, int>
  > intint_adjlist_types;
  
  
  
typedef mpl::list< 
  adjacency_list_BC< vecBC,  vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  listBC, bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< listBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< listBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< listBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  vecBC,  directedS, int, int>,
  adjacency_list_BC< vecBC,  listBC, directedS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, directedS, int, int>,
  adjacency_list_BC< listBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< listBC, listBC, directedS, int, int>,
  adjacency_list_BC< listBC, poolBC, directedS, int, int>,
  adjacency_list_BC< poolBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< poolBC, listBC, directedS, int, int>,
  adjacency_list_BC< poolBC, poolBC, directedS, int, int> > intint_adjlistBC_nosets_types;

typedef mpl::list< 
  adjacency_list_BC< vecBC,  vecBC,  undirectedS, int, int>,
  adjacency_list_BC< vecBC,  listBC, undirectedS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, undirectedS, int, int>,
  adjacency_list_BC< listBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< listBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< listBC, poolBC, undirectedS, int, int>,
  adjacency_list_BC< poolBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< poolBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< poolBC, poolBC, undirectedS, int, int> > intint_adjlistBC_nosets_undir_types;
  
  
  
typedef mpl::list< 
  adjacency_list_BC< vecBC,  vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  listBC, bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< vecBC,  vecBC,  directedS, int, int>,
  adjacency_list_BC< vecBC,  listBC, directedS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, directedS, int, int>,
  adjacency_list_BC< vecBC,  vecBC,  undirectedS, int, int>,
  adjacency_list_BC< vecBC,  listBC, undirectedS, int, int>,
  adjacency_list_BC< vecBC,  poolBC, undirectedS, int, int> > intint_adjlistBC_vec_types;
  
typedef mpl::list< 
  adjacency_list_BC< listBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< listBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< listBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< listBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< listBC, listBC, directedS, int, int>,
  adjacency_list_BC< listBC, poolBC, directedS, int, int>,
  adjacency_list_BC< listBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< listBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< listBC, poolBC, undirectedS, int, int> > intint_adjlistBC_list_types;
  
typedef mpl::list< 
  adjacency_list_BC< poolBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< poolBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< poolBC, listBC, directedS, int, int>,
  adjacency_list_BC< poolBC, poolBC, directedS, int, int>,
  adjacency_list_BC< poolBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< poolBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< poolBC, poolBC, undirectedS, int, int> > intint_adjlistBC_pool_types;
  
typedef mpl::list< 
  adjacency_list_BC< setBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< setBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< setBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< setBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< setBC, listBC, directedS, int, int>,
  adjacency_list_BC< setBC, poolBC, directedS, int, int>,
  adjacency_list_BC< setBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< setBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< setBC, poolBC, undirectedS, int, int> > intint_adjlistBC_set_types;
  
typedef mpl::list< 
  adjacency_list_BC< unordered_setBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< unordered_setBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< unordered_setBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< unordered_setBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< unordered_setBC, listBC, directedS, int, int>,
  adjacency_list_BC< unordered_setBC, poolBC, directedS, int, int>,
  adjacency_list_BC< unordered_setBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< unordered_setBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< unordered_setBC, poolBC, undirectedS, int, int> > intint_adjlistBC_unordered_set_types;
  
typedef mpl::list< 
  adjacency_list_BC< multisetBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< multisetBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< multisetBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< multisetBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< multisetBC, listBC, directedS, int, int>,
  adjacency_list_BC< multisetBC, poolBC, directedS, int, int>,
  adjacency_list_BC< multisetBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< multisetBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< multisetBC, poolBC, undirectedS, int, int> > intint_adjlistBC_multiset_types;
  
typedef mpl::list< 
  adjacency_list_BC< unordered_multisetBC, vecBC,  bidirectionalS, int, int>,
  adjacency_list_BC< unordered_multisetBC, listBC, bidirectionalS, int, int>,
  adjacency_list_BC< unordered_multisetBC, poolBC, bidirectionalS, int, int>,
  adjacency_list_BC< unordered_multisetBC, vecBC,  directedS, int, int>,
  adjacency_list_BC< unordered_multisetBC, listBC, directedS, int, int>,
  adjacency_list_BC< unordered_multisetBC, poolBC, directedS, int, int>,
  adjacency_list_BC< unordered_multisetBC, vecBC,  undirectedS, int, int>,
  adjacency_list_BC< unordered_multisetBC, listBC, undirectedS, int, int>,
  adjacency_list_BC< unordered_multisetBC, poolBC, undirectedS, int, int> > intint_adjlistBC_unordered_multiset_types;

  

template <typename Seq1, typename Seq2>
struct join_seqs {
  typedef typename mpl::copy< Seq1, mpl::front_inserter< Seq2 > >::type type;
};


// These types are currently supported, but should be since having sets for storing out-edges could be useful:

// typedef mpl::list< adjacency_list_BC< setBC, vecBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< setBC, listBC, bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< setBC, poolBC, bidirectionalS, int, int> > intint_graphtest_types;
// 
// typedef mpl::list< adjacency_list_BC< multisetBC, vecBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< multisetBC, listBC, bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< multisetBC, poolBC, bidirectionalS, int, int> > intint_graphtest_types;
// 
// typedef mpl::list< adjacency_list_BC< unordered_setBC, vecBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< unordered_setBC, listBC, bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< unordered_setBC, poolBC, bidirectionalS, int, int> > intint_graphtest_types;
// 
// typedef mpl::list< adjacency_list_BC< unordered_multisetBC, vecBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< unordered_multisetBC, listBC, bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< unordered_multisetBC, poolBC, bidirectionalS, int, int> > intint_graphtest_types;



// These types should all trigger an "disallowed vertex list" compile-time error:

// typedef mpl::list< adjacency_list_BC< vecBC, setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, unordered_setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, unordered_multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, multisetBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, unordered_setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< vecBC, unordered_multisetBC,  directedS, int, int> > intint_graphtest_types;

// typedef mpl::list< adjacency_list_BC< listBC, setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, unordered_setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, unordered_multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, multisetBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, unordered_setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< listBC, unordered_multisetBC,  directedS, int, int> > intint_graphtest_types;

// typedef mpl::list< adjacency_list_BC< poolBC, setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, unordered_setBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, unordered_multisetBC,  bidirectionalS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, multisetBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, unordered_setBC,  directedS, int, int> > intint_graphtest_types;
// typedef mpl::list< adjacency_list_BC< poolBC, unordered_multisetBC,  directedS, int, int> > intint_graphtest_types;


// typedef join_seqs< intint_adjlistBC_nosets_types, intint_adjlistBC_nosets_undir_types >::type intint_graphtest_types;


typedef 
  join_seqs< intint_adjlistBC_vec_types,
  join_seqs< intint_adjlistBC_list_types, 
  join_seqs< intint_adjlistBC_pool_types, 
  join_seqs< intint_adjlistBC_set_types, 
  join_seqs< intint_adjlistBC_unordered_set_types,
  join_seqs< intint_adjlistBC_multiset_types, 
  join_seqs< intint_adjlistBC_unordered_multiset_types,
  intint_adjlist_types
  >::type >::type >::type >::type >::type >::type >::type intint_graphtest_types;

  


template <typename Graph>
typename enable_if< is_vertex_list_graph< Graph >,
void >::type check_graph_vertex_count(const Graph& g, std::size_t expected_count) {
  BOOST_CHECK_EQUAL( num_vertices(g), expected_count );
};

template <typename Graph>
typename disable_if< is_vertex_list_graph< Graph >,
void >::type check_graph_vertex_count(const Graph&, std::size_t) { };



template <typename Graph>
typename enable_if< is_vertex_list_graph< Graph >,
void >::type intint_check_graph_vertex_values(const Graph& g, std::size_t expected_count, const int* ref_values) {
  typedef typename graph_traits<Graph>::vertex_iterator VertexIter;
  VertexIter vi, vi_end;
  BOOST_CHECK_NO_THROW( tie(vi, vi_end) = vertices(g) );
  std::vector<int> vp_list;
  for(; vi != vi_end; ++vi)
    vp_list.push_back( g[*vi] );
  std::sort(vp_list.begin(), vp_list.end());
  BOOST_CHECK_EQUAL( vp_list.size(), expected_count );
  for(std::size_t i = 0; i < vp_list.size(); ++i) {
    BOOST_CHECK_EQUAL( vp_list[i], *(ref_values++) );
  };
};

template <typename Graph>
typename disable_if< is_vertex_list_graph< Graph >,
void >::type intint_check_graph_vertex_values(const Graph&, std::size_t, const int*) { };






template <typename Graph>
typename enable_if< is_edge_list_graph< Graph >,
void >::type check_graph_edge_count(const Graph& g, std::size_t expected_count) {
  BOOST_CHECK_EQUAL( num_edges(g), expected_count );
};

template <typename Graph>
typename disable_if< is_edge_list_graph< Graph >,
void >::type check_graph_edge_count(const Graph&, std::size_t) { };



template <typename Graph>
typename enable_if< is_edge_list_graph< Graph >,
void >::type intint_check_graph_edge_values(const Graph& g, std::size_t expected_count, const int* ref_values) {
  typedef typename graph_traits<Graph>::edge_iterator EdgeIter;
  EdgeIter ei, ei_end;
  BOOST_CHECK_NO_THROW( tie(ei, ei_end) = edges(g) );
  std::vector<int> ep_list;
  for(; ei != ei_end; ++ei)
    ep_list.push_back( g[*ei] );
  std::sort(ep_list.begin(), ep_list.end());
  BOOST_CHECK_EQUAL( ep_list.size(), expected_count );
  for(std::size_t i = 0; i < ep_list.size(); ++i) {
    BOOST_CHECK_EQUAL( ep_list[i], *(ref_values++) );
  };
};

template <typename Graph>
typename disable_if< is_edge_list_graph< Graph >,
void >::type intint_check_graph_edge_values(const Graph&, std::size_t, const int*) { };




template <typename Graph>
typename enable_if< is_bidirectional_graph< Graph >,
void >::type check_graph_in_degree(const Graph& g, 
                                   typename graph_traits<Graph>::vertex_descriptor v,
                                   std::size_t expected_count) {
  BOOST_CHECK_EQUAL( in_degree(v,g), expected_count );
};

template <typename Graph>
typename disable_if< is_bidirectional_graph< Graph >,
void >::type check_graph_in_degree(const Graph&, typename graph_traits<Graph>::vertex_descriptor, std::size_t) {};



template <typename Graph>
typename enable_if< is_bidirectional_graph< Graph >,
void >::type check_graph_in_edge_values(const Graph& g, 
                                        typename graph_traits<Graph>::vertex_descriptor v,
                                        std::size_t expected_count,
                                        const int* edge_values, const int* vertex_values) {
  typedef typename graph_traits<Graph>::in_edge_iterator InEdgeIter;
  BOOST_CHECK_EQUAL( in_degree(v,g), expected_count );
  InEdgeIter ei, ei_end;
  BOOST_CHECK_NO_THROW( tie(ei, ei_end) = in_edges(v, g) );
  
  std::vector<int> e_list;
  std::vector<int> vp_list;
  for(; ei != ei_end; ++ei) {
    BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
    e_list.push_back(g[*ei]);
    vp_list.push_back(g[source(*ei,g)]);
  };
  std::sort(e_list.begin(), e_list.end());
  for(std::size_t i = 0; i < e_list.size(); ++i)
    BOOST_CHECK_EQUAL( e_list[i], *(edge_values++) );
  std::sort(vp_list.begin(), vp_list.end());
  for(std::size_t i = 0; i < vp_list.size(); ++i)
    BOOST_CHECK_EQUAL( vp_list[i], *(vertex_values++) );
};

template <typename Graph>
typename disable_if< is_bidirectional_graph< Graph >,
void >::type check_graph_in_edge_values(const Graph&, 
                                        typename graph_traits<Graph>::vertex_descriptor,
                                        std::size_t, const int*, const int*) { };



template <typename Graph>
typename enable_if< is_incidence_graph< Graph >,
void >::type check_graph_out_degree(const Graph& g, 
                                   typename graph_traits<Graph>::vertex_descriptor u,
                                   std::size_t expected_count) {
  BOOST_CHECK_EQUAL( out_degree(u,g), expected_count );
};

template <typename Graph>
typename disable_if< is_incidence_graph< Graph >,
void >::type check_graph_out_degree(const Graph&, typename graph_traits<Graph>::vertex_descriptor, std::size_t) {};


template <typename Graph>
typename enable_if< is_incidence_graph< Graph >,
void >::type check_graph_out_edge_values(const Graph& g, 
                                         typename graph_traits<Graph>::vertex_descriptor u,
                                         std::size_t expected_count,
                                         const int* edge_values, const int* vertex_values) {
  typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIter;
  BOOST_CHECK_EQUAL( out_degree(u,g), expected_count );
  OutEdgeIter ei, ei_end;
  BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(u,g) );
  
  std::vector<int> e_list;
  std::vector<int> vp_list;
  for(; ei != ei_end; ++ei) {
    BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
    e_list.push_back(g[*ei]);
    vp_list.push_back(g[target(*ei,g)]);
  };
  std::sort(e_list.begin(), e_list.end());
  for(std::size_t i = 0; i < e_list.size(); ++i)
    BOOST_CHECK_EQUAL( e_list[i], *(edge_values++) );
  std::sort(vp_list.begin(), vp_list.end());
  for(std::size_t i = 0; i < vp_list.size(); ++i)
    BOOST_CHECK_EQUAL( vp_list[i], *(vertex_values++) );
};

template <typename Graph>
typename disable_if< is_incidence_graph< Graph >,
void >::type check_graph_out_edge_values(const Graph&, 
                                         typename graph_traits<Graph>::vertex_descriptor,
                                         std::size_t, const int*, const int*) { };




template <typename Graph>
void intint_check_fullbranch_integrity(const Graph& g, typename graph_traits<Graph>::vertex_descriptor u) {
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;
  typedef typename graph_traits<Graph>::out_edge_iterator OutEdgeIter;
  
  if(out_degree(u, g) == 0)
    return;
  
  BOOST_CHECK_EQUAL( out_degree(u, g), 4 );
  OutEdgeIter ei, ei_end;
  for(tie(ei,ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
    BOOST_CHECK_EQUAL( g[*ei], g[source(*ei, g)] * 1000 + g[target(*ei, g)]);
    intint_check_fullbranch_integrity(g, target(*ei, g));
  };
};





BOOST_AUTO_TEST_CASE_TEMPLATE( intint_bgl_mutable_graph_test, Graph, intint_graphtest_types )
{
  typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename graph_traits<Graph>::edge_descriptor Edge;
  
  Graph g;
  Vertex v_root = Vertex(); TEST_PRINT_REACHED_MARKER
  
  /* MutableGraph */
  BOOST_CHECK_NO_THROW( v_root = add_vertex(g) ); TEST_PRINT_REACHED_MARKER
  check_graph_vertex_count(g, 1); TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_NO_THROW( remove_vertex(v_root,g) ); TEST_PRINT_REACHED_MARKER
  check_graph_vertex_count(g, 0); TEST_PRINT_REACHED_MARKER
  
  /* MutablePropertyGraph */
  BOOST_CHECK_NO_THROW( v_root = add_vertex(1, g) ); TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_EQUAL( g[v_root], 1 ); TEST_PRINT_REACHED_MARKER
  g[v_root] = 1; TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_EQUAL( g[v_root], 1 ); TEST_PRINT_REACHED_MARKER
  
  /* MutableGraph */
  int vp_rc[] = {2,3,4,5};
  int ep_rc[] = {1002,1003,1004,1005};
  Vertex v_rc[4];
  Edge e_rc[4]; 
  for(int i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( v_rc[i] = add_vertex(g) ); TEST_PRINT_REACHED_MARKER
    g[ v_rc[i] ] = vp_rc[i];
    BOOST_CHECK_EQUAL( g[ v_rc[i] ], vp_rc[i] ); TEST_PRINT_REACHED_MARKER
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW( tie(e_rc[i],edge_added_success) = add_edge(v_root, v_rc[i], g) ); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK( edge_added_success ); TEST_PRINT_REACHED_MARKER
    g[ e_rc[i] ] = ep_rc[i];
    BOOST_CHECK_EQUAL( g[ e_rc[i] ], ep_rc[i] ); TEST_PRINT_REACHED_MARKER
  };
  check_graph_vertex_count(g, 5); TEST_PRINT_REACHED_MARKER
  
  /* MutablePropertyGraph */
  int vp_rc1c[] = {6,7,8,9};
  int ep_rc1c[] = {2006,2007,2008,2009};
  Vertex v_rc1c[4];
  Edge e_rc1c[4];
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( v_rc1c[i] = add_vertex(vp_rc1c[i], g) ); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK_EQUAL( g[ v_rc1c[i] ], vp_rc1c[i] ); TEST_PRINT_REACHED_MARKER
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW( tie(e_rc1c[i],edge_added_success) = add_edge(v_rc[0], v_rc1c[i], ep_rc1c[i], g) ); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK( edge_added_success ); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK_EQUAL( g[ e_rc1c[i] ], ep_rc1c[i] ); TEST_PRINT_REACHED_MARKER
  };
  check_graph_vertex_count(g, 9); TEST_PRINT_REACHED_MARKER
  
  
  BOOST_CHECK_EQUAL( g[v_root], 1 ); TEST_PRINT_REACHED_MARKER
  {
    /* IncidenceGraph */
    {
      const int e_vals[] = {1002, 1003, 1004, 1005};
      const int v_vals[] = {2, 3, 4, 5};
      check_graph_out_edge_values(g, v_root, 4, e_vals, v_vals); TEST_PRINT_REACHED_MARKER
    };
    
    /* BidirectionalGraph */
    {
      const int e_vals[] = {1002};
      const int v_vals[] = {1};
      check_graph_in_edge_values(g, v_rc[0], 1, e_vals, v_vals); TEST_PRINT_REACHED_MARKER
    };
    
    /* IncidenceGraph */
    {
      const int e_vals[] = {2006, 2007, 2008, 2009};
      const int v_vals[] = {6, 7, 8, 9};
      check_graph_out_edge_values(g, v_rc[0], 4, e_vals, v_vals); TEST_PRINT_REACHED_MARKER
    };
  };
  
  /* MutablePropertyGraph (with rvalue-ref) */
  int vp_rc2c[] = {10,11,12,13};
  int ep_rc2c[] = {3010,3011,3012,3013};
  Vertex v_rc2c[4];
  Edge e_rc2c[4];
  for(std::size_t i = 0; i < 4; ++i) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    BOOST_CHECK_NO_THROW( v_rc2c[i] = add_vertex(std::move(vp_rc2c[i]), g) ); TEST_PRINT_REACHED_MARKER
#else
    BOOST_CHECK_NO_THROW( v_rc2c[i] = add_vertex(vp_rc2c[i], g) ); TEST_PRINT_REACHED_MARKER
#endif
    bool edge_added_success = false;
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    BOOST_CHECK_NO_THROW( tie(e_rc2c[i],edge_added_success) = add_edge(v_rc[1], v_rc2c[i], ep_rc2c[i], g) ); TEST_PRINT_REACHED_MARKER
#else
    BOOST_CHECK_NO_THROW( tie(e_rc2c[i],edge_added_success) = add_edge(v_rc[1], v_rc2c[i], ep_rc2c[i], g) ); TEST_PRINT_REACHED_MARKER
#endif
    BOOST_CHECK( edge_added_success ); TEST_PRINT_REACHED_MARKER
  };
  check_graph_vertex_count(g, 13); TEST_PRINT_REACHED_MARKER
  
  {
    /* IncidenceGraph */
    const int e_vals[] = {3010, 3011, 3012, 3013};
    const int v_vals[] = {10, 11, 12, 13};
    check_graph_out_edge_values(g, v_rc[1], 4, e_vals, v_vals); TEST_PRINT_REACHED_MARKER
  };
  
  
  /* Copying function */
  intint_check_fullbranch_integrity(g, v_root);
  {
    Graph* p_g_cpy = NULL;
    BOOST_CHECK_NO_THROW( p_g_cpy = new Graph(g) ); TEST_PRINT_REACHED_MARKER
    intint_check_fullbranch_integrity(*p_g_cpy, *(vertices(*p_g_cpy).first)); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK_NO_THROW( delete p_g_cpy ); TEST_PRINT_REACHED_MARKER
  };
  
  {
    Graph g_cpy;
    BOOST_CHECK_NO_THROW( g_cpy = g ); TEST_PRINT_REACHED_MARKER
    intint_check_fullbranch_integrity(g_cpy, *(vertices(g_cpy).first)); TEST_PRINT_REACHED_MARKER
  };
  
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  {
    Graph* p_g_mv = NULL;
    BOOST_CHECK_NO_THROW( p_g_mv = new Graph(std::move(g)) ); TEST_PRINT_REACHED_MARKER
    intint_check_fullbranch_integrity(*p_g_mv, *(vertices(*p_g_mv).first)); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK_NO_THROW( g = std::move(*p_g_mv) ); TEST_PRINT_REACHED_MARKER
    intint_check_fullbranch_integrity(g, *(vertices(g).first)); TEST_PRINT_REACHED_MARKER
    BOOST_CHECK_NO_THROW( delete p_g_mv ); TEST_PRINT_REACHED_MARKER
    v_root = *(vertices(g).first); TEST_PRINT_REACHED_MARKER
  };
#endif
  
  /* MutableGraph */
  BOOST_CHECK_NO_THROW( clear_vertex(v_rc[0],g) ); TEST_PRINT_REACHED_MARKER
  
  /* IncidenceGraph */
  check_graph_out_degree(g, v_rc[0], 0); TEST_PRINT_REACHED_MARKER
  check_graph_out_degree(g, v_root, 3); TEST_PRINT_REACHED_MARKER
  
  /* BidirectionalGraph */
  check_graph_in_degree(g, v_rc[0], 0); TEST_PRINT_REACHED_MARKER
  check_graph_in_degree(g, v_rc1c[0], 0); TEST_PRINT_REACHED_MARKER
  check_graph_in_degree(g, v_rc1c[1], 0); TEST_PRINT_REACHED_MARKER
  check_graph_in_degree(g, v_rc1c[2], 0); TEST_PRINT_REACHED_MARKER
  check_graph_in_degree(g, v_rc1c[3], 0); TEST_PRINT_REACHED_MARKER
  
  /* VertexListGraph */
  check_graph_vertex_count(g, 13); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  /* EdgeListGraph */
  check_graph_edge_count(g, 7); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1003, 1004, 1005, 3010, 3011, 3012, 3013};
    intint_check_graph_edge_values(g, 7, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  
  
  /* MutableGraph */
  BOOST_CHECK_NO_THROW( remove_edge(v_rc[1], v_rc2c[2], g) ); TEST_PRINT_REACHED_MARKER
  
  /* VertexListGraph */
  check_graph_vertex_count(g, 13); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  /* EdgeListGraph */
  check_graph_edge_count(g, 6); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1003, 1004, 1005, 3010, 3011, 3013};
    intint_check_graph_edge_values(g, 6, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  
  
  /* MutableGraph */
  std::pair<Edge, bool> last_e_of_rc2;
  BOOST_CHECK_NO_THROW( last_e_of_rc2 = edge(v_rc[1], v_rc2c[3], g) ); TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_EQUAL( last_e_of_rc2.second, true ); TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_NO_THROW( remove_edge(last_e_of_rc2.first, g) ); TEST_PRINT_REACHED_MARKER
  
  /* VertexListGraph */
  check_graph_vertex_count(g, 13); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  /* EdgeListGraph */
  check_graph_edge_count(g, 5); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1003, 1004, 1005, 3010, 3011};
    intint_check_graph_edge_values(g, 5, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  
  
  /* MutableGraph */
  BOOST_CHECK_NO_THROW( clear_vertex(v_rc2c[0], g) ); TEST_PRINT_REACHED_MARKER
  BOOST_CHECK_NO_THROW( remove_vertex(v_rc2c[0], g) ); TEST_PRINT_REACHED_MARKER
  
  /* VertexListGraph */
  check_graph_vertex_count(g, 12); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13};
    intint_check_graph_vertex_values(g, 12, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  /* EdgeListGraph */
  check_graph_edge_count(g, 4); TEST_PRINT_REACHED_MARKER
  {
    const int ref_values[] = {1003, 1004, 1005, 3011};
    intint_check_graph_edge_values(g, 4, ref_values); TEST_PRINT_REACHED_MARKER
  };
  
  
};




