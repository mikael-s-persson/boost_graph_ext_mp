
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
#include <boost/graph/linked_tree_BC.hpp>
#include <boost/graph/pooled_adjacency_list.hpp>

// #include <boost/graph/d_ary_bf_tree.hpp>
// #include <boost/graph/d_ary_cob_tree.hpp>



#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE boost_graph_trees
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/front_inserter.hpp>


using namespace boost;


// typedef mpl::list< 
// //   d_ary_bf_tree<int, 4, int>, 
// //   d_ary_cob_tree<int, 4, int>, 
// 
// 
// // All of this works:
//   linked_tree_BC<vecBC, vecBC, bidirectionalS, int, int>,
//   linked_tree_BC<poolBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<listBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<setBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<multisetBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<unordered_setBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<unordered_multisetBC, vecBC, bidirectionalS, int, int>, 
//   linked_tree_BC<vecBC, listBC, bidirectionalS, int, int>, 
//   linked_tree_BC<poolBC, listBC, bidirectionalS, int, int>, 
//   linked_tree_BC<listBC, listBC, bidirectionalS, int, int>,
//   linked_tree_BC<setBC, listBC, bidirectionalS, int, int>, 
//   linked_tree_BC<multisetBC, listBC, bidirectionalS, int, int>, 
//   linked_tree_BC<unordered_setBC, listBC, bidirectionalS, int, int>,
//   linked_tree_BC<unordered_multisetBC, listBC, bidirectionalS, int, int>,
//   linked_tree_BC<vecBC, poolBC, bidirectionalS, int, int>, 
//   linked_tree_BC<poolBC, poolBC, bidirectionalS, int, int>, 
//   linked_tree_BC<listBC, poolBC, bidirectionalS, int, int>,
//   linked_tree_BC<setBC, poolBC, bidirectionalS, int, int>, 
//   linked_tree_BC<multisetBC, poolBC, bidirectionalS, int, int>, 
//   linked_tree_BC<unordered_setBC, poolBC, bidirectionalS, int, int>,
//   linked_tree_BC<unordered_multisetBC, poolBC, bidirectionalS, int, int>,
//   tree_storage<int, int>::type,
//   pooled_adjacency_list<bidirectionalS, int, int > > intint_treetest_types;
  
  
  
  
  
typedef mpl::list< 
//   d_ary_bf_tree<int, 4, int>, 
//   d_ary_cob_tree<int, 4, int>, 
  tree_storage<int, int>::type,
  pooled_adjacency_list<bidirectionalS, int, int > > intint_othertrees_types;
  
  
typedef mpl::list< 
  linked_tree_BC< vecBC,  vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< vecBC,  listBC, bidirectionalS, int, int>,
  linked_tree_BC< vecBC,  poolBC, bidirectionalS, int, int>,
  linked_tree_BC< vecBC,  vecBC,  directedS, int, int>,
  linked_tree_BC< vecBC,  listBC, directedS, int, int>,
  linked_tree_BC< vecBC,  poolBC, directedS, int, int> > intint_ltreeBC_vec_types;
  
typedef mpl::list< 
  linked_tree_BC< listBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< listBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< listBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< listBC, vecBC,  directedS, int, int>,
  linked_tree_BC< listBC, listBC, directedS, int, int>,
  linked_tree_BC< listBC, poolBC, directedS, int, int> > intint_ltreeBC_list_types;
  
typedef mpl::list< 
  linked_tree_BC< poolBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< poolBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< poolBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< poolBC, vecBC,  directedS, int, int>,
  linked_tree_BC< poolBC, listBC, directedS, int, int>,
  linked_tree_BC< poolBC, poolBC, directedS, int, int> > intint_ltreeBC_pool_types;
  
typedef mpl::list< 
  linked_tree_BC< setBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< setBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< setBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< setBC, vecBC,  directedS, int, int>,
  linked_tree_BC< setBC, listBC, directedS, int, int>,
  linked_tree_BC< setBC, poolBC, directedS, int, int> > intint_ltreeBC_set_types;
  
typedef mpl::list< 
  linked_tree_BC< unordered_setBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< unordered_setBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< unordered_setBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< unordered_setBC, vecBC,  directedS, int, int>,
  linked_tree_BC< unordered_setBC, listBC, directedS, int, int>,
  linked_tree_BC< unordered_setBC, poolBC, directedS, int, int> > intint_ltreeBC_unordered_set_types;
  
typedef mpl::list< 
  linked_tree_BC< multisetBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< multisetBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< multisetBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< multisetBC, vecBC,  directedS, int, int>,
  linked_tree_BC< multisetBC, listBC, directedS, int, int>,
  linked_tree_BC< multisetBC, poolBC, directedS, int, int> > intint_ltreeBC_multiset_types;
  
typedef mpl::list< 
  linked_tree_BC< unordered_multisetBC, vecBC,  bidirectionalS, int, int>,
  linked_tree_BC< unordered_multisetBC, listBC, bidirectionalS, int, int>,
  linked_tree_BC< unordered_multisetBC, poolBC, bidirectionalS, int, int>,
  linked_tree_BC< unordered_multisetBC, vecBC,  directedS, int, int>,
  linked_tree_BC< unordered_multisetBC, listBC, directedS, int, int>,
  linked_tree_BC< unordered_multisetBC, poolBC, directedS, int, int> > intint_ltreeBC_unordered_multiset_types;

  
template <typename Seq1, typename Seq2>
struct join_seqs {
  typedef typename mpl::copy< Seq1, mpl::front_inserter< Seq2 > >::type type;
};
  
  
typedef 
  join_seqs< intint_ltreeBC_vec_types,
  join_seqs< intint_ltreeBC_list_types, 
  join_seqs< intint_ltreeBC_pool_types, 
  join_seqs< intint_ltreeBC_set_types, 
  join_seqs< intint_ltreeBC_unordered_set_types,
  join_seqs< intint_ltreeBC_multiset_types, 
  join_seqs< intint_ltreeBC_unordered_multiset_types,
  intint_othertrees_types
  >::type >::type >::type >::type >::type >::type >::type intint_treetest_types;

  
  
  


template <typename TreeType>
typename enable_if< is_vertex_list_graph< TreeType >,
void >::type intint_do_final_vertex_check(const TreeType& g) {
  typedef typename graph_traits<TreeType>::vertex_iterator VertexIter;
  std::vector<int> all_vertices;
  VertexIter vi = VertexIter();
  VertexIter vi_end = VertexIter();
  BOOST_CHECK_NO_THROW( tie(vi,vi_end) = vertices(g) ); 
  for(; vi != vi_end; ++vi)
    all_vertices.push_back(g[*vi]);
  std::sort(all_vertices.begin(), all_vertices.end());
  BOOST_CHECK_EQUAL( all_vertices[0], 1 );
  BOOST_CHECK_EQUAL( all_vertices[1], 2 );
  BOOST_CHECK_EQUAL( all_vertices[2], 4 );
  BOOST_CHECK_EQUAL( all_vertices[3], 5 );
  BOOST_CHECK_EQUAL( all_vertices[4], 6 );
  BOOST_CHECK_EQUAL( all_vertices[5], 7 );
  BOOST_CHECK_EQUAL( all_vertices[6], 8 );
  BOOST_CHECK_EQUAL( all_vertices[7], 9 );
};

template <typename TreeType>
typename disable_if< is_vertex_list_graph< TreeType >,
void >::type intint_do_final_vertex_check(const TreeType&) { };



template <typename TreeType>
typename enable_if< is_vertex_list_graph< TreeType >,
void >::type intint_check_tree_vertex_count(const TreeType& g, std::size_t expected_count) {
  BOOST_CHECK_EQUAL( num_vertices(g), expected_count );
};

template <typename TreeType>
typename disable_if< is_vertex_list_graph< TreeType >,
void >::type intint_check_tree_vertex_count(const TreeType&, std::size_t) { };



template <typename TreeType>
typename enable_if< is_bidirectional_graph< TreeType >,
void >::type intint_do_in_edge_check(const TreeType& g, 
                                     typename graph_traits<TreeType>::vertex_descriptor v,
                                     int parent_value, int cur_value) {
  typedef typename graph_traits<TreeType>::vertex_descriptor Vertex;
  typedef typename graph_traits<TreeType>::in_edge_iterator InEdgeIter;
  
  Vertex u = Vertex();
  BOOST_CHECK_NO_THROW( u = parent_vertex(v, g) );
  BOOST_CHECK_EQUAL( g[u], parent_value );
  
  InEdgeIter iei = InEdgeIter();
  InEdgeIter iei_end = InEdgeIter();
  BOOST_CHECK_NO_THROW( tie(iei,iei_end) = in_edges(v,g) );
  for(; iei != iei_end; ++iei) {
    BOOST_CHECK_EQUAL( g[*iei], parent_value * 1000 + cur_value );
    BOOST_CHECK_EQUAL( g[source(*iei,g)], parent_value );
    BOOST_CHECK_EQUAL( g[target(*iei,g)], cur_value );
  };
};

template <typename TreeType>
typename disable_if< is_bidirectional_graph< TreeType >,
void >::type intint_do_in_edge_check(const TreeType&, typename graph_traits<TreeType>::vertex_descriptor, int, int) {};



template <typename TreeType>
typename enable_if< is_adjacency_matrix< TreeType >,
void >::type intint_do_edge_check(const TreeType& g, 
                                  typename graph_traits<TreeType>::vertex_descriptor u,
                                  typename graph_traits<TreeType>::vertex_descriptor v,
                                  int u_value, int v_value) {
  typedef typename graph_traits<TreeType>::edge_descriptor Edge;
  Edge e = Edge();
  BOOST_CHECK_NO_THROW( e = edge(u, v, g).first );
  BOOST_CHECK_EQUAL( g[e], u_value * 1000 + v_value );
  BOOST_CHECK_EQUAL( g[source(e, g)], u_value );
  BOOST_CHECK_EQUAL( g[target(e, g)], v_value );
};

template <typename TreeType>
typename disable_if< is_adjacency_matrix< TreeType >,
void >::type intint_do_edge_check(const TreeType&, typename graph_traits<TreeType>::vertex_descriptor, typename graph_traits<TreeType>::vertex_descriptor, int, int) {};




BOOST_AUTO_TEST_CASE_TEMPLATE( intint_treetest, TreeType, intint_treetest_types )
{
  typedef typename graph_traits<TreeType>::vertex_descriptor Vertex;
  typedef typename graph_traits<TreeType>::edge_descriptor Edge;
  typedef typename graph_traits<TreeType>::out_edge_iterator OutEdgeIter;
  typedef typename tree_traits<TreeType>::child_vertex_iterator ChildVertIter;
  
  TreeType g;
  
  /* MutableTree, Graph */
  Vertex v_root = graph_traits<TreeType>::null_vertex();
  BOOST_CHECK_NO_THROW( v_root = create_root(g) );
  intint_check_tree_vertex_count(g, 1);
  
  /* MutableTree */
  BOOST_CHECK_NO_THROW( remove_branch(v_root,g) );
  intint_check_tree_vertex_count(g, 0);
  
  /* MutableTree */
  int vp_r = 1;
  BOOST_CHECK_NO_THROW( v_root = create_root(vp_r, g) );
  intint_check_tree_vertex_count(g, 1);
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  
  /* MutablePropertyTree */
  std::vector<int> props;
  BOOST_CHECK_NO_THROW( remove_branch(v_root, back_inserter(props), g) );
  intint_check_tree_vertex_count(g, 0);
  BOOST_CHECK_EQUAL( props.size(), 1 );
  BOOST_CHECK_EQUAL( props[0], 1 );
  props.clear();
  
  /* MutablePropertyTree */
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  BOOST_CHECK_NO_THROW( v_root = create_root(std::move(vp_r), g) );
#else
  BOOST_CHECK_NO_THROW( v_root = create_root(vp_r, g) );
#endif
  intint_check_tree_vertex_count(g, 1);
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  
  /* MutablePropertyTree */
  int vp_rc[] = {2,3,4,5};
  int ep_rc[] = {1002,1003,1004,1005};
  Vertex v_rc[4];
  Edge e_rc[4];
  for(int i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( tie(v_rc[i],e_rc[i]) = add_child_vertex(v_root,vp_rc[i],ep_rc[i],g) );
  };
  intint_check_tree_vertex_count(g, 5);
  
  /* Tree(IncidenceGraph) */
  BOOST_CHECK_EQUAL( out_degree(v_root,g), 4 );   
  
  /* MutablePropertyTree */
  int vp_rc1c[] = {6,7,8,9};
  int ep_rc1c[] = {2006,2007,2008,2009};
  Vertex v_rc1c[4];
  Edge e_rc1c[4];
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g) );
  };
  intint_check_tree_vertex_count(g, 9);
  
  /* Tree */
  BOOST_CHECK_NO_THROW( v_root = get_root_vertex(g) );
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
      e_list.push_back(g[*ei]);
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1002);
    BOOST_CHECK_EQUAL( e_list[1], 1003);
    BOOST_CHECK_EQUAL( e_list[2], 1004);
    BOOST_CHECK_EQUAL( e_list[3], 1005);
    
    /* Tree */
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 2);
    BOOST_CHECK_EQUAL( vp_list[1], 3);
    BOOST_CHECK_EQUAL( vp_list[2], 4);
    BOOST_CHECK_EQUAL( vp_list[3], 5);
    
    /* Tree */
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if(g[*cvi] == 2) {
        
        /* Tree(IncidenceGraph) */
        BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
          e_list2.push_back(g[*ei]);
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 2006);
        BOOST_CHECK_EQUAL( e_list2[1], 2007);
        BOOST_CHECK_EQUAL( e_list2[2], 2008);
        BOOST_CHECK_EQUAL( e_list2[3], 2009);
        
      };
    };
  };
  
  /* MutablePropertyTree */
  int vp_rc2c[] = {10,11,12,13};
  int ep_rc2c[] = {3010,3011,3012,3013};
  Vertex v_rc2c[4];
  Edge e_rc2c[4];
  for(std::size_t i = 0; i < 4; ++i) {
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    BOOST_CHECK_NO_THROW( tie(v_rc2c[i],e_rc2c[i]) = add_child_vertex(v_rc[1],std::move(vp_rc2c[i]),std::move(ep_rc2c[i]),g) );
#else
    BOOST_CHECK_NO_THROW( tie(v_rc2c[i],e_rc2c[i]) = add_child_vertex(v_rc[1],vp_rc2c[i],ep_rc2c[i],g) );
#endif
  };
  intint_check_tree_vertex_count(g, 13);
  
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(v_rc[1],g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
      e_list.push_back(g[*ei]);
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 3010);
    BOOST_CHECK_EQUAL( e_list[1], 3011);
    BOOST_CHECK_EQUAL( e_list[2], 3012);
    BOOST_CHECK_EQUAL( e_list[3], 3013);
    
    /* Tree */
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_rc[1],g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 10);
    BOOST_CHECK_EQUAL( vp_list[1], 11);
    BOOST_CHECK_EQUAL( vp_list[2], 12);
    BOOST_CHECK_EQUAL( vp_list[3], 13);
  };
  
  /* MutableTree */
  BOOST_CHECK_NO_THROW( remove_branch(v_rc[0],g) );
  intint_check_tree_vertex_count(g, 8);
  
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
      e_list.push_back(g[*ei]);
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1003);
    BOOST_CHECK_EQUAL( e_list[1], 1004);
    BOOST_CHECK_EQUAL( e_list[2], 1005);
    
    /* Tree */
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 3);
    BOOST_CHECK_EQUAL( vp_list[1], 4);
    BOOST_CHECK_EQUAL( vp_list[2], 5);
    
    /* Tree */
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if(g[*cvi] == 3) {
        
        /* Tree(IncidenceGraph) */
        BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
          e_list2.push_back(g[*ei]);
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 3010);
        BOOST_CHECK_EQUAL( e_list2[1], 3011);
        BOOST_CHECK_EQUAL( e_list2[2], 3012);
        BOOST_CHECK_EQUAL( e_list2[3], 3013);
        
      };
    };
  };
  
  /* MutablePropertyTree */
  BOOST_CHECK_NO_THROW( remove_branch(v_rc[1], back_inserter(props), g) );
  BOOST_CHECK_EQUAL( props.size(), 5 );
  BOOST_CHECK_EQUAL( props[0], 3);  // the first vertex should be the root of the branch.
  std::sort(props.begin(), props.end());
  BOOST_CHECK_EQUAL( props[1], 10);
  BOOST_CHECK_EQUAL( props[2], 11);
  BOOST_CHECK_EQUAL( props[3], 12);
  BOOST_CHECK_EQUAL( props[4], 13);
  intint_check_tree_vertex_count(g, 3);
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
      e_list.push_back(g[*ei]);
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1004);
    BOOST_CHECK_EQUAL( e_list[1], 1005);
    
    /* Tree */
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 4);
    BOOST_CHECK_EQUAL( vp_list[1], 5);
    
  };
  
  
  /* MutablePropertyTree */
  BOOST_CHECK_NO_THROW( tie(v_rc[0],e_rc[0]) = add_child_vertex(v_root,vp_rc[0],ep_rc[0],g) );
  for(std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW( tie(v_rc1c[i],e_rc1c[i]) = add_child_vertex(v_rc[0],vp_rc1c[i],ep_rc1c[i],g) );
  };
  intint_check_tree_vertex_count(g, 8);
  
  
  BOOST_CHECK_EQUAL( g[v_root], 1 );
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei, ei_end;
    BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(v_root,g) );
    std::vector<int> e_list;
    for(; ei != ei_end; ++ei) {
      BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
      e_list.push_back(g[*ei]);
    };
    std::sort(e_list.begin(), e_list.end());
    BOOST_CHECK_EQUAL( e_list[0], 1002);
    BOOST_CHECK_EQUAL( e_list[1], 1004);
    BOOST_CHECK_EQUAL( e_list[2], 1005);
    
    /* Tree */
    ChildVertIter cvi, cvi_end;
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector<int> vp_list;
    for(; cvi != cvi_end; ++cvi)
      vp_list.push_back(g[*cvi]);
    std::sort(vp_list.begin(), vp_list.end());
    BOOST_CHECK_EQUAL( vp_list[0], 2);
    BOOST_CHECK_EQUAL( vp_list[1], 4);
    BOOST_CHECK_EQUAL( vp_list[2], 5);
    
    /* Tree */
    BOOST_CHECK_NO_THROW( tie(cvi,cvi_end) = child_vertices(v_root,g) );
    std::vector< Vertex > v_list;
    for(; cvi != cvi_end; ++cvi) {
      if(g[*cvi] == 2) {
        
        /* Tree(IncidenceGraph) */
        BOOST_CHECK_NO_THROW( tie(ei,ei_end) = out_edges(*cvi,g) );
        std::vector<int> e_list2;
        for(; ei != ei_end; ++ei) {
          BOOST_CHECK_EQUAL( g[*ei], (g[source(*ei,g)] * 1000 + g[target(*ei,g)]) );
          e_list2.push_back(g[*ei]);
        };
        std::sort(e_list2.begin(), e_list2.end());
        BOOST_CHECK_EQUAL( e_list2[0], 2006);
        BOOST_CHECK_EQUAL( e_list2[1], 2007);
        BOOST_CHECK_EQUAL( e_list2[2], 2008);
        BOOST_CHECK_EQUAL( e_list2[3], 2009);
        
      };
    };
  };
  
  
  /* VertexListGraph */
  intint_do_final_vertex_check(g);
  /* BidirectionalGraph */
  intint_do_in_edge_check(g, v_rc1c[1], 2, 7);
  /* AdjacencyMatrix */
  intint_do_edge_check(g, v_rc[0], v_rc1c[2], 2, 8);
  
};




