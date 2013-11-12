
// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#define PRINT_LINE_NUM_REACHED std::cout << __LINE__ << " reached!" << std::endl;


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



template <typename TreeType>
void intint_check_fullbranch_integrity(const TreeType& g, typename graph_traits<TreeType>::vertex_descriptor u) {
  typedef typename graph_traits<TreeType>::out_edge_iterator OutEdgeIter;
  
  if(out_degree(u, g) == 0)
    return;
  
  BOOST_CHECK_EQUAL( out_degree(u, g), 4 );
  OutEdgeIter ei, ei_end;
  for(tie(ei,ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
    BOOST_CHECK_EQUAL( g[*ei], g[source(*ei, g)] * 1000 + g[target(*ei, g)]);
    intint_check_fullbranch_integrity(g, target(*ei, g));
  };
};





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
  
  /* Copying function */
  intint_check_fullbranch_integrity(g, v_root); 
  {
    TreeType* p_g_cpy = NULL;
    BOOST_CHECK_NO_THROW( p_g_cpy = new TreeType(g) ); 
    intint_check_fullbranch_integrity(*p_g_cpy, get_root_vertex(*p_g_cpy)); 
    BOOST_CHECK_NO_THROW( delete p_g_cpy ); 
  };
  
  {
    TreeType g_cpy; 
    BOOST_CHECK_NO_THROW( g_cpy = g ); 
    intint_check_fullbranch_integrity(g_cpy, get_root_vertex(g_cpy)); 
  };
  
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
  {
    TreeType* p_g_mv = NULL;
    BOOST_CHECK_NO_THROW( p_g_mv = new TreeType(std::move(g)) ); 
    intint_check_fullbranch_integrity(*p_g_mv, get_root_vertex(*p_g_mv)); 
    BOOST_CHECK_NO_THROW( g = std::move(*p_g_mv) ); 
    intint_check_fullbranch_integrity(g, get_root_vertex(g)); 
    BOOST_CHECK_NO_THROW( delete p_g_mv ); 
    v_root = get_root_vertex(g); 
  };
#endif
  
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




