


BOOST_AUTO_TEST_CASE_TEMPLATE( PROPMAP_GRAPHTEST_NAME, Graph, PROPMAP_GRAPHTEST_TYPES )
{
  
  typedef typename graph_traits< Graph >::vertex_descriptor Vertex;
  typedef typename graph_traits< Graph >::edge_descriptor Edge;
  
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::vname_t >::type VNameMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::vdistance_t >::type VDistMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::ename_t >::type ENameMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::eweight_t >::type EWeightMap;
  
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::vname_t >::const_type VNameCMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::vdistance_t >::const_type VDistCMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::ename_t >::const_type ENameCMap;
  typedef typename property_map< Graph, PROPMAP_GRAPHTEST_MAPS::eweight_t >::const_type EWeightCMap;
  
  Graph g;
  
  VNameMap vname_map;
  BOOST_CHECK_NO_THROW( vname_map = get(PROPMAP_GRAPHTEST_MAPS::vname, g) );
  VDistMap dist_map;
  BOOST_CHECK_NO_THROW( dist_map = get(PROPMAP_GRAPHTEST_MAPS::vdistance, g) );
  ENameMap ename_map;
  BOOST_CHECK_NO_THROW( ename_map = get(PROPMAP_GRAPHTEST_MAPS::ename, g) );
  EWeightMap weight_map;
  BOOST_CHECK_NO_THROW( weight_map = get(PROPMAP_GRAPHTEST_MAPS::eweight, g) );
  
  const Graph& cg = g;
  
  VNameCMap vname_cmap;
  BOOST_CHECK_NO_THROW( vname_cmap = get(PROPMAP_GRAPHTEST_MAPS::vname, cg) );
  VDistCMap dist_cmap;
  BOOST_CHECK_NO_THROW( dist_cmap = get(PROPMAP_GRAPHTEST_MAPS::vdistance, cg) );
  ENameCMap ename_cmap;
  BOOST_CHECK_NO_THROW( ename_cmap = get(PROPMAP_GRAPHTEST_MAPS::ename, cg) );
  EWeightCMap weight_cmap;
  BOOST_CHECK_NO_THROW( weight_cmap = get(PROPMAP_GRAPHTEST_MAPS::eweight, cg) );
  
  
  Vertex v_root = create_root(g);
  Vertex v1;
  Edge e1;
  tie(v1, e1) = add_child_vertex(v_root,g);
  
  BOOST_CHECK_NO_THROW( put(vname_map, v_root, "root-vertex") );
  BOOST_CHECK_EQUAL( get(vname_map, v_root), "root-vertex" );
  BOOST_CHECK_EQUAL( get(vname_cmap, v_root), "root-vertex" );
  
  BOOST_CHECK_NO_THROW( put(dist_map, v_root, 42) );
  BOOST_CHECK_EQUAL( get(dist_map, v_root), 42 );
  BOOST_CHECK_EQUAL( get(dist_cmap, v_root), 42 );
  
  
  BOOST_CHECK_NO_THROW( put(ename_map, e1, "root-edge") );
  BOOST_CHECK_EQUAL( get(ename_map, e1), "root-edge" );
  BOOST_CHECK_EQUAL( get(ename_cmap, e1), "root-edge" );
  
  BOOST_CHECK_NO_THROW( put(weight_map, e1, 69) );
  BOOST_CHECK_EQUAL( get(weight_map, e1), 69 );
  BOOST_CHECK_EQUAL( get(weight_cmap, e1), 69 );
  
  
  BOOST_CHECK_NO_THROW( vname_map[v1] = "child-vertex" );
  BOOST_CHECK_EQUAL( vname_map[v1], "child-vertex" );
  BOOST_CHECK_EQUAL( vname_cmap[v1], "child-vertex" );
  BOOST_CHECK_EQUAL( get(vname_map, v_root), "root-vertex" );
  BOOST_CHECK_EQUAL( get(vname_cmap, v_root), "root-vertex" );
  
  BOOST_CHECK_NO_THROW( dist_map[v1] = 21 );
  BOOST_CHECK_EQUAL( dist_map[v1], 21 );
  BOOST_CHECK_EQUAL( dist_cmap[v1], 21 );
  BOOST_CHECK_EQUAL( get(dist_map, v_root), 42 );
  BOOST_CHECK_EQUAL( get(dist_cmap, v_root), 42 );
  
  
  BOOST_CHECK_NO_THROW( ename_map[e1] = "root-edge-again" );
  BOOST_CHECK_EQUAL( ename_map[e1], "root-edge-again" );
  BOOST_CHECK_EQUAL( ename_cmap[e1], "root-edge-again" );
  
  BOOST_CHECK_NO_THROW( weight_map[e1] = 96 );
  BOOST_CHECK_EQUAL( weight_map[e1], 96 );
  BOOST_CHECK_EQUAL( weight_cmap[e1], 96 );
  
};




