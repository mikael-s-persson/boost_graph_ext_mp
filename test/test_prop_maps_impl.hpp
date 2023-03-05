

TYPED_TEST_SUITE(PropMapGraphTest, PropMapGraphTestTypes);

TYPED_TEST(PropMapGraphTest, AllCases) {
  using Graph = TypeParam;
  using Vertex = typename graph_traits<Graph>::vertex_descriptor;
  using Edge = typename graph_traits<Graph>::edge_descriptor;

  using VNameMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::vname_t>::type;
  using VDistMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::vdistance_t>::type;
  using ENameMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::ename_t>::type;
  using EWeightMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::eweight_t>::type;

  using VNameCMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::vname_t>::const_type;
  using VDistCMap =
      typename property_map<Graph,
                            PROPMAP_GRAPHTEST_MAPS::vdistance_t>::const_type;
  using ENameCMap =
      typename property_map<Graph, PROPMAP_GRAPHTEST_MAPS::ename_t>::const_type;
  using EWeightCMap =
      typename property_map<Graph,
                            PROPMAP_GRAPHTEST_MAPS::eweight_t>::const_type;

  Graph g;

  VNameMap vname_map;
  EXPECT_NO_THROW(vname_map = get(PROPMAP_GRAPHTEST_MAPS::vname, g));
  VDistMap dist_map;
  EXPECT_NO_THROW(dist_map = get(PROPMAP_GRAPHTEST_MAPS::vdistance, g));
  ENameMap ename_map;
  EXPECT_NO_THROW(ename_map = get(PROPMAP_GRAPHTEST_MAPS::ename, g));
  EWeightMap weight_map;
  EXPECT_NO_THROW(weight_map = get(PROPMAP_GRAPHTEST_MAPS::eweight, g));

  const Graph& cg = g;

  VNameCMap vname_cmap;
  EXPECT_NO_THROW(vname_cmap = get(PROPMAP_GRAPHTEST_MAPS::vname, cg));
  VDistCMap dist_cmap;
  EXPECT_NO_THROW(dist_cmap = get(PROPMAP_GRAPHTEST_MAPS::vdistance, cg));
  ENameCMap ename_cmap;
  EXPECT_NO_THROW(ename_cmap = get(PROPMAP_GRAPHTEST_MAPS::ename, cg));
  EWeightCMap weight_cmap;
  EXPECT_NO_THROW(weight_cmap = get(PROPMAP_GRAPHTEST_MAPS::eweight, cg));

  Vertex v_root = create_root(g);
  Vertex v1;
  Edge e1;
  tie(v1, e1) = add_child_vertex(v_root, g);

  EXPECT_NO_THROW(put(vname_map, v_root, "root-vertex"));
  EXPECT_EQ(get(vname_map, v_root), "root-vertex");
  EXPECT_EQ(get(vname_cmap, v_root), "root-vertex");

  EXPECT_NO_THROW(put(dist_map, v_root, 42));
  EXPECT_EQ(get(dist_map, v_root), 42);
  EXPECT_EQ(get(dist_cmap, v_root), 42);

  EXPECT_NO_THROW(put(ename_map, e1, "root-edge"));
  EXPECT_EQ(get(ename_map, e1), "root-edge");
  EXPECT_EQ(get(ename_cmap, e1), "root-edge");

  EXPECT_NO_THROW(put(weight_map, e1, 69));
  EXPECT_EQ(get(weight_map, e1), 69);
  EXPECT_EQ(get(weight_cmap, e1), 69);

  EXPECT_NO_THROW(vname_map[v1] = "child-vertex");
  EXPECT_EQ(vname_map[v1], "child-vertex");
  EXPECT_EQ(vname_cmap[v1], "child-vertex");
  EXPECT_EQ(get(vname_map, v_root), "root-vertex");
  EXPECT_EQ(get(vname_cmap, v_root), "root-vertex");

  EXPECT_NO_THROW(dist_map[v1] = 21);
  EXPECT_EQ(dist_map[v1], 21);
  EXPECT_EQ(dist_cmap[v1], 21);
  EXPECT_EQ(get(dist_map, v_root), 42);
  EXPECT_EQ(get(dist_cmap, v_root), 42);

  EXPECT_NO_THROW(ename_map[e1] = "root-edge-again");
  EXPECT_EQ(ename_map[e1], "root-edge-again");
  EXPECT_EQ(ename_cmap[e1], "root-edge-again");

  EXPECT_NO_THROW(weight_map[e1] = 96);
  EXPECT_EQ(weight_map[e1], 96);
  EXPECT_EQ(weight_cmap[e1], 96);
}
