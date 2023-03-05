
// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define PRINT_LINE_NUM_REACHED \
  std::cout << __LINE__ << " reached!" << std::endl;

template <typename TreeType>
std::enable_if_t<is_vertex_list_graph<TreeType>::value>
intint_do_final_vertex_check(const TreeType &g) {
  using VertexIter = typename graph_traits<TreeType>::vertex_iterator;
  std::vector<int> all_vertices;
  VertexIter vi = VertexIter();
  VertexIter vi_end = VertexIter();
  EXPECT_NO_THROW(tie(vi, vi_end) = vertices(g));
  for (; vi != vi_end; ++vi) {
    all_vertices.push_back(g[*vi]);
  }
  std::sort(all_vertices.begin(), all_vertices.end());
  EXPECT_EQ(all_vertices[0], 1);
  EXPECT_EQ(all_vertices[1], 2);
  EXPECT_EQ(all_vertices[2], 4);
  EXPECT_EQ(all_vertices[3], 5);
  EXPECT_EQ(all_vertices[4], 6);
  EXPECT_EQ(all_vertices[5], 7);
  EXPECT_EQ(all_vertices[6], 8);
  EXPECT_EQ(all_vertices[7], 9);
}

template <typename TreeType>
std::enable_if_t<!is_vertex_list_graph<TreeType>::value>
intint_do_final_vertex_check(const TreeType & /*unused*/) {}

template <typename TreeType>
std::enable_if_t<is_vertex_list_graph<TreeType>::value>
intint_check_tree_vertex_count(const TreeType &g, std::size_t expected_count) {
  EXPECT_EQ(num_vertices(g), expected_count);
}

template <typename TreeType>
std::enable_if_t<!is_vertex_list_graph<TreeType>::value>
intint_check_tree_vertex_count(const TreeType & /*unused*/,
                               std::size_t /*unused*/) {}

template <typename TreeType>
std::enable_if_t<is_bidirectional_graph<TreeType>::value>
intint_do_in_edge_check(const TreeType &g,
                        typename graph_traits<TreeType>::vertex_descriptor v,
                        int parent_value, int cur_value) {
  using Vertex = typename graph_traits<TreeType>::vertex_descriptor;
  using InEdgeIter = typename graph_traits<TreeType>::in_edge_iterator;

  Vertex u = Vertex();
  EXPECT_NO_THROW(u = parent_vertex(v, g));
  EXPECT_EQ(g[u], parent_value);

  InEdgeIter iei = InEdgeIter();
  InEdgeIter iei_end = InEdgeIter();
  EXPECT_NO_THROW(tie(iei, iei_end) = in_edges(v, g));
  for (; iei != iei_end; ++iei) {
    EXPECT_EQ(g[*iei], parent_value * 1000 + cur_value);
    EXPECT_EQ(g[source(*iei, g)], parent_value);
    EXPECT_EQ(g[target(*iei, g)], cur_value);
  }
}

template <typename TreeType>
std::enable_if_t<!is_bidirectional_graph<TreeType>::value>
intint_do_in_edge_check(
    const TreeType & /*unused*/,
    typename graph_traits<TreeType>::vertex_descriptor /*unused*/,
    int /*unused*/, int /*unused*/) {}

template <typename TreeType>
std::enable_if_t<is_adjacency_matrix<TreeType>::value>
intint_do_edge_check(const TreeType &g,
                     typename graph_traits<TreeType>::vertex_descriptor u,
                     typename graph_traits<TreeType>::vertex_descriptor v,
                     int u_value, int v_value) {
  using Edge = typename graph_traits<TreeType>::edge_descriptor;
  Edge e = Edge();
  EXPECT_NO_THROW(e = edge(u, v, g).first);
  EXPECT_EQ(g[e], u_value * 1000 + v_value);
  EXPECT_EQ(g[source(e, g)], u_value);
  EXPECT_EQ(g[target(e, g)], v_value);
}

template <typename TreeType>
std::enable_if_t<!is_adjacency_matrix<TreeType>::value> intint_do_edge_check(
    const TreeType & /*unused*/,
    typename graph_traits<TreeType>::vertex_descriptor /*unused*/,
    typename graph_traits<TreeType>::vertex_descriptor /*unused*/,
    int /*unused*/, int /*unused*/) {}

template <typename TreeType>
void intint_check_fullbranch_integrity(
    const TreeType& g, typename graph_traits<TreeType>::vertex_descriptor u) {
  using OutEdgeIter = typename graph_traits<TreeType>::out_edge_iterator;

  if (out_degree(u, g) == 0) {
    return;
  }

  EXPECT_EQ(out_degree(u, g), 4);
  OutEdgeIter ei;
  OutEdgeIter ei_end;
  for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
    EXPECT_EQ(g[*ei], g[source(*ei, g)] * 1000 + g[target(*ei, g)]);
    intint_check_fullbranch_integrity(g, target(*ei, g));
  }
}

TYPED_TEST_SUITE(IntIntTreeTest, IntIntTreeTestTypes);

TYPED_TEST(IntIntTreeTest, AllCases) {
  using TreeType = TypeParam;
  using Vertex = typename graph_traits<TreeType>::vertex_descriptor;
  using Edge = typename graph_traits<TreeType>::edge_descriptor;
  using OutEdgeIter = typename graph_traits<TreeType>::out_edge_iterator;
  using ChildVertIter = typename tree_traits<TreeType>::child_vertex_iterator;

  TreeType g;

  /* MutableTree, Graph */
  Vertex v_root = graph_traits<TreeType>::null_vertex();
  EXPECT_NO_THROW(v_root = create_root(g));
  intint_check_tree_vertex_count(g, 1);

  /* MutableTree */
  EXPECT_NO_THROW(remove_branch(v_root, g));
  intint_check_tree_vertex_count(g, 0);

  /* MutableTree */
  int vp_r = 1;
  EXPECT_NO_THROW(v_root = create_root(vp_r, g));
  intint_check_tree_vertex_count(g, 1);
  EXPECT_EQ(g[v_root], 1);

  /* MutablePropertyTree */
  std::vector<int> props;
  EXPECT_NO_THROW(remove_branch(v_root, back_inserter(props), g));
  intint_check_tree_vertex_count(g, 0);
  EXPECT_EQ(props.size(), 1);
  EXPECT_EQ(props[0], 1);
  props.clear();

  /* MutablePropertyTree */
  EXPECT_NO_THROW(v_root = create_root(std::move(vp_r), g));
  intint_check_tree_vertex_count(g, 1);
  EXPECT_EQ(g[v_root], 1);

  /* MutablePropertyTree */
  std::array<int, 4> vp_rc = {2, 3, 4, 5};
  std::array<int, 4> ep_rc = {1002, 1003, 1004, 1005};
  std::array<Vertex, 4> v_rc;
  std::array<Edge, 4> e_rc;
  for (int i = 0; i < 4; ++i) {
    EXPECT_NO_THROW(tie(v_rc[i], e_rc[i]) =
                        add_child_vertex(v_root, vp_rc[i], ep_rc[i], g));
  }
  intint_check_tree_vertex_count(g, 5);

  /* Tree(IncidenceGraph) */
  EXPECT_EQ(out_degree(v_root, g), 4);

  /* MutablePropertyTree */
  std::array<int, 4> vp_rc1c = {6, 7, 8, 9};
  std::array<int, 4> ep_rc1c = {2006, 2007, 2008, 2009};
  std::array<Vertex, 4> v_rc1c;
  std::array<Edge, 4> e_rc1c;
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_NO_THROW(tie(v_rc1c[i], e_rc1c[i]) =
                        add_child_vertex(v_rc[0], vp_rc1c[i], ep_rc1c[i], g));
  }
  intint_check_tree_vertex_count(g, 9);

  /* Tree */
  EXPECT_NO_THROW(v_root = get_root_vertex(g));
  EXPECT_EQ(g[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(v_root, g));
    std::vector<int> e_list;
    for (; ei != ei_end; ++ei) {
      EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
      e_list.push_back(g[*ei]);
    }
    std::sort(e_list.begin(), e_list.end());
    EXPECT_EQ(e_list[0], 1002);
    EXPECT_EQ(e_list[1], 1003);
    EXPECT_EQ(e_list[2], 1004);
    EXPECT_EQ(e_list[3], 1005);

    /* Tree */
    ChildVertIter cvi;
    ChildVertIter cvi_end;
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<int> vp_list;
    for (; cvi != cvi_end; ++cvi) {
      vp_list.push_back(g[*cvi]);
    }
    std::sort(vp_list.begin(), vp_list.end());
    EXPECT_EQ(vp_list[0], 2);
    EXPECT_EQ(vp_list[1], 3);
    EXPECT_EQ(vp_list[2], 4);
    EXPECT_EQ(vp_list[3], 5);

    /* Tree */
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<Vertex> v_list;
    for (; cvi != cvi_end; ++cvi) {
      if (g[*cvi] == 2) {

        /* Tree(IncidenceGraph) */
        EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(*cvi, g));
        std::vector<int> e_list2;
        for (; ei != ei_end; ++ei) {
          EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
          e_list2.push_back(g[*ei]);
        }
        std::sort(e_list2.begin(), e_list2.end());
        EXPECT_EQ(e_list2[0], 2006);
        EXPECT_EQ(e_list2[1], 2007);
        EXPECT_EQ(e_list2[2], 2008);
        EXPECT_EQ(e_list2[3], 2009);
      }
    }
  }

  /* MutablePropertyTree */
  std::array<int, 4> vp_rc2c = {10, 11, 12, 13};
  std::array<int, 4> ep_rc2c = {3010, 3011, 3012, 3013};
  std::array<Vertex, 4> v_rc2c;
  std::array<Edge, 4> e_rc2c;
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_NO_THROW(tie(v_rc2c[i], e_rc2c[i]) =
                        add_child_vertex(v_rc[1], std::move(vp_rc2c[i]),
                                         std::move(ep_rc2c[i]), g));
  }
  intint_check_tree_vertex_count(g, 13);

  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(v_rc[1], g));
    std::vector<int> e_list;
    for (; ei != ei_end; ++ei) {
      EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
      e_list.push_back(g[*ei]);
    }
    std::sort(e_list.begin(), e_list.end());
    EXPECT_EQ(e_list[0], 3010);
    EXPECT_EQ(e_list[1], 3011);
    EXPECT_EQ(e_list[2], 3012);
    EXPECT_EQ(e_list[3], 3013);

    /* Tree */
    ChildVertIter cvi;
    ChildVertIter cvi_end;
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_rc[1], g));
    std::vector<int> vp_list;
    for (; cvi != cvi_end; ++cvi) {
      vp_list.push_back(g[*cvi]);
    }
    std::sort(vp_list.begin(), vp_list.end());
    EXPECT_EQ(vp_list[0], 10);
    EXPECT_EQ(vp_list[1], 11);
    EXPECT_EQ(vp_list[2], 12);
    EXPECT_EQ(vp_list[3], 13);
  }

  /* Copying function */
  intint_check_fullbranch_integrity(g, v_root);
  {
    TreeType* p_g_cpy = nullptr;
    EXPECT_NO_THROW(p_g_cpy = new TreeType(g));
    intint_check_fullbranch_integrity(*p_g_cpy, get_root_vertex(*p_g_cpy));
    EXPECT_NO_THROW(delete p_g_cpy);
  }

  {
    TreeType g_cpy;
    EXPECT_NO_THROW(g_cpy = g);
    intint_check_fullbranch_integrity(g_cpy, get_root_vertex(g_cpy));
  }

  {
    TreeType* p_g_mv = nullptr;
    EXPECT_NO_THROW(p_g_mv = new TreeType(std::move(g)));
    intint_check_fullbranch_integrity(*p_g_mv, get_root_vertex(*p_g_mv));
    EXPECT_NO_THROW(g = std::move(*p_g_mv));
    intint_check_fullbranch_integrity(g, get_root_vertex(g));
    EXPECT_NO_THROW(delete p_g_mv);
    v_root = get_root_vertex(g);
  }

  /* MutableTree */
  EXPECT_NO_THROW(remove_branch(v_rc[0], g));
  intint_check_tree_vertex_count(g, 8);

  EXPECT_EQ(g[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(v_root, g));
    std::vector<int> e_list;
    for (; ei != ei_end; ++ei) {
      EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
      e_list.push_back(g[*ei]);
    }
    std::sort(e_list.begin(), e_list.end());
    EXPECT_EQ(e_list[0], 1003);
    EXPECT_EQ(e_list[1], 1004);
    EXPECT_EQ(e_list[2], 1005);

    /* Tree */
    ChildVertIter cvi;
    ChildVertIter cvi_end;
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<int> vp_list;
    for (; cvi != cvi_end; ++cvi) {
      vp_list.push_back(g[*cvi]);
    }
    std::sort(vp_list.begin(), vp_list.end());
    EXPECT_EQ(vp_list[0], 3);
    EXPECT_EQ(vp_list[1], 4);
    EXPECT_EQ(vp_list[2], 5);

    /* Tree */
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<Vertex> v_list;
    for (; cvi != cvi_end; ++cvi) {
      if (g[*cvi] == 3) {

        /* Tree(IncidenceGraph) */
        EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(*cvi, g));
        std::vector<int> e_list2;
        for (; ei != ei_end; ++ei) {
          EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
          e_list2.push_back(g[*ei]);
        }
        std::sort(e_list2.begin(), e_list2.end());
        EXPECT_EQ(e_list2[0], 3010);
        EXPECT_EQ(e_list2[1], 3011);
        EXPECT_EQ(e_list2[2], 3012);
        EXPECT_EQ(e_list2[3], 3013);
      }
    }
  }

  /* MutablePropertyTree */
  EXPECT_NO_THROW(remove_branch(v_rc[1], back_inserter(props), g));
  EXPECT_EQ(props.size(), 5);
  EXPECT_EQ(props[0],
            3); // the first vertex should be the root of the branch.
  std::sort(props.begin(), props.end());
  EXPECT_EQ(props[1], 10);
  EXPECT_EQ(props[2], 11);
  EXPECT_EQ(props[3], 12);
  EXPECT_EQ(props[4], 13);
  intint_check_tree_vertex_count(g, 3);

  EXPECT_EQ(g[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(v_root, g));
    std::vector<int> e_list;
    for (; ei != ei_end; ++ei) {
      EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
      e_list.push_back(g[*ei]);
    }
    std::sort(e_list.begin(), e_list.end());
    EXPECT_EQ(e_list[0], 1004);
    EXPECT_EQ(e_list[1], 1005);

    /* Tree */
    ChildVertIter cvi;
    ChildVertIter cvi_end;
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<int> vp_list;
    for (; cvi != cvi_end; ++cvi) {
      vp_list.push_back(g[*cvi]);
    }
    std::sort(vp_list.begin(), vp_list.end());
    EXPECT_EQ(vp_list[0], 4);
    EXPECT_EQ(vp_list[1], 5);
  }

  /* MutablePropertyTree */
  EXPECT_NO_THROW(tie(v_rc[0], e_rc[0]) =
                      add_child_vertex(v_root, vp_rc[0], ep_rc[0], g));
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_NO_THROW(tie(v_rc1c[i], e_rc1c[i]) =
                        add_child_vertex(v_rc[0], vp_rc1c[i], ep_rc1c[i], g));
  }
  intint_check_tree_vertex_count(g, 8);

  EXPECT_EQ(g[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(v_root, g));
    std::vector<int> e_list;
    for (; ei != ei_end; ++ei) {
      EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
      e_list.push_back(g[*ei]);
    }
    std::sort(e_list.begin(), e_list.end());
    EXPECT_EQ(e_list[0], 1002);
    EXPECT_EQ(e_list[1], 1004);
    EXPECT_EQ(e_list[2], 1005);

    /* Tree */
    ChildVertIter cvi;
    ChildVertIter cvi_end;
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<int> vp_list;
    for (; cvi != cvi_end; ++cvi) {
      vp_list.push_back(g[*cvi]);
    }
    std::sort(vp_list.begin(), vp_list.end());
    EXPECT_EQ(vp_list[0], 2);
    EXPECT_EQ(vp_list[1], 4);
    EXPECT_EQ(vp_list[2], 5);

    /* Tree */
    EXPECT_NO_THROW(tie(cvi, cvi_end) = child_vertices(v_root, g));
    std::vector<Vertex> v_list;
    for (; cvi != cvi_end; ++cvi) {
      if (g[*cvi] == 2) {

        /* Tree(IncidenceGraph) */
        EXPECT_NO_THROW(tie(ei, ei_end) = out_edges(*cvi, g));
        std::vector<int> e_list2;
        for (; ei != ei_end; ++ei) {
          EXPECT_EQ(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
          e_list2.push_back(g[*ei]);
        }
        std::sort(e_list2.begin(), e_list2.end());
        EXPECT_EQ(e_list2[0], 2006);
        EXPECT_EQ(e_list2[1], 2007);
        EXPECT_EQ(e_list2[2], 2008);
        EXPECT_EQ(e_list2[3], 2009);
      }
    }
  }

  /* VertexListGraph */
  intint_do_final_vertex_check(g);
  /* BidirectionalGraph */
  intint_do_in_edge_check(g, v_rc1c[1], 2, 7);
  /* AdjacencyMatrix */
  intint_do_edge_check(g, v_rc[0], v_rc1c[2], 2, 8);
}
