// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#define PRINT_LINE_NUM_REACHED \
  std::cout << __LINE__ << " reached!" << std::endl;

template <typename Graph>
typename enable_if<is_vertex_list_graph<Graph>, void>::type
check_graph_vertex_count(const Graph& g, std::size_t expected_count) {
  BOOST_CHECK_EQUAL(num_vertices(g), expected_count);
}

template <typename Graph>
typename disable_if<is_vertex_list_graph<Graph>, void>::type
check_graph_vertex_count(const Graph& /*unused*/, std::size_t /*unused*/) {}

template <typename Graph>
typename enable_if<is_vertex_list_graph<Graph>, void>::type
intint_check_graph_vertex_values(const Graph& g, std::size_t expected_count,
                                 const int* ref_values) {
  using VertexIter = typename graph_traits<Graph>::vertex_iterator;
  VertexIter vi;
  VertexIter vi_end;
  BOOST_CHECK_NO_THROW(tie(vi, vi_end) = vertices(g));
  std::vector<int> vp_list;
  for (; vi != vi_end; ++vi) {
    vp_list.push_back(g[*vi]);
  }
  std::sort(vp_list.begin(), vp_list.end());
  BOOST_CHECK_EQUAL(vp_list.size(), expected_count);
  for (int& i : vp_list) {
    BOOST_CHECK_EQUAL(i, *(ref_values++));
  }
}

template <typename Graph>
typename disable_if<is_vertex_list_graph<Graph>, void>::type
intint_check_graph_vertex_values(const Graph& /*unused*/,
                                 std::size_t /*unused*/,
                                 const int* /*unused*/) {}

template <typename Graph>
typename enable_if<is_edge_list_graph<Graph>, void>::type
check_graph_edge_count(const Graph& g, std::size_t expected_count) {
  BOOST_CHECK_EQUAL(num_edges(g), expected_count);
}

template <typename Graph>
typename disable_if<is_edge_list_graph<Graph>, void>::type
check_graph_edge_count(const Graph& /*unused*/, std::size_t /*unused*/) {}

template <typename Graph>
typename enable_if<is_edge_list_graph<Graph>, void>::type
intint_check_graph_edge_values(const Graph& g, std::size_t expected_count,
                               const int* ref_values) {
  using EdgeIter = typename graph_traits<Graph>::edge_iterator;
  EdgeIter ei;
  EdgeIter ei_end;
  BOOST_CHECK_NO_THROW(tie(ei, ei_end) = edges(g));
  std::vector<int> ep_list;
  for (; ei != ei_end; ++ei) {
    ep_list.push_back(g[*ei]);
  }
  std::sort(ep_list.begin(), ep_list.end());
  BOOST_CHECK_EQUAL(ep_list.size(), expected_count);
  for (int& i : ep_list) {
    BOOST_CHECK_EQUAL(i, *(ref_values++));
  }
}

template <typename Graph>
typename disable_if<is_edge_list_graph<Graph>, void>::type
intint_check_graph_edge_values(const Graph& /*unused*/, std::size_t /*unused*/,
                               const int* /*unused*/) {}

template <typename Graph>
typename enable_if<is_bidirectional_graph<Graph>, void>::type
check_graph_in_degree(const Graph& g,
                      typename graph_traits<Graph>::vertex_descriptor v,
                      std::size_t expected_count) {
  BOOST_CHECK_EQUAL(in_degree(v, g), expected_count);
}

template <typename Graph>
typename disable_if<is_bidirectional_graph<Graph>, void>::type
check_graph_in_degree(
    const Graph& /*unused*/,
    typename graph_traits<Graph>::vertex_descriptor /*unused*/,
    std::size_t /*unused*/) {}

template <typename Graph>
typename enable_if<is_bidirectional_graph<Graph>, void>::type
check_graph_in_edge_values(const Graph& g,
                           typename graph_traits<Graph>::vertex_descriptor v,
                           std::size_t expected_count, const int* edge_values,
                           const int* vertex_values) {
  using InEdgeIter = typename graph_traits<Graph>::in_edge_iterator;
  BOOST_CHECK_EQUAL(in_degree(v, g), expected_count);
  InEdgeIter ei;
  InEdgeIter ei_end;
  BOOST_CHECK_NO_THROW(tie(ei, ei_end) = in_edges(v, g));

  std::vector<int> e_list;
  std::vector<int> vp_list;
  for (; ei != ei_end; ++ei) {
    if (is_undirected_graph<Graph>::type::value) {
      BOOST_CHECK((g[*ei] == g[source(*ei, g)] * 1000 + g[target(*ei, g)]) ||
                  (g[*ei] == g[source(*ei, g)] + g[target(*ei, g)] * 1000));
    } else {
      BOOST_CHECK_EQUAL(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
    }
    e_list.push_back(g[*ei]);
    vp_list.push_back(g[source(*ei, g)]);
  }
  std::sort(e_list.begin(), e_list.end());
  for (int& i : e_list) {
    BOOST_CHECK_EQUAL(i, *(edge_values++));
  }
  std::sort(vp_list.begin(), vp_list.end());
  for (int& i : vp_list) {
    BOOST_CHECK_EQUAL(i, *(vertex_values++));
  }
}

template <typename Graph>
typename disable_if<is_bidirectional_graph<Graph>, void>::type
check_graph_in_edge_values(
    const Graph& /*unused*/,
    typename graph_traits<Graph>::vertex_descriptor /*unused*/,
    std::size_t /*unused*/, const int* /*unused*/, const int* /*unused*/) {}

template <typename Graph>
typename enable_if<is_incidence_graph<Graph>, void>::type
check_graph_out_degree(const Graph& g,
                       typename graph_traits<Graph>::vertex_descriptor u,
                       std::size_t expected_count) {
  BOOST_CHECK_EQUAL(out_degree(u, g), expected_count);
}

template <typename Graph>
typename disable_if<is_incidence_graph<Graph>, void>::type
check_graph_out_degree(
    const Graph& /*unused*/,
    typename graph_traits<Graph>::vertex_descriptor /*unused*/,
    std::size_t /*unused*/) {}

template <typename Graph>
typename enable_if<is_incidence_graph<Graph>, void>::type
check_graph_out_edge_values(const Graph& g,
                            typename graph_traits<Graph>::vertex_descriptor u,
                            std::size_t expected_count, const int* edge_values,
                            const int* vertex_values) {
  using OutEdgeIter = typename graph_traits<Graph>::out_edge_iterator;
  BOOST_CHECK_EQUAL(out_degree(u, g), expected_count);
  OutEdgeIter ei;
  OutEdgeIter ei_end;
  BOOST_CHECK_NO_THROW(tie(ei, ei_end) = out_edges(u, g));

  std::vector<int> e_list;
  std::vector<int> vp_list;
  for (; ei != ei_end; ++ei) {
    if (is_undirected_graph<Graph>::type::value) {
      BOOST_CHECK((g[*ei] == g[source(*ei, g)] * 1000 + g[target(*ei, g)]) ||
                  (g[*ei] == g[source(*ei, g)] + g[target(*ei, g)] * 1000));
    } else {
      BOOST_CHECK_EQUAL(g[*ei], (g[source(*ei, g)] * 1000 + g[target(*ei, g)]));
    }
    e_list.push_back(g[*ei]);
    vp_list.push_back(g[target(*ei, g)]);
  }
  std::sort(e_list.begin(), e_list.end());
  for (int& i : e_list) {
    BOOST_CHECK_EQUAL(i, *(edge_values++));
  }
  std::sort(vp_list.begin(), vp_list.end());
  for (int& i : vp_list) {
    BOOST_CHECK_EQUAL(i, *(vertex_values++));
  }
}

template <typename Graph>
typename disable_if<is_incidence_graph<Graph>, void>::type
check_graph_out_edge_values(
    const Graph& /*unused*/,
    typename graph_traits<Graph>::vertex_descriptor /*unused*/,
    std::size_t /*unused*/, const int* /*unused*/, const int* /*unused*/) {}

template <typename Graph>
typename disable_if<is_undirected_graph<Graph>, void>::type
intint_check_fullbranch_integrity(
    const Graph& g, typename graph_traits<Graph>::vertex_descriptor u) {
  using OutEdgeIter = typename graph_traits<Graph>::out_edge_iterator;

  if (out_degree(u, g) == 0) {
    return;
  }

  BOOST_CHECK_EQUAL(out_degree(u, g), 4);
  OutEdgeIter ei;
  OutEdgeIter ei_end;
  for (tie(ei, ei_end) = out_edges(u, g); ei != ei_end; ++ei) {
    BOOST_CHECK_EQUAL(g[*ei], g[source(*ei, g)] * 1000 + g[target(*ei, g)]);
    intint_check_fullbranch_integrity(g, target(*ei, g));
  }
}

template <typename Graph>
typename enable_if<is_undirected_graph<Graph>, void>::type
intint_check_fullbranch_integrity(
    const Graph& g,
    typename graph_traits<Graph>::vertex_descriptor /*unused*/) {
  using OutEdgeIter = typename graph_traits<Graph>::out_edge_iterator;
  using VertexIter = typename graph_traits<Graph>::vertex_iterator;

  VertexIter vi;
  VertexIter vi_end;
  BOOST_CHECK_NO_THROW(tie(vi, vi_end) = vertices(g));
  for (; vi != vi_end; ++vi) {
    OutEdgeIter ei;
    OutEdgeIter ei_end;
    for (tie(ei, ei_end) = out_edges(*vi, g); ei != ei_end; ++ei) {
      BOOST_CHECK((g[*ei] == g[source(*ei, g)] * 1000 + g[target(*ei, g)]) ||
                  (g[*ei] == g[source(*ei, g)] + g[target(*ei, g)] * 1000));
    }
  }
}

BOOST_AUTO_TEST_CASE_TEMPLATE(intint_bgl_mutable_graph_test, Graph,
                              intint_graphtest_types) {
  using Vertex = typename graph_traits<Graph>::vertex_descriptor;
  using Edge = typename graph_traits<Graph>::edge_descriptor;

  Graph g;
  Vertex v_root = Vertex();

  /* MutableGraph */
  BOOST_CHECK_NO_THROW(v_root = add_vertex(g));
  check_graph_vertex_count(g, 1);
  BOOST_CHECK_NO_THROW(remove_vertex(v_root, g));
  check_graph_vertex_count(g, 0);

  /* MutablePropertyGraph */
  BOOST_CHECK_NO_THROW(v_root = add_vertex(1, g));
  BOOST_CHECK_EQUAL(g[v_root], 1);
  g[v_root] = 1;
  BOOST_CHECK_EQUAL(g[v_root], 1);

  /* MutableGraph */
  std::array<int, 4> vp_rc = {2, 3, 4, 5};
  std::array<int, 4> ep_rc = {1002, 1003, 1004, 1005};
  std::array<Vertex, 4> v_rc;
  std::array<Edge, 4> e_rc;
  for (int i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW(v_rc[i] = add_vertex(g));
    g[v_rc[i]] = vp_rc[i];
    BOOST_CHECK_EQUAL(g[v_rc[i]], vp_rc[i]);
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW(tie(e_rc[i], edge_added_success) =
                             add_edge(v_root, v_rc[i], g));
    BOOST_CHECK(edge_added_success);
    g[e_rc[i]] = ep_rc[i];
    BOOST_CHECK_EQUAL(g[e_rc[i]], ep_rc[i]);
  }
  check_graph_vertex_count(g, 5);

  /* MutablePropertyGraph */
  std::array<int, 4> vp_rc1c = {6, 7, 8, 9};
  std::array<int, 4> ep_rc1c = {2006, 2007, 2008, 2009};
  std::array<Vertex, 4> v_rc1c;
  std::array<Edge, 4> e_rc1c;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW(v_rc1c[i] = add_vertex(vp_rc1c[i], g));
    BOOST_CHECK_EQUAL(g[v_rc1c[i]], vp_rc1c[i]);
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW(tie(e_rc1c[i], edge_added_success) =
                             add_edge(v_rc[0], v_rc1c[i], ep_rc1c[i], g));
    BOOST_CHECK(edge_added_success);
    BOOST_CHECK_EQUAL(g[e_rc1c[i]], ep_rc1c[i]);
  }
  check_graph_vertex_count(g, 9);

  BOOST_CHECK_EQUAL(g[v_root], 1);
  {
    /* IncidenceGraph */
    {
      const std::array<int, 4> e_vals = {1002, 1003, 1004, 1005};
      const std::array<int, 4> v_vals = {2, 3, 4, 5};
      check_graph_out_edge_values(g, v_root, 4, e_vals.data(), v_vals.data());
    }

    /* BidirectionalGraph */
    if (is_undirected_graph<Graph>::type::value) {
      const std::array<int, 5> e_vals = {1002, 2006, 2007, 2008, 2009};
      const std::array<int, 5> v_vals = {1, 6, 7, 8, 9};
      check_graph_in_edge_values(g, v_rc[0], 5, e_vals.data(), v_vals.data());
    } else if (is_bidirectional_graph<Graph>::type::value) {
      const std::array<int, 1> e_vals = {1002};
      const std::array<int, 1> v_vals = {1};
      check_graph_in_edge_values(g, v_rc[0], 1, e_vals.data(), v_vals.data());
    }

    /* IncidenceGraph */
    if (is_undirected_graph<Graph>::type::value) {
      const std::array<int, 5> e_vals = {1002, 2006, 2007, 2008, 2009};
      const std::array<int, 5> v_vals = {1, 6, 7, 8, 9};
      check_graph_out_edge_values(g, v_rc[0], 5, e_vals.data(), v_vals.data());
    } else {
      const std::array<int, 4> e_vals = {2006, 2007, 2008, 2009};
      const std::array<int, 4> v_vals = {6, 7, 8, 9};
      check_graph_out_edge_values(g, v_rc[0], 4, e_vals.data(), v_vals.data());
    }
  }

  /* MutablePropertyGraph (with rvalue-ref) */
  std::array<int, 4> vp_rc2c = {10, 11, 12, 13};
  std::array<int, 4> ep_rc2c = {3010, 3011, 3012, 3013};
  std::array<Vertex, 4> v_rc2c;
  std::array<Edge, 4> e_rc2c;
  for (std::size_t i = 0; i < 4; ++i) {
    BOOST_CHECK_NO_THROW(v_rc2c[i] = add_vertex(std::move(vp_rc2c[i]), g));
    bool edge_added_success = false;
    BOOST_CHECK_NO_THROW(tie(e_rc2c[i], edge_added_success) =
                             add_edge(v_rc[1], v_rc2c[i], ep_rc2c[i], g));
    BOOST_CHECK(edge_added_success);
  }
  check_graph_vertex_count(g, 13);

  if (is_undirected_graph<Graph>::type::value) {
    const std::array<int, 5> e_vals = {1003, 3010, 3011, 3012, 3013};
    const std::array<int, 5> v_vals = {1, 10, 11, 12, 13};
    check_graph_out_edge_values(g, v_rc[1], 5, e_vals.data(), v_vals.data());
  } else {
    /* IncidenceGraph */
    const std::array<int, 4> e_vals = {3010, 3011, 3012, 3013};
    const std::array<int, 4> v_vals = {10, 11, 12, 13};
    check_graph_out_edge_values(g, v_rc[1], 4, e_vals.data(), v_vals.data());
  }

  /* Copying function */
  intint_check_fullbranch_integrity(g, v_root);
  {
    Graph* p_g_cpy = nullptr;
    BOOST_CHECK_NO_THROW(p_g_cpy = new Graph(g));
    intint_check_fullbranch_integrity(*p_g_cpy, *(vertices(*p_g_cpy).first));
    BOOST_CHECK_NO_THROW(delete p_g_cpy);
  }

  {
    Graph g_cpy;
    BOOST_CHECK_NO_THROW(g_cpy = g);
    intint_check_fullbranch_integrity(g_cpy, *(vertices(g_cpy).first));
  }

  {
    Graph* p_g_mv = nullptr;
    BOOST_CHECK_NO_THROW(p_g_mv = new Graph(std::move(g)));
    intint_check_fullbranch_integrity(*p_g_mv, *(vertices(*p_g_mv).first));
    BOOST_CHECK_NO_THROW(g = std::move(*p_g_mv));
    intint_check_fullbranch_integrity(g, *(vertices(g).first));
    BOOST_CHECK_NO_THROW(delete p_g_mv);
    v_root = *(vertices(g).first);
  }

  /* MutableGraph */
  BOOST_CHECK_NO_THROW(clear_vertex(v_rc[0], g));

  if (is_undirected_graph<Graph>::type::value) {
    /* IncidenceGraph */
    check_graph_out_degree(g, v_rc[0], 0);
    check_graph_out_degree(g, v_root, 3);

    /* BidirectionalGraph */
    check_graph_in_degree(g, v_rc[0], 0);
    check_graph_in_degree(g, v_rc1c[0], 0);
    check_graph_in_degree(g, v_rc1c[1], 0);
    check_graph_in_degree(g, v_rc1c[2], 0);
    check_graph_in_degree(g, v_rc1c[3], 0);
  } else {
    /* IncidenceGraph */
    check_graph_out_degree(g, v_rc[0], 0);
    check_graph_out_degree(g, v_root, 3);

    /* BidirectionalGraph */
    check_graph_in_degree(g, v_rc[0], 0);
    check_graph_in_degree(g, v_rc1c[0], 0);
    check_graph_in_degree(g, v_rc1c[1], 0);
    check_graph_in_degree(g, v_rc1c[2], 0);
    check_graph_in_degree(g, v_rc1c[3], 0);
  }

  /* VertexListGraph */
  check_graph_vertex_count(g, 13);
  {
    const std::array<int, 13> ref_values = {1, 2, 3,  4,  5,  6, 7,
                                            8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values.data());
  }

  /* EdgeListGraph */
  check_graph_edge_count(g, 7);
  {
    const std::array<int, 7> ref_values = {1003, 1004, 1005, 3010,
                                           3011, 3012, 3013};
    intint_check_graph_edge_values(g, 7, ref_values.data());
  }

  /* MutableGraph */
  BOOST_CHECK_NO_THROW(remove_edge(v_rc[1], v_rc2c[2], g));

  /* VertexListGraph */
  check_graph_vertex_count(g, 13);
  {
    const std::array<int, 13> ref_values = {1, 2, 3,  4,  5,  6, 7,
                                            8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values.data());
  }

  /* EdgeListGraph */
  check_graph_edge_count(g, 6);
  {
    const std::array<int, 6> ref_values = {1003, 1004, 1005, 3010, 3011, 3013};
    intint_check_graph_edge_values(g, 6, ref_values.data());
  }

  /* MutableGraph */
  std::pair<Edge, bool> last_e_of_rc2;
  BOOST_CHECK_NO_THROW(last_e_of_rc2 = edge(v_rc[1], v_rc2c[3], g));
  BOOST_CHECK_EQUAL(last_e_of_rc2.second, true);
  BOOST_CHECK_NO_THROW(remove_edge(last_e_of_rc2.first, g));

  /* VertexListGraph */
  check_graph_vertex_count(g, 13);
  {
    const std::array<int, 13> ref_values = {1, 2, 3,  4,  5,  6, 7,
                                            8, 9, 10, 11, 12, 13};
    intint_check_graph_vertex_values(g, 13, ref_values.data());
  }

  /* EdgeListGraph */
  check_graph_edge_count(g, 5);
  {
    const std::array<int, 5> ref_values = {1003, 1004, 1005, 3010, 3011};
    intint_check_graph_edge_values(g, 5, ref_values.data());
  }

  /* MutableGraph */
  BOOST_CHECK_NO_THROW(clear_vertex(v_rc2c[0], g));
  BOOST_CHECK_NO_THROW(remove_vertex(v_rc2c[0], g));

  /* VertexListGraph */
  check_graph_vertex_count(g, 12);
  {
    const std::array<int, 12> ref_values = {1, 2, 3, 4,  5,  6,
                                            7, 8, 9, 11, 12, 13};
    intint_check_graph_vertex_values(g, 12, ref_values.data());
  }

  /* EdgeListGraph */
  check_graph_edge_count(g, 4);
  {
    const std::array<int, 4> ref_values = {1003, 1004, 1005, 3011};
    intint_check_graph_edge_values(g, 4, ref_values.data());
  }
}
