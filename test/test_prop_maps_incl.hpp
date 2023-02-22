// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>

#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/graph/linked_tree_BC.hpp>

#include <boost/graph/tree_adaptor.hpp>

#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE boost_graph_prop_maps
#include <boost/mpl/copy.hpp>
#include <boost/mpl/front_inserter.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost;

using VertexPropTest =
    property<vertex_name_t, std::string, property<vertex_distance_t, int>>;

using EdgePropTest =
    property<edge_name_t, std::string, property<edge_weight_t, int>>;

using GraphPropTest = property<graph_name_t, std::string>;

struct PropMapMaps {
  using vname_t = vertex_name_t;
  using vdistance_t = vertex_distance_t;
  using ename_t = edge_name_t;
  using eweight_t = edge_weight_t;
  using gname_t = graph_name_t;

  static const vname_t vname = vertex_name;
  static const vdistance_t vdistance = vertex_distance;
  static const ename_t ename = edge_name;
  static const eweight_t eweight = edge_weight;
  static const gname_t gname = graph_name;
};

struct VertexBundleTest {
  std::string name;
  int distance;
};

struct EdgeBundleTest {
  std::string name;
  int weight;
};

struct GraphBundleTest {
  std::string name;
};

struct BundleMaps {
  using vname_t = std::string VertexBundleTest::*;
  using vdistance_t = int VertexBundleTest::*;
  using ename_t = std::string EdgeBundleTest::*;
  using eweight_t = int EdgeBundleTest::*;
  using gname_t = std::string GraphBundleTest::*;

  static const vname_t vname;
  static const vdistance_t vdistance;
  static const ename_t ename;
  static const eweight_t eweight;
  static const gname_t gname;
};

const BundleMaps::vname_t BundleMaps::vname = &VertexBundleTest::name;
const BundleMaps::vdistance_t BundleMaps::vdistance =
    &VertexBundleTest::distance;
const BundleMaps::ename_t BundleMaps::ename = &EdgeBundleTest::name;
const BundleMaps::eweight_t BundleMaps::eweight = &EdgeBundleTest::weight;
const BundleMaps::gname_t BundleMaps::gname = &GraphBundleTest::name;
