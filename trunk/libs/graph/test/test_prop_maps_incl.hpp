// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)


#include <iostream>

#include <boost/graph/linked_tree_BC.hpp>
#include <boost/graph/adjacency_list_BC.hpp>

#include <boost/graph/tree_adaptor.hpp>


#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE boost_graph_prop_maps
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/front_inserter.hpp>


using namespace boost;



typedef property< vertex_name_t, std::string, 
        property< vertex_distance_t, int > > VertexPropTest;

typedef property< edge_name_t, std::string, 
        property< edge_weight_t, int > > EdgePropTest;

typedef property< graph_name_t, std::string > GraphPropTest;

struct PropMapMaps {
  typedef vertex_name_t vname_t;
  typedef vertex_distance_t vdistance_t;
  typedef edge_name_t ename_t;
  typedef edge_weight_t eweight_t;
  typedef graph_name_t gname_t;
  
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
  typedef std::string VertexBundleTest::* vname_t;
  typedef int VertexBundleTest::* vdistance_t;
  typedef std::string EdgeBundleTest::* ename_t;
  typedef int EdgeBundleTest::* eweight_t;
  typedef std::string GraphBundleTest::* gname_t;
  
  static const vname_t vname;
  static const vdistance_t vdistance;
  static const ename_t ename;
  static const eweight_t eweight;
  static const gname_t gname;
};

const BundleMaps::vname_t BundleMaps::vname = &VertexBundleTest::name;
const BundleMaps::vdistance_t BundleMaps::vdistance = &VertexBundleTest::distance;
const BundleMaps::ename_t BundleMaps::ename = &EdgeBundleTest::name;
const BundleMaps::eweight_t BundleMaps::eweight = &EdgeBundleTest::weight;
const BundleMaps::gname_t BundleMaps::gname = &GraphBundleTest::name;





