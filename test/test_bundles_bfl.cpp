// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/graph/bfl_d_ary_tree.hpp>

#include "test_prop_maps_incl.hpp"

template <typename VertexProp, typename EdgeProp> struct graphtype_list {

  using types = ::testing::Types<bfl_d_ary_tree<4, VertexProp, EdgeProp>>;
};

using PropMapGraphTestTypes =
    graphtype_list<VertexBundleTest, EdgeBundleTest>::types;

template <typename T> class PropMapGraphTest : public ::testing::Test {};

#define PROPMAP_GRAPHTEST_MAPS BundleMaps

#include "test_prop_maps_impl.hpp"

#undef PROPMAP_GRAPHTEST_MAPS
