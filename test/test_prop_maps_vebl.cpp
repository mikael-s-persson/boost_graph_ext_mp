// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/graph/vebl_d_ary_tree.hpp>

#include "test_prop_maps_incl.hpp"

template <typename VertexProp, typename EdgeProp>
struct graphtype_list {
  
  typedef mpl::list< 
    vebl_d_ary_tree< 4, VertexProp, EdgeProp> > types;
  
};

typedef graphtype_list<VertexPropTest, EdgePropTest>::types propmap_graphtest_types;
typedef graphtype_list<VertexBundleTest, EdgeBundleTest>::types bundlemap_graphtest_types;


#define PROPMAP_GRAPHTEST_NAME propmap_graphtest
#define PROPMAP_GRAPHTEST_TYPES propmap_graphtest_types
#define PROPMAP_GRAPHTEST_MAPS PropMapMaps

#include "test_prop_maps_impl.hpp"

#undef PROPMAP_GRAPHTEST_NAME
#undef PROPMAP_GRAPHTEST_TYPES
#undef PROPMAP_GRAPHTEST_MAPS



#define PROPMAP_GRAPHTEST_NAME bundlemap_graphtest
#define PROPMAP_GRAPHTEST_TYPES bundlemap_graphtest_types
#define PROPMAP_GRAPHTEST_MAPS BundleMaps

#include "test_prop_maps_impl.hpp"

#undef PROPMAP_GRAPHTEST_NAME
#undef PROPMAP_GRAPHTEST_TYPES
#undef PROPMAP_GRAPHTEST_MAPS











