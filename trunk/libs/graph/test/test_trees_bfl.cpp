// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/graph/bfl_d_ary_tree.hpp>

#include "test_trees_incl.hpp"



template <typename VertexProp, typename EdgeProp>
struct graphtype_list {
  
  typedef mpl::list< 
    bfl_d_ary_tree< 4, VertexProp, EdgeProp> > types;
  
};

typedef graphtype_list<int, int>::types intint_treetest_types;

#include "test_trees_impl.hpp"






