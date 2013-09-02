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

#define BOOST_TEST_MODULE boost_graph_trees
#include <boost/test/unit_test.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/front_inserter.hpp>


using namespace boost;

