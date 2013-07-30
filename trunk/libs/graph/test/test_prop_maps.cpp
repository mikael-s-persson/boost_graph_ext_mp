
/*
 *    Copyright 2013 Sven Mikael Persson
 *
 *    THIS SOFTWARE IS DISTRIBUTED UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE v3 (GPLv3).
 *
 *    This file is part of ReaK.
 *
 *    ReaK is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    ReaK is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with ReaK (as LICENSE in the root folder).  
 *    If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/tree_adaptor.hpp>
#include <boost/graph/linked_tree.hpp>
#include <boost/graph/pooled_adjacency_list.hpp>

// #include <boost/graph/d_ary_bf_tree.hpp>
// #include <boost/graph/d_ary_cob_tree.hpp>



#define BOOST_TEST_DYN_LINK

#define BOOST_TEST_MODULE boost_graph_prop_maps
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>


typedef boost::property< boost::vertex_name_t, std::string, 
        boost::property< boost::vertex_distance_t, int > > VertexPropTest;

typedef boost::property< boost::edge_name_t, std::string, 
        boost::property< boost::edge_weight_t, int > > EdgePropTest;

typedef boost::property< boost::graph_name_t, std::string > GraphPropTest;


typedef boost::mpl::list< 
  boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::listS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::listS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::listS, boost::listS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::listS, boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::listS, boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::vecS, boost::vecS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::vecS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::vecS, boost::listS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::adjacency_list< boost::listS, boost::listS, boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS>,
  boost::linked_tree<boost::vecS, boost::vecS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::linked_tree<boost::listS, boost::vecS, VertexPropTest, EdgePropTest, GraphPropTest>, 
  boost::linked_tree<boost::vecS, boost::listS, VertexPropTest, EdgePropTest, GraphPropTest>, 
  boost::linked_tree<boost::listS, boost::listS, VertexPropTest, EdgePropTest, GraphPropTest>,
  boost::pooled_adjacency_list<boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest >,
  boost::pooled_adjacency_list<boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest >,
  boost::pooled_adjacency_list<boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest >,
  boost::pooled_adjacency_list<boost::bidirectionalS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS >,
  boost::pooled_adjacency_list<boost::undirectedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS >,
  boost::pooled_adjacency_list<boost::directedS, VertexPropTest, EdgePropTest, GraphPropTest, boost::vecS > > intint_propmaptest_types;
  


BOOST_AUTO_TEST_CASE_TEMPLATE( intint_propmaptest, Graph, intint_propmaptest_types )
{
  
};




