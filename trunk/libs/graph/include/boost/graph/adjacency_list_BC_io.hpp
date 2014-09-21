//=======================================================================
// Copyright 2014 Sven Mikael Persson
// Author: Sven Mikael Persson
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================

#ifndef BOOST_GRAPH_ADJACENCY_LIST_BC_IO_HPP
#define BOOST_GRAPH_ADJACENCY_LIST_BC_IO_HPP

#include <boost/graph/adjacency_list_io.hpp>

// Method read to parse an adjacency list from an input stream. Examples:
// cin >> read( G );
// cin >> read( G, NodePropertySubset(), EdgepropertySubset() );
//
// Method write to print an adjacency list to an output stream. Examples:
// cout << write( G );
// cout << write( G, NodePropertySubset(), EdgepropertySubset() );

namespace boost {

// graph printer
//=========================================================================
// user methods

/// graph parser for given subsets of internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP, class VPS, class EPS>
GraphParser<adjacency_list_BC<EL,VL,D,VP,EP>,VP,EP,VPS,EPS> 
read( adjacency_list_BC<EL,VL,D,VP,EP,GP>& g, VPS vps, EPS eps )
{
  return GraphParser<adjacency_list_BC<EL,VL,D,VP,EP,GP>,VP,EP,VPS,EPS>(&g);
}

/// graph parser for all internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP>
GraphParser<adjacency_list_BC<EL,VL,D,VP,EP>,VP,EP,VP,EP> 
read( adjacency_list_BC<EL,VL,D,VP,EP>& g )
{
  return GraphParser<adjacency_list_BC<EL,VL,D,VP,EP>,VP,EP,VP,EP>(&g);
}


/// write the graph with given property subsets
template<class EL, class VL, class D, class VP, class EP, class VPS, class EPS>
GraphPrinter<adjacency_list_BC<EL,VL,D,VP,EP>,VPS,EPS> 
write( const adjacency_list_BC<EL,VL,D,VP,EP>& g, VPS, EPS )
{
  return GraphPrinter<adjacency_list_BC<EL,VL,D,VP,EP>,VPS,EPS>(g);
}

/// write the graph with all internal vertex and edge properties
template<class EL, class VL, class D, class VP, class EP>
GraphPrinter<adjacency_list_BC<EL,VL,D,VP,EP>,VP,EP> 
write( const adjacency_list_BC<EL,VL,D,VP,EP>& g )
{
  return GraphPrinter<adjacency_list_BC<EL,VL,D,VP,EP>,VP,EP>(g);
}

// user methods
//=========================================================================
}// boost
#endif
