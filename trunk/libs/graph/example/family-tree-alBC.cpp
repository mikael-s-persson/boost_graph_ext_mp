//=======================================================================
// Copyright 2014 Sven Mikael Persson, 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/tuple/tuple.hpp>

int main()
{
  using namespace boost;
  
  typedef adjacency_list_BC<vecBC, vecBC, directedS, std::string> Graph;
  typedef graph_traits< Graph >::vertex_descriptor Vertex;
  Graph family;
  
  Vertex Jeanie   = add_vertex("Jeanie",   family);
  Vertex Debbie   = add_vertex("Debbie",   family);
  Vertex Rick     = add_vertex("Rick",     family);
  Vertex John     = add_vertex("John",     family);
  Vertex Amanda   = add_vertex("Amanda",   family);
  Vertex Margaret = add_vertex("Margaret", family);
  Vertex Benjamin = add_vertex("Benjamin", family);
  
  add_edge(Jeanie, Debbie, family);
  add_edge(Jeanie, Rick,   family);
  add_edge(Jeanie, John,   family);
  add_edge(Debbie, Amanda, family);
  add_edge(Rick, Margaret, family);
  add_edge(John, Benjamin, family);
  
  graph_traits< Graph >::vertex_iterator i, end;
  graph_traits< Graph >::adjacency_iterator ai, a_end;

  for (tie(i, end) = vertices(family); i != end; ++i) {
    std::cout << family[*i];
    tie(ai, a_end) = adjacent_vertices(*i, family);
    if (ai == a_end)
      std::cout << " has no children";
    else
      std::cout << " is the parent of ";
    for (; ai != a_end; ++ai) {
      std::cout << family[*ai];
      if (next(ai) != a_end)
        std::cout << ", ";
    }
    std::cout << std::endl;
  }
  return EXIT_SUCCESS;
}
