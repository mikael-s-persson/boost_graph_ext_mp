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
#include <boost/graph/vebl_d_ary_tree.hpp>
#include <boost/tuple/tuple.hpp>

int main()
{
  using namespace boost;
  
  typedef vebl_d_ary_tree<3, std::string> Graph;
  typedef graph_traits< Graph >::vertex_descriptor Vertex;
  Graph family;
  
  Vertex Jeanie = create_root("Jeanie", family);
  
  Vertex Debbie = add_child_vertex(Jeanie, "Debbie", family).first;
  Vertex Rick   = add_child_vertex(Jeanie, "Rick", family).first;
  Vertex John   = add_child_vertex(Jeanie, "John", family).first;
  Vertex Amanda   = add_child_vertex(Debbie, "Amanda", family).first;
  Vertex Margaret = add_child_vertex(Rick, "Margaret", family).first;
  Vertex Benjamin = add_child_vertex(John, "Benjamin", family).first;
  
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

