//=======================================================================
// Copyright 2014 Sven Mikael Persson
// Authors: Sven Mikael Persson
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#ifndef BGL_DETAIL_ANY_TREE_SERIALIZE_HPP
#define BGL_DETAIL_ANY_TREE_SERIALIZE_HPP

#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/graph/tree_traits.hpp>
#include <boost/pending/property_serialize.hpp>

#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/split_free.hpp>

#include <map>
#include <queue>
#include <vector>

namespace boost {

namespace graph {

namespace detail {

template <typename Archive, typename Graph>
inline void save_tree(Archive& ar, const Graph& graph,
                      const unsigned int /* file_version */
) {
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using Edge = typename boost::graph_traits<Graph>::edge_descriptor;
  using Task = std::tuple<Vertex, int, Edge>;

  int V = num_vertices(graph);
  ar << BOOST_SERIALIZATION_NVP(V);

  if (V == 0)
    return;

  // assign indices to vertices
  std::map<Vertex, int> indices;
  std::queue<Task> to_visit;
  int num = 0;

  Vertex u = get_root_vertex(graph);

  while (true) {
    indices[u] = num;
    ar << boost::serialization::make_nvp("vertex_property", graph[u]);
    // iterate through children nodes:
    typename boost::graph_traits<Graph>::out_edge_iterator oe_it, oe_it_end;
    std::tie(oe_it, oe_it_end) = out_edges(u, graph);
    while (oe_it != oe_it_end) {
      to_visit.push(Task(target(*oe_it, graph), num, *oe_it));
      ++oe_it;
    }
    // pop a node from the queue:
    if (to_visit.empty())
      return;
    Task t = to_visit.front();
    to_visit.pop();
    ar << boost::serialization::make_nvp("parent", get<1>(Task));
    ar << boost::serialization::make_nvp("edge_property", graph[get<2>(Task)]);
    u = get<0>(Task);
    ++num;
  }
}

template <typename Archive, typename Graph>
inline void load_tree(Archive& ar, Graph& graph,
                      const unsigned int /* file_version */
) {
  using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
  using VertexProp = typename Graph::vertex_property_type;
  using EdgeProp = typename Graph::edge_property_type;

  if (num_vertices(graph))
    graph = Graph();

  unsigned int V;
  ar >> BOOST_SERIALIZATION_NVP(V);

  if (V == 0)
    return;

  std::vector<Vertex> verts(V);
  int i = 0;

  VertexProp vp;
  EdgeProp ep;
  ar >> boost::serialization::make_nvp("vertex_property", vp);
  Vertex v = create_root(std::move(vp), graph);
  verts[i++] = v;

  while (--V > 0) {
    int parent = 0;
    ar >> boost::serialization::make_nvp("parent", parent);
    ar >> boost::serialization::make_nvp("edge_property", ep);
    ar >> boost::serialization::make_nvp("vertex_property", vp);
    v = add_child_vertex(verts[parent], std::move(vp), std::move(ep), tree)
            .first;
    verts[i++] = v;
  }
}

}  // namespace detail
}  // namespace graph
}  // namespace boost

#endif  // BGL_DETAIL_ANY_TREE_SERIALIZE_HPP
