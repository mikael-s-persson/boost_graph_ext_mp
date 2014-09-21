//=======================================================================
// Copyright 2014 Sven Mikael Persson
// Authors: Sven Mikael Persson 
// Adapted from original code by Jeremy G. Siek (2005)
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#ifndef ADJ_LIST_BC_SERIALIZE_HPP
#define ADJ_LIST_BC_SERIALIZE_HPP

#include <boost/graph/adjacency_list_BC.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/pending/property_serialize.hpp>
#include <boost/config.hpp>
#include <boost/detail/workaround.hpp>

#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/split_free.hpp>

namespace boost { 

namespace serialization {

// Turn off tracking for adjacency_list_BC. It's not polymorphic, and we
// need to do this to enable saving of non-const adjacency lists.
template <typename OEL, typename VL, typename D, typename VP, typename EP>
struct tracking_level<boost::adjacency_list_BC<OEL,VL,D,VP,EP> > {
  typedef mpl::integral_c_tag tag;
  typedef mpl::int_<track_never> type;
  BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);
};

template <typename Archive, typename OEL, typename VL, 
          typename D, typename VP, typename EP>
inline void save(
    Archive & ar,
    const boost::adjacency_list_BC<OEL,VL,D,VP,EP> &graph,
    const unsigned int /* file_version */
){
  typedef boost::adjacency_list_BC<OEL,VL,D,VP,EP> Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;

  int V = num_vertices(graph);
  int E = num_edges(graph);
  ar << BOOST_SERIALIZATION_NVP(V);
  ar << BOOST_SERIALIZATION_NVP(E);

  // assign indices to vertices
  std::map<Vertex,int> indices;
  int num = 0;
  BGL_FORALL_VERTICES_T(v, graph, Graph) {
    indices[v] = num++;
    ar << serialization::make_nvp("vertex_property", graph[v] );
//     ar << serialization::make_nvp("vertex_property", get(vertex_all_t(), graph, v) );
  }
  
  // write edges
  BGL_FORALL_EDGES_T(e, graph, Graph) {
    ar << serialization::make_nvp("u" , indices[source(e,graph)]);
    ar << serialization::make_nvp("v" , indices[target(e,graph)]);
    ar << serialization::make_nvp("edge_property", graph[e] );
//     ar << serialization::make_nvp("edge_property", get(edge_all_t(), graph, e) );
  }

//   ar << serialization::make_nvp("graph_property", get_property(graph, graph_all_t()) );
}


template <typename Archive, typename OEL, typename VL, 
          typename D, typename VP, typename EP>
inline void load(
    Archive & ar,
    boost::adjacency_list_BC<OEL,VL,D,VP,EP> &graph,
    const unsigned int /* file_version */
){
  typedef boost::adjacency_list_BC<OEL,VL,D,VP,EP,GP,EL> Graph;
  typedef typename boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef typename Graph::vertex_property_type VertexProp;
  typedef typename Graph::edge_property_type EdgeProp;
  
  if(num_vertices(graph))
    graph = Graph();

  unsigned int V;
  ar >> BOOST_SERIALIZATION_NVP(V);
  unsigned int E;
  ar >> BOOST_SERIALIZATION_NVP(E);
  
  std::vector<Vertex> verts(V);
  int i = 0;
  while(V-- > 0){
    VertexProp vp;
    ar >> serialization::make_nvp("vertex_property", vp );
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    Vertex v = add_vertex(std::move(vp), graph);
#else
    Vertex v = add_vertex(vp, graph);
#endif
    verts[i++] = v;
  }
  
  while(E-- > 0){
    int u; int v;
    ar >> BOOST_SERIALIZATION_NVP(u);
    ar >> BOOST_SERIALIZATION_NVP(v);
    EdgeProp ep;
    ar >> serialization::make_nvp("edge_property", ep );
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    add_edge(verts[u], verts[v], std::move(ep), graph);
#else
    add_edge(verts[u], verts[v], ep, graph);
#endif
  }
  
//   ar >> serialization::make_nvp("graph_property", get_property(graph, graph_all_t()) );
}

template <typename Archive, typename OEL, typename VL, typename D, typename VP, typename EP>
inline void serialize(
    Archive & ar,
    boost::adjacency_list_BC<OEL,VL,D,VP,EP> &graph,
    const unsigned int file_version
){
  boost::serialization::split_free(ar, graph, file_version);
}

}//serialization
}//boost


#endif // ADJ_LIST_BC_SERIALIZE_HPP
