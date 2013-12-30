// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file tree_concepts.hpp
 * 
 * This library provides concept classes to verify that a data structure has a tree interface, 
 * i.e., models the concepts of trees as used in the Boost Graph Library. The tree concepts 
 * should be regarded as a special kind of graph and reuses much of the general graph traits.
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */

#ifndef BOOST_TREE_CONCEPTS_HPP
#define BOOST_TREE_CONCEPTS_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>
#include <boost/concept_check.hpp>

#include <boost/graph/tree_traits.hpp>

#include <utility>

namespace boost {

/**
 * This concept defines the requirements to fulfill in order to model a tree 
 * as used in the Boost Graph Library.
 * 
 * Required Concepts:
 * 
 * TreeType should model the IncidenceGraphConcept.
 * 
 * Valid expressions (vertex v; child_vertex_iterator cv_it, cv_it_end; TreeType tree):
 * 
 * v = get_root_vertex(tree);  The root vertex of the tree can be obtained.
 * 
 * tie(cv_it, cv_it_end) = child_vertices(v,tree);  The iterator range for child vertices of a given parent vertex (v) in a tree can be obtained.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct TreeConcept {
  typename graph_traits<TreeType>::vertex_descriptor v;
  typename tree_traits<TreeType>::child_vertex_iterator cv_it, cv_it_end;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((IncidenceGraphConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(TreeConcept) 
  {
    v = get_root_vertex(tree);
    tie(cv_it, cv_it_end) = child_vertices(v,tree);
  };
  
};

/**
 * This concept defines the requirements to fulfill in order to model a bidirectional tree 
 * as used in the Boost Graph Library.
 * 
 * Required Concepts:
 * 
 * TreeType should model the TreeConcept.
 * 
 * Valid expressions (vertex v; TreeType tree):
 * 
 * u = parent_vertex(v, tree);  The parent vertex of a given vertex can be obtained via a call to "parent_vertex" free-function.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct BidirectionalTreeConcept {
  typename graph_traits<TreeType>::vertex_descriptor v;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((TreeConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(BidirectionalTreeConcept) 
  {
    v = parent_vertex(v, tree);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a mutable tree 
 * as used in the Boost Graph Library.
 * 
 * Required concepts:
 * 
 * TreeType should model the TreeConcept.
 * 
 * Valid expressions (vertex u, v; edge e; TreeType tree):
 * 
 * v = create_root(tree);  The root vertex of the tree can be created.
 * 
 * tie(v, e) = add_child_vertex(u, tree);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) in a tree.
 * 
 * remove_branch(v, tree);  The entire branch (or sub-tree) below a given vertex (including the vertex itself) can be removed from the tree.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct MutableTreeConcept {
  typename graph_traits<TreeType>::vertex_descriptor u, v;
  typename graph_traits<TreeType>::edge_descriptor e;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((TreeConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(MutableTreeConcept) 
  {
    v = create_root(tree);
    tie(v, e) = add_child_vertex(u, tree);
    remove_branch(v, tree);
  };
  
};


/**
 * This concept defines the requirements to fulfill in order to model a mutable property-tree 
 * as used in the Boost Graph Library. A mutable property-tree is essentially a mutable tree 
 * whose mutating functions take or deliver the vertex- or edge- property values associated with 
 * the vertices or edges. During removal of a branch, all the vertex-properties are collected into 
 * an output iterator (e.g., back-inserter). During additions of child nodes, the corresponding 
 * vertex-properties can be used to initialize the new vertex and edge directly. This not only 
 * makes such a tree easier to use (not having to manually collect or set vertex properties before 
 * or after the mutation), but it can also be necessary in some situations. One typical use-case 
 * is when re-balancing a branch of a tree, which often results in collecting properties 
 * of vertices (to preserve them), clearing the branch, re-adding the vertices in a new arrangement, and 
 * restoring their original properties.
 * 
 * Required concepts:
 * 
 * TreeType should model the MutableTreeConcept.
 * 
 * Valid expressions (vertex u, v; edge e; TreeType tree; vertex_property_type vp; edge_property_type ep; ):
 * 
 * v = create_root(vp, g);  The root vertex of the tree can be created with a given vertex-property value (vp).
 * 
 * tie(v,e) = add_child_vertex(u, vp, g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) and with a given vertex-property value (vp) in a tree.
 * 
 * tie(v,e) = add_child_vertex(u, vp, ep, g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u), with a given vertex-property value (vp) and  with a given edge-property value (ep) in a tree.
 * 
 * remove_branch(v, back_inserter(vp_vect), g);  The entire branch (or sub-tree) below a given vertex (including the vertex itself) can be removed from the tree, with its vertex-property values put on an output iterator (e.g., a back-inserter).
 * 
 * C++11 only:
 * 
 * v = create_root(std::move(vp), g);  The root vertex of the tree can be created with a given vertex-property value (vp) (to move in).
 *
 * tie(v,e) = add_child_vertex(u, std::move(vp), g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) and with a given vertex-property value (vp) (to move in) in a tree.
 *
 * tie(v,e) = add_child_vertex(u, std::move(vp), std::move(ep), g);  A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u), with a given vertex-property value (vp) (to move in) and  with a given edge-property value (ep) (to move in) in a tree.
 * 
 * \tparam TreeType The tree type to be checked for this concept.
 */
template <typename TreeType>
struct MutablePropertyTreeConcept {
  typename graph_traits<TreeType>::vertex_descriptor u, v;
  typename graph_traits<TreeType>::edge_descriptor e;
  typename TreeType::vertex_property_type vp;
  typename TreeType::edge_property_type ep;
  TreeType tree;
  
  BOOST_CONCEPT_ASSERT((TreeConcept<TreeType>));
  
  BOOST_CONCEPT_USAGE(MutablePropertyTreeConcept) 
  {
    v = create_root(vp, tree);
    tie(v,e) = add_child_vertex(u, vp, tree);
    tie(v,e) = add_child_vertex(u, vp, ep, tree);
    std::vector< typename TreeType::vertex_property_type > vp_vect;
    remove_branch(v, back_inserter(vp_vect), tree);
    
#ifndef BOOST_NO_CXX11_RVALUE_REFERENCES
    v = create_root(std::move(vp), tree);
    tie(v,e) = add_child_vertex(u, std::move(vp), tree);
    tie(v,e) = add_child_vertex(u, std::move(vp), std::move(ep), tree);
#endif
  };
  
};


};


#endif


















