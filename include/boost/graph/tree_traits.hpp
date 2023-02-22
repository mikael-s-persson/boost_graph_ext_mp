// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file tree_traits.hpp
 * 
 * This library provides traits for tree types, i.e., models the concepts of trees as used 
 * in the Boost Graph Library. The tree traits are an add-on to graph traits.
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */

#ifndef BOOST_TREE_TRAITS_HPP
#define BOOST_TREE_TRAITS_HPP

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/properties.hpp>

#include <boost/mpl/bool.hpp>
#include <boost/mpl/if.hpp>

namespace boost {

namespace detail {

template <typename A>
struct tree_traits_try_child_iter {
  using type = void;
};

template <typename Graph, typename Enable = void>
struct tree_traits_get_child_iter {
  using type = typename Graph::adjacency_iterator;
};

template <typename Graph>
struct tree_traits_get_child_iter<
    Graph, typename tree_traits_try_child_iter<
               typename Graph::child_vertex_iterator>::type> {
  using type = typename Graph::child_vertex_iterator;
};

}  // namespace detail

/**
 * This traits class defines a number of nested types associated to a tree structure.
 */
template <typename TreeType>
struct tree_traits {
  /** This type describes iterators to iterate through child vertices of a vertex. */
  using child_vertex_iterator =
      typename detail::tree_traits_get_child_iter<TreeType>::type;
};

/**
 * 
 */
template <typename VertexDescriptor, typename EdgeDescriptor,
          typename StorageTag>
struct tree_storage {
  using type = typename StorageTag::template bind<VertexDescriptor,
                                                  EdgeDescriptor>::type;
};

/**
 * 
 */
template <typename StorageTag>
struct tree_storage_traits {
  using is_rand_access = boost::mpl::false_;
  using is_bidir = boost::mpl::true_;
  using is_directed = boost::mpl::true_;

  using directed_category = typename boost::mpl::if_<
      is_bidir, boost::bidirectional_tag,
      typename boost::mpl::if_<is_directed, boost::directed_tag,
                               boost::undirected_tag>::type>::type;

  using edge_parallel_category = boost::disallow_parallel_edge_tag;

  using vertices_size_type = std::size_t;
  using vertex_descriptor = typename StorageTag::vertex_descriptor;
  using edges_size_type = std::size_t;
  using edge_descriptor = typename StorageTag::edge_descriptor;
};

}  // namespace boost

#endif
