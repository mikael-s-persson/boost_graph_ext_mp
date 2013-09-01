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
#include <boost/concept_check.hpp>

namespace boost {


  namespace detail {
    
    template <typename A> struct tree_traits_try_child_iter {typedef void type;};
    
    template <typename Graph, typename Enable = void>
    struct tree_traits_get_child_iter {
      typedef typename Graph::adjacency_iterator type;
    };
    
    template <typename Graph>
    struct tree_traits_get_child_iter<Graph, 
             typename tree_traits_try_child_iter<typename Graph::child_vertex_iterator>::type> {
      typedef typename Graph::child_vertex_iterator type;
    };
    
  };

/**
 * This traits class defines a number of nested types associated to a tree structure.
 */
template <typename TreeType>
struct tree_traits {
  /** This type describes iterators to iterate through child vertices of a vertex. */
  typedef typename detail::tree_traits_get_child_iter<TreeType>::type child_vertex_iterator;
};



};

#endif
