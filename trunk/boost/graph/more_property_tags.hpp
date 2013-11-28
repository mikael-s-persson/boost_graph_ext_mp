// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file more_property_tags.hpp
 *
 * This library adds a few property-tags (useful in a path-planning / shortest-path code).
 * 
 * \author Sven Mikael Persson
 * \date May 2013
 */

#ifndef BOOST_MORE_PROPERTY_TAGS_HPP
#define BOOST_MORE_PROPERTY_TAGS_HPP

#include <boost/graph/properties.hpp>

namespace boost {

  enum vertex_heuristic_t { vertex_heuristic };
  enum vertex_rhs_t { vertex_rhs };
  enum vertex_key_t { vertex_key };
  enum vertex_position_t { vertex_position };
  enum vertex_second_bundle_t { vertex_second_bundle };
  
  BOOST_INSTALL_PROPERTY(vertex, heuristic);
  BOOST_INSTALL_PROPERTY(vertex, rhs);
  BOOST_INSTALL_PROPERTY(vertex, key);
  BOOST_INSTALL_PROPERTY(vertex, position);
  BOOST_INSTALL_PROPERTY(vertex, second_bundle);

};

#endif

