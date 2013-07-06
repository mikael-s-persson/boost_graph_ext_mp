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

  BOOST_DEF_PROPERTY(vertex, heuristic);
  BOOST_DEF_PROPERTY(vertex, rhs);
  BOOST_DEF_PROPERTY(vertex, key);
  BOOST_DEF_PROPERTY(vertex, position);
  BOOST_DEF_PROPERTY(vertex, second_bundle);

};

#endif

