//=======================================================================
// Copyright 2014 Sven Mikael Persson
// Authors: Sven Mikael Persson 
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
#ifndef BFL_D_ARY_TREE_SERIALIZE_HPP
#define BFL_D_ARY_TREE_SERIALIZE_HPP

#include <boost/graph/bfl_d_ary_tree.hpp>

#include <boost/graph/detail/any_tree_serialize.hpp>

namespace boost { 

namespace serialization {

// Turn off tracking for bfl_d_ary_tree. It's not polymorphic, and we
// need to do this to enable saving of non-const adjacency lists.
template <std::size_t Arity, typename VP, typename EP>
struct tracking_level<boost::bfl_d_ary_tree<Arity,VP,EP> > {
  typedef mpl::integral_c_tag tag;
  typedef mpl::int_<track_never> type;
  BOOST_STATIC_CONSTANT(int, value = tracking_level::type::value);
};

template <typename Archive, std::size_t Arity, typename VP, typename EP>
inline void save(
    Archive & ar,
    const boost::bfl_d_ary_tree<Arity,VP,EP> &graph,
    const unsigned int file_version
){
  boost::graph::detail::save_tree(ar,graph,file_version);
}

template <typename Archive, std::size_t Arity, typename VP, typename EP>
inline void load(
    Archive & ar,
    boost::bfl_d_ary_tree<Arity,VP,EP> &graph,
    const unsigned int file_version
){
  boost::graph::detail::load_tree(ar,graph,file_version);
}

template <typename Archive, std::size_t Arity, typename VP, typename EP>
inline void serialize(
    Archive & ar,
    boost::bfl_d_ary_tree<Arity,VP,EP> &graph,
    const unsigned int file_version
){
  boost::serialization::split_free(ar, graph, file_version);
}

}//serialization
}//boost


#endif // BFL_D_ARY_TREE_SERIALIZE_HPP
