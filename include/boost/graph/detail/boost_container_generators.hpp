// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file boost_container_generators.hpp
 * 
 * This library declares a few useful boost-container generators for graph classes like adjacency_list or linked_tree. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_BOOST_CONTAINER_GENERATORS_HPP
#define BOOST_BOOST_CONTAINER_GENERATORS_HPP

#include <boost/config.hpp>

#include <boost/variant.hpp>

#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/container/list.hpp>
#include <boost/container/vector.hpp>

#include <boost/container/set.hpp>
#include <boost/unordered_set.hpp>

namespace boost {

namespace graph::detail {

/*************************************************************************
 *        pool-container values: validity and retrieval
 * **********************************************************************/

struct BC_hole_desc {
  std::size_t value;
  explicit BC_hole_desc(std::size_t aValue) : value(aValue) {}

  BC_hole_desc() : BC_hole_desc((std::numeric_limits<std::size_t>::max)()) {}

  bool operator==(BC_hole_desc rhs) const { return (rhs.value == this->value); }
  bool operator!=(BC_hole_desc rhs) const { return (rhs.value != this->value); }
};

template <typename ValueType>
bool BC_is_elem_valid(const ValueType& /*unused*/) {
  return true;
}

template <typename ValueType>
bool BC_is_elem_valid(const variant<ValueType, BC_hole_desc>& vp) {
  return (vp.which() == 0);
}

template <typename Container>
struct BC_is_not_hole {
  const Container* p_cont;
  explicit BC_is_not_hole(const Container* aPCont) : p_cont(aPCont) {}
  BC_is_not_hole() : BC_is_not_hole(nullptr) {}
  template <typename Desc>
  bool operator()(Desc d) {
    return ((*p_cont)[d].which() == 0);
  }
};

template <typename Container>
struct BC_edge_is_not_hole {
  const Container* p_cont;
  explicit BC_edge_is_not_hole(const Container* aPCont) : p_cont(aPCont) {}
  BC_edge_is_not_hole() : BC_edge_is_not_hole(nullptr) {}
  template <typename Desc>
  bool operator()(Desc e) {
    return ((*p_cont)[e.edge_id].which() == 0);
  }
};

template <typename ValueType>
ValueType& BC_get_value(ValueType& vp) {
  return vp;
}

template <typename ValueType>
ValueType& BC_get_value(variant<ValueType, BC_hole_desc>& vp) {
  return get<ValueType>(vp);
}

template <typename ValueType>
const ValueType& BC_get_value(const variant<ValueType, BC_hole_desc>& vp) {
  return get<ValueType>(vp);
}

template <typename ValueType>
struct BC_pooled_vector {
  using value_type = ValueType;
  using size_type = std::size_t;
  using stored_type = variant<ValueType, BC_hole_desc>;
  using container_type = ::boost::container::vector<stored_type>;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  ::boost::container::vector<stored_type> m_data;
  BC_hole_desc m_first_hole;
  std::size_t m_num_elements{0};

  BC_pooled_vector() : m_data() {}

  void swap(BC_pooled_vector& rhs) {
    using std::swap;
    m_data.swap(rhs.m_data);
    swap(m_first_hole, rhs.m_first_hole);
    swap(m_num_elements, rhs.m_num_elements);
  }

  std::size_t size() const { return m_num_elements; }
  std::size_t capacity() const { return m_data.capacity(); }
  void clear() {
    m_data.clear();
    m_first_hole = BC_hole_desc();
    m_num_elements = 0;
  }
  iterator begin() { return m_data.begin(); }
  const_iterator begin() const { return m_data.begin(); }
  iterator end() { return m_data.end(); }
  const_iterator end() const { return m_data.end(); }
};

template <typename Container>
std::size_t BC_get_size(const Container& cont) {
  return cont.size();
}

template <typename Container>
std::size_t BC_get_size(Container* cont) {
  return cont->size();
}

template <typename Container>
std::size_t BC_get_capacity(const Container& cont) {
  return cont.capacity();
}

template <typename Container>
std::size_t BC_get_capacity(Container* cont) {
  return cont->capacity();
}

template <typename Container>
void BC_clear_all(Container& cont) {
  cont.clear();
}

template <typename Container>
void BC_clear_all(Container* cont) {
  cont->clear();
}

}  // namespace graph::detail

struct vecBC {};
struct poolBC {};
struct listBC {};
struct setBC {};
struct multisetBC {};
struct mapBC {};
struct multimapBC {};
struct unordered_setBC {};
struct unordered_multisetBC {};
struct unordered_mapBC {};
struct unordered_multimapBC {};

template <class Selector, class ValueType>
struct BC_container_gen {};

template <class ValueType>
struct BC_container_gen<vecBC, ValueType> {
  using type = ::boost::container::vector<ValueType>;
};

template <class ValueType>
struct BC_container_gen<poolBC, ValueType> {
  using type = graph::detail::BC_pooled_vector<ValueType>;
};

template <class ValueType>
struct BC_container_gen<listBC, ValueType> {
  using type = ::boost::container::list<ValueType>;
};

template <class ValueType>
struct BC_container_gen<setBC, ValueType> {
  using type = ::boost::container::set<ValueType>;
};

template <class ValueType>
struct BC_container_gen<multisetBC, ValueType> {
  using type = ::boost::container::multiset<ValueType>;
};

template <class ValueType>
struct BC_container_gen<mapBC, ValueType> {
  using type = ::boost::container::set<ValueType>;
};

template <class ValueType>
struct BC_container_gen<multimapBC, ValueType> {
  using type = ::boost::container::multiset<ValueType>;
};

template <class ValueType>
struct BC_container_gen<unordered_setBC, ValueType> {
  using type = ::boost::unordered_set<ValueType>;
};

template <class ValueType>
struct BC_container_gen<unordered_multisetBC, ValueType> {
  using type = ::boost::unordered_multiset<ValueType>;
};

template <class ValueType>
struct BC_container_gen<unordered_mapBC, ValueType> {
  using type = ::boost::unordered_set<ValueType>;
};

template <class ValueType>
struct BC_container_gen<unordered_multimapBC, ValueType> {
  using type = ::boost::unordered_multiset<ValueType>;
};

namespace detail {

template <class Selector>
struct parallel_edge_BC_traits {
  using type = allow_parallel_edge_tag;
};

template <>
struct parallel_edge_BC_traits<setBC> {
  using type = disallow_parallel_edge_tag;
};

template <>
struct parallel_edge_BC_traits<unordered_setBC> {
  using type = disallow_parallel_edge_tag;
};

// mapBC is obsolete, replaced with setBC
template <>
struct parallel_edge_BC_traits<mapBC> {
  using type = disallow_parallel_edge_tag;
};

template <>
struct parallel_edge_BC_traits<unordered_mapBC> {
  using type = disallow_parallel_edge_tag;
};

template <typename Selector>
struct is_random_access_BC : mpl::false_ {};

template <>
struct is_random_access_BC<vecBC> : mpl::true_ {};

}  // namespace detail

namespace graph::detail {

template <typename DescType>
struct BC_null_desc {
  static DescType value() { return DescType(); }
};

template <typename DescType>
struct BC_null_desc<DescType*> {
  static DescType* value() { return NULL; }
};

template <>
struct BC_null_desc<std::size_t> {
  static std::size_t value() {
    return (std::numeric_limits<std::size_t>::max)();
  }
};

inline bool BC_desc_less_than(std::size_t lhs, std::size_t rhs) {
  return (lhs < rhs);
}

template <typename Iter>
bool BC_desc_less_than(Iter lhs, Iter rhs) {
  return (&(*lhs) < &(*rhs));
}

inline std::size_t BC_desc_get_hash(std::size_t d) {
  ::boost::hash<std::size_t> hasher;
  return hasher(d);
}

template <typename Iter>
std::size_t BC_desc_get_hash(Iter it) {
  using ValueType = typename Iter::value_type;
  ::boost::hash<ValueType*> hasher;
  return hasher(&(*it));
}

struct BC_desc_hasher {
  std::size_t operator()(std::size_t d) const {
    ::boost::hash<std::size_t> hasher;
    return hasher(d);
  }
  template <typename Iter>
  std::size_t operator()(Iter it) const {
    using ValueType = typename Iter::value_type;
    ::boost::hash<ValueType*> hasher;
    return hasher(&(*it));
  }
};

/*************************************************************************
 *                      edge descriptors
 * **********************************************************************/

template <typename Vertex, typename EdgeRawDesc>
struct BC_edge_desc {
  using edge_id_type = EdgeRawDesc;
  using source_descriptor = Vertex;

  Vertex source;
  edge_id_type edge_id;

  explicit BC_edge_desc(
      Vertex aSrc, edge_id_type aEdgeId = BC_null_desc<edge_id_type>::value())
      : source(aSrc), edge_id(aEdgeId) {}

  BC_edge_desc() : BC_edge_desc(BC_null_desc<Vertex>::value()) {}

  static BC_edge_desc null_value() { return BC_edge_desc(); }
};

template <typename Vertex, typename EdgeRawDesc>
bool operator==(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
                const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  if ((lhs.source == BC_null_desc<Vertex>::value()) &&
      (rhs.source == BC_null_desc<Vertex>::value())) {
    return true;
  }
  return (lhs.source == rhs.source) && (lhs.edge_id == rhs.edge_id);
}
template <typename Vertex, typename EdgeRawDesc>
bool operator!=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
                const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  if ((lhs.source == BC_null_desc<Vertex>::value()) &&
      (rhs.source == BC_null_desc<Vertex>::value())) {
    return false;
  }
  return (lhs.source != rhs.source) || (lhs.edge_id != rhs.edge_id);
}
template <typename Vertex, typename EdgeRawDesc>
bool operator>(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
               const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (BC_desc_less_than(rhs.source, lhs.source) ||
          ((lhs.source == rhs.source) &&
           (BC_desc_less_than(rhs.edge_id, lhs.edge_id))));
}
template <typename Vertex, typename EdgeRawDesc>
bool operator>=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
                const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (BC_desc_less_than(rhs.source, lhs.source) ||
          ((lhs.source == rhs.source) &&
           (!BC_desc_less_than(lhs.edge_id, rhs.edge_id))));
}

template <typename Vertex, typename EdgeRawDesc>
bool operator<(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
               const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (BC_desc_less_than(lhs.source, rhs.source) ||
          ((lhs.source == rhs.source) &&
           (BC_desc_less_than(lhs.edge_id, rhs.edge_id))));
}
template <typename Vertex, typename EdgeRawDesc>
bool operator<=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs,
                const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (BC_desc_less_than(lhs.source, rhs.source) ||
          ((lhs.source == rhs.source) &&
           (!BC_desc_less_than(rhs.edge_id, lhs.edge_id))));
}

template <typename EdgeDesc>
struct BC_undir_edge_desc : EdgeDesc {

  bool is_reversed;

  explicit BC_undir_edge_desc(const EdgeDesc& aBase, bool aIsReversed = false)
      : EdgeDesc(aBase), is_reversed(aIsReversed) {}

  BC_undir_edge_desc() : BC_undir_edge_desc(EdgeDesc()) {}

  static BC_undir_edge_desc<EdgeDesc> null_value() {
    return BC_undir_edge_desc<EdgeDesc>(EdgeDesc::null_value());
  }
};

/*************************************************************************
 *        iterator / descriptor translation functions for all container types
 * **********************************************************************/

template <typename Container>
typename Container::iterator BC_desc_to_iterator(Container& c, std::size_t d) {
  return c.begin() + d;
}

template <typename Container>
typename Container::const_iterator BC_desc_to_iterator(const Container& c,
                                                       std::size_t d) {
  return c.begin() + d;
}

template <typename Container, typename Iter>
Iter BC_desc_to_iterator(const Container& /*unused*/, Iter it) {
  return it;
}

template <typename ValueType>
std::size_t BC_iterator_to_desc(
    const ::boost::container::vector<ValueType>& c,
    typename ::boost::container::vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t BC_iterator_to_desc(
    const ::boost::container::vector<ValueType>& c,
    typename ::boost::container::vector<ValueType>::const_iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t BC_iterator_to_desc(
    const BC_pooled_vector<ValueType>& c,
    typename BC_pooled_vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t BC_iterator_to_desc(
    const BC_pooled_vector<ValueType>& c,
    typename BC_pooled_vector<ValueType>::const_iterator it) {
  return it - c.begin();
}

template <typename Container, typename Iter>
Iter BC_iterator_to_desc(const Container& /*unused*/, Iter it) {
  return it;
}

template <typename Container>
typename Container::iterator BC_get_begin_iter(Container& c) {
  return c.begin();
}
template <typename Container>
typename Container::iterator BC_get_begin_iter(Container* c) {
  return BC_get_begin_iter(*c);
}

template <typename ValueType>
std::size_t BC_get_begin_desc(const ::boost::container::vector<ValueType>& c) {
  return 0;
}
template <typename ValueType>
std::size_t BC_get_begin_desc(const BC_pooled_vector<ValueType>& c) {
  return 0;
}
template <typename Container>
typename Container::iterator BC_get_begin_desc(Container& c) {
  return c.begin();
}
template <typename Container>
typename Container::iterator BC_get_begin_desc(Container* c) {
  return BC_get_begin_desc(*c);
}

template <typename ValueType>
std::size_t BC_get_end_desc(const ::boost::container::vector<ValueType>& c) {
  return c.size();
}
template <typename ValueType>
std::size_t BC_get_end_desc(const BC_pooled_vector<ValueType>& c) {
  return c.m_data.size();
}
template <typename Container>
typename Container::iterator BC_get_end_desc(Container& c) {
  return c.end();
}
template <typename Container>
typename Container::iterator BC_get_end_desc(Container* c) {
  return BC_get_end_desc(*c);
}
template <typename Container>
typename Container::iterator BC_get_end_iter(Container& c) {
  return c.end();
}
template <typename Container>
typename Container::iterator BC_get_end_iter(Container* c) {
  return BC_get_end_iter(*c);
}

/*************************************************************************
 *        descriptor / value translation functions for all container types
 * **********************************************************************/

template <typename ValueType>
ValueType& BC_desc_to_value(::boost::container::vector<ValueType>& c,
                            std::size_t d) {
  return c[d];
}
template <typename ValueType>
const ValueType& BC_desc_to_value(
    const ::boost::container::vector<ValueType>& c, std::size_t d) {
  return c[d];
}

template <typename ValueType>
ValueType& BC_desc_to_value(BC_pooled_vector<ValueType>& c, std::size_t d) {
  return get<ValueType>(c.m_data[d]);
}
template <typename ValueType>
const ValueType& BC_desc_to_value(const BC_pooled_vector<ValueType>& c,
                                  std::size_t d) {
  return get<ValueType>(c.m_data[d]);
}

template <typename Container>
typename Container::value_type& BC_desc_to_value(
    const Container& /*unused*/, typename Container::iterator it) {
  return *it;
}

template <typename Container, typename Desc>
typename Container::value_type& BC_desc_to_value(Container* p_c, Desc d) {
  return BC_desc_to_value(*p_c, d);
}

template <typename Container>
struct BC_select_descriptor {
  using type = typename Container::iterator;
};

template <typename Container>
struct BC_select_descriptor<Container*> {
  using type = typename Container::iterator;
};

template <typename ValueType>
struct BC_select_descriptor<::boost::container::vector<ValueType>> {
  using type = std::size_t;
};

template <typename ValueType>
struct BC_select_descriptor<BC_pooled_vector<ValueType>> {
  using type = std::size_t;
};

/* Dummy "ignore" output iterator that is used in remove-branch functions. */
struct ignore_output_iter {
  struct value_type {
    template <typename T>
    value_type& operator=(const T& /*unused*/) {
      return *this;
    }
    template <typename T>
    value_type& operator=(T&& /*unused*/) {
      return *this;
    }
  };
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::bidirectional_iterator_tag;

  ignore_output_iter& operator++() { return *this; }
  ignore_output_iter& operator++(int) { return *this; }
  ignore_output_iter& operator--() { return *this; }
  ignore_output_iter& operator--(int) { return *this; }

  value_type operator*() const { return {}; }
};

}  // namespace graph::detail

}  // namespace boost

#endif
