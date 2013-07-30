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

#include <boost/container/vector.hpp>
#include <boost/container/list.hpp>

#include <boost/container/set.hpp>
#include <boost/unordered_set.hpp>


namespace boost {
  
namespace graph {

namespace detail {
  
  
/*************************************************************************
 *        pool-container values: validity and retrieval
 * **********************************************************************/
  
  struct BC_hole_desc {
    std::size_t value;
    explicit BC_hole_desc(std::size_t aValue = (std::numeric_limits<std::size_t>::max)()) : value(aValue) { };
    
    bool operator==(BC_hole_desc rhs) const { return (rhs.value == this->value); };
    bool operator!=(BC_hole_desc rhs) const { return (rhs.value != this->value); };
  };
  
  template <typename ValueType>
  bool BC_is_elem_valid(const ValueType&) { return true; };
  
  template <typename ValueType>
  bool BC_is_elem_valid(const variant< ValueType, BC_hole_desc >& vp) { 
    return (vp.which() == 0);
  };
  
  template <typename Container>
  struct BC_is_not_hole {
    const Container* p_cont;
    explicit BC_is_not_hole(const Container* aPCont = NULL) : p_cont(aPCont) { };
    template <typename Desc>
    bool operator()(Desc d) {
      return ((*p_cont)[d].which() == 0);
    };
  };
  
  template <typename Container>
  struct BC_edge_is_not_hole {
    const Container* p_cont;
    explicit BC_edge_is_not_hole(const Container* aPCont = NULL) : p_cont(aPCont) { };
    template <typename Desc>
    bool operator()(Desc e) {
      return ((*p_cont)[e.edge_id].which() == 0);
    };
  };
  
  template <typename ValueType>
  ValueType& BC_get_value(ValueType& vp) { return vp; };
  
  template <typename ValueType>
  ValueType& BC_get_value(variant< ValueType, BC_hole_desc >& vp) { 
    return get<ValueType>(vp);
  };
  
  template <typename ValueType>
  const ValueType& BC_get_value(const variant< ValueType, BC_hole_desc >& vp) { 
    return get<ValueType>(vp);
  };
  
  
  template <typename ValueType>
  struct BC_pooled_vector {
    typedef ValueType value_type;
    typedef variant< ValueType, BC_hole_desc > stored_type;
    typedef ::boost::container::vector<stored_type> container_type;
    typedef typename container_type::iterator iterator;
    
    ::boost::container::vector<stored_type> m_data;
    BC_hole_desc m_first_hole;
    std::size_t m_num_elements;
    
    BC_pooled_vector() : m_data(), m_first_hole(), m_num_elements() { };
    
    std::size_t size() const { return m_num_elements; };
    std::size_t capacity() const { return m_data.capacity(); };
    void clear() { 
      m_data.clear();
      m_first_hole = BC_hole_desc();
      m_num_elements = 0;
    };
  };
  
}; // detail

}; // graph
  
  
  struct vecBC { };
  struct poolBC { };
  struct listBC { };
  struct setBC { };
  struct multisetBC { };
  struct mapBC { };
  struct multimapBC { };
  struct unordered_setBC { };
  struct unordered_multisetBC { };
  struct unordered_mapBC { };
  struct unordered_multimapBC { };
  
  
  template <class Selector, class ValueType>
  struct BC_container_gen { };
  
  template <class ValueType>
  struct BC_container_gen<vecBC, ValueType> {
    typedef ::boost::container::vector<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<poolBC, ValueType> {
    typedef BC_pooled_vector<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<listBC, ValueType> {
    typedef ::boost::container::list<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<setBC, ValueType> {
    typedef ::boost::container::set<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<multisetBC, ValueType> {
    typedef ::boost::container::multiset<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<mapBC, ValueType> {
    typedef ::boost::container::set<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<multimapBC, ValueType> {
    typedef ::boost::container::multiset<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<unordered_setBC, ValueType> {
    typedef ::boost::unordered_set<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<unordered_multisetBC, ValueType> {
    typedef ::boost::unordered_multiset<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<unordered_mapBC, ValueType> {
    typedef ::boost::unordered_set<ValueType> type;
  };
  
  template <class ValueType>
  struct BC_container_gen<unordered_multimapBC, ValueType> {
    typedef ::boost::unordered_multiset<ValueType> type;
  };
  
  
  
namespace graph {

namespace detail {
  
  
  template <typename DescType>
  struct BC_null_desc {
    static DescType value() { return DescType(); };
  };
  
  template <typename DescType>
  struct BC_null_desc<DescType*> {
    static DescType* value() { return NULL; };
  };
  
  template <>
  struct BC_null_desc<std::size_t> {
    static std::size_t value() { return (std::numeric_limits<std::size_t>::max)(); };
  };
  
  
  bool BC_desc_less_than(std::size_t lhs, std::size_t rhs) {
    return (lhs < rhs);
  };
  
  template <typename Iter>
  bool BC_desc_less_than(Iter lhs, Iter rhs) {
    return (&(*lhs) < &(*rhs));
  };
  
  
/*************************************************************************
 *                      edge descriptors
 * **********************************************************************/
  
  template <typename Vertex, typename EdgeRawDesc>
  struct BC_edge_desc {
    typedef EdgeRawDesc edge_id_type;
    typedef Vertex source_descriptor;
    
    Vertex source;
    edge_id_type edge_id;
    
    BC_edge_desc(Vertex aSrc = BC_null_desc<Vertex>::value(), 
                 edge_id_type aEdgeId = BC_null_desc<edge_id_type>::value()) : 
                 source(aSrc), edge_id(aEdgeId) { };
    
    static BC_edge_desc null_value() { return BC_edge_desc(); };
  };
  
  template <typename Vertex, typename EdgeRawDesc>
  bool operator==(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                  const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (lhs.source == rhs.source) && (lhs.edge_id == rhs.edge_id);
  };
  template <typename Vertex, typename EdgeRawDesc>
  bool operator!=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                  const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (lhs.source != rhs.source) || (lhs.edge_id != rhs.edge_id);
  };
  template <typename Vertex, typename EdgeRawDesc>
  bool operator >(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                  const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (BC_desc_less_than(rhs.source, lhs.source) ||
           ((lhs.source == rhs.source) && ( BC_desc_less_than(rhs.edge_id, lhs.edge_id) )));
  };
  template <typename Vertex, typename EdgeRawDesc>
  bool operator >=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                   const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (BC_desc_less_than(rhs.source, lhs.source) ||
           ((lhs.source == rhs.source) && ( ! BC_desc_less_than(lhs.edge_id, rhs.edge_id) )));
  };
  
  template <typename Vertex, typename EdgeRawDesc>
  bool operator <(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                  const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (BC_desc_less_than(lhs.source, rhs.source) ||
           ((lhs.source == rhs.source) && ( BC_desc_less_than(lhs.edge_id, rhs.edge_id) )));
  };
  template <typename Vertex, typename EdgeRawDesc>
  bool operator <=(const BC_edge_desc<Vertex, EdgeRawDesc>& lhs, 
                   const BC_edge_desc<Vertex, EdgeRawDesc>& rhs) { 
    return (BC_desc_less_than(lhs.source, rhs.source) ||
           ((lhs.source == rhs.source) && ( ! BC_desc_less_than(rhs.edge_id, lhs.edge_id) )));
  };
  
  
  
/*************************************************************************
 *        iterator / descriptor translation functions for all container types
 * **********************************************************************/
  
  template <typename Container>
  typename Container::iterator BC_desc_to_iterator(Container& c, std::size_t d) {
    return c.begin() + d;
  };
  
  template <typename Container>
  typename Container::const_iterator BC_desc_to_iterator(const Container& c, std::size_t d) {
    return c.begin() + d;
  };
  
  template <typename Container, typename Iter>
  Iter BC_desc_to_iterator(const Container&, Iter it) { return it; };
  
  template <typename ValueType>
  std::size_t BC_iterator_to_desc( ::boost::container::vector<ValueType>& c, typename ::boost::container::vector<ValueType>::iterator it) {
    return it - c.begin();
  };
  
  template <typename ValueType>
  std::size_t BC_iterator_to_desc(const ::boost::container::vector<ValueType>& c, typename ::boost::container::vector<ValueType>::const_iterator it) {
    return it - c.begin();
  };
  
  template <typename Container, typename Iter>
  Iter BC_iterator_to_desc(const Container&, Iter it) { return it; };
  
  
  
  template <typename ValueType>
  std::size_t BC_get_begin_desc(const ::boost::container::vector<ValueType>& c) {
    return 0;
  };
  template <typename ValueType>
  std::size_t BC_get_begin_desc(const BC_pooled_vector<ValueType>& c) {
    return 0;
  };
  template <typename Container>
  typename Container::iterator BC_get_begin_desc(Container& c) {
    return c.begin();
  };
  template <typename Container>
  typename Container::iterator BC_get_begin_desc(Container* c) {
    return BC_get_begin_desc(*c);
  };
  
  template <typename ValueType>
  std::size_t BC_get_end_desc(const ::boost::container::vector<ValueType>& c) {
    return c.size();
  };
  template <typename ValueType>
  std::size_t BC_get_end_desc(const BC_pooled_vector<ValueType>& c) {
    return c.m_data.size();
  };
  template <typename Container>
  typename Container::iterator BC_get_end_desc(Container& c) {
    return c.end();
  };
  template <typename Container>
  typename Container::iterator BC_get_end_desc(Container* c) {
    return BC_get_end_desc(*c);
  };
  
  
/*************************************************************************
 *        descriptor / value translation functions for all container types
 * **********************************************************************/
  
  template <typename ValueType>
  ValueType& BC_desc_to_value(::boost::container::vector<ValueType>& c, std::size_t d) {
    return c[d];
  };
  template <typename ValueType>
  const ValueType& BC_desc_to_value(const ::boost::container::vector<ValueType>& c, std::size_t d) {
    return c[d];
  };
  
  template <typename ValueType>
  ValueType& BC_desc_to_value(BC_pooled_vector<ValueType>& c, std::size_t d) {
    return get<ValueType>(c.m_data[d]);
  };
  template <typename ValueType>
  const ValueType& BC_desc_to_value(const BC_pooled_vector<ValueType>& c, std::size_t d) {
    return get<ValueType>(c.m_data[d]);
  };
  
  template <typename Container>
  typename Container::value_type& BC_desc_to_value(const Container&, typename Container::iterator it) { 
    return *it;
  };
  
  template <typename Container, typename Desc>
  typename Container::value_type& BC_desc_to_value(Container* p_c, Desc d) { 
    return BC_desc_to_value(*p_c, d);
  };
  
  
  
  template <typename ListS, typename Iter>
  struct BC_select_descriptor {
    typedef Iter type;
  };
  
  template <typename Iter>
  struct BC_select_descriptor<vecBC, Iter> {
    typedef std::size_t type;
  };
  
  template <typename Iter>
  struct BC_select_descriptor<poolBC, Iter> {
    typedef std::size_t type;
  };
  
  
  template <typename ListS, typename ValueType>
  struct BC_select_value_type {
    typedef ValueType type;
  };
  
  template <typename ValueType>
  struct BC_select_value_type<poolBC, ValueType> {
    typedef variant< ValueType, BC_hole_desc > type;
  };
  
  
}; // namespace detail

}; // namespace graph
  
};


#endif


















