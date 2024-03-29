// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file adjlistBC_containers.hpp
 * 
 * This library implements adjacency-list data structures that underpin the adjacency_list_v2 template. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_ADJLISTBC_CONTAINERS_HPP
#define BOOST_ADJLISTBC_CONTAINERS_HPP

#include <boost/config.hpp>

#include <boost/graph/detail/adjlistBC_iterators.hpp>
#include <boost/graph/detail/boost_container_generators.hpp>

#include <boost/graph/graph_selectors.hpp>  // for directedS, undirectedS, bidirectionalS.

#include <iterator>
#include <utility>

#include <boost/mpl/and.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/not.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace boost::graph::detail {

/*************************************************************************
 *        value-types for vertices and edges
 * **********************************************************************/

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_vertex_stored_type;  // forward declaration.

// this is fine because boost-containers allow incomplete types:
template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_vertex_config {
  using stored_type =
      adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>;
  using container = typename BC_container_gen<VertexListS, stored_type>::type;
  using value_type = typename container::value_type;
  using descriptor = typename BC_select_descriptor<container>::type;

  static descriptor null_vertex() { return BC_null_desc<descriptor>::value(); }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_edge_stored_type {
  using self = adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                          VertexProperties, EdgeProperties>;
  using VConfig = adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                          VertexProperties, EdgeProperties>;
  using vertex_descriptor = typename VConfig::descriptor;

  vertex_descriptor target;
  mutable EdgeProperties data;

  explicit adjlistBC_edge_stored_type(vertex_descriptor aTarget)
      : target(aTarget), data() {}
  adjlistBC_edge_stored_type(vertex_descriptor aTarget,
                             const EdgeProperties& aData)
      : target(aTarget), data(aData) {}
  adjlistBC_edge_stored_type()
      : adjlistBC_edge_stored_type(VConfig::null_vertex()) {}
  adjlistBC_edge_stored_type(vertex_descriptor aTarget, EdgeProperties&& aData)
      : target(aTarget), data(std::move(aData)) {}

  bool operator<(const self& rhs) const {
    return BC_desc_less_than(this->target, rhs.target);
  }
  bool operator<=(const self& rhs) const {
    return !BC_desc_less_than(rhs.target, this->target);
  }
  bool operator>(const self& rhs) const {
    return BC_desc_less_than(rhs.target, this->target);
  }
  bool operator>=(const self& rhs) const {
    return !BC_desc_less_than(this->target, rhs.target);
  }
  bool operator==(const self& rhs) const {
    return (this->target == rhs.target);
  }
  bool operator!=(const self& rhs) const {
    return (this->target != rhs.target);
  }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
void swap(adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& lhs,
          adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& rhs) {
  using std::swap;
  swap(lhs.target, rhs.target);
  swap(lhs.data, rhs.data);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
std::size_t hash_value(
    const adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& ep) {
  return BC_desc_get_hash(ep.target);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_edge_config {
  using vertex_descriptor =
      typename adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                       VertexProperties,
                                       EdgeProperties>::descriptor;

  using stored_type =
      adjlistBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>;
  using container = typename BC_container_gen<OutEdgeListS, stored_type>::type;
  using value_type = typename container::value_type;
  using raw_descriptor = typename BC_select_descriptor<container>::type;
  using container_ptr = typename mpl::if_<
      mpl::and_<is_same<vertex_descriptor, std::size_t>,
                mpl::not_<is_same<raw_descriptor, std::size_t>>>,
      container*, container>::type;
  using descriptor = BC_edge_desc<vertex_descriptor, raw_descriptor>;

  static descriptor null_edge() { return BC_null_desc<descriptor>::value(); }
};

/*************************************************************************
 *        vertex values, as stored in the vertex containers.
 * **********************************************************************/

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_vertex_stored_type {
  using self =
      adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>;
  using directed_tag = DirectedS;
  using Config = adjlistBC_edge_config<VertexListS, OutEdgeListS, DirectedS,
                                       VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_container_ptr = typename Config::container_ptr;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_container = ::boost::container::vector<edge_descriptor>;
  using in_edge_iterator = typename in_edge_container::iterator;

  VertexProperties data;
  edge_container_ptr out_edges;
  in_edge_container in_edges;

  adjlistBC_vertex_stored_type() : data(), out_edges(), in_edges() {}
  explicit adjlistBC_vertex_stored_type(const VertexProperties& aData)
      : data(aData), out_edges(), in_edges() {}
  explicit adjlistBC_vertex_stored_type(VertexProperties&& aData)
      : data(std::move(aData)), out_edges(), in_edges() {}
  void swap(self& rhs) {
    using std::swap;
    swap(data, rhs.data);
    swap(out_edges, rhs.out_edges);
    swap(in_edges, rhs.in_edges);
  }
  self& operator=(self rhs) {
    this->swap(rhs);
    return *this;
  }
};

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                    VertexProperties, EdgeProperties> {
  using self =
      adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                   VertexProperties, EdgeProperties>;
  using directed_tag = directedS;
  using Config = adjlistBC_edge_config<VertexListS, OutEdgeListS, directedS,
                                       VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_container_ptr = typename Config::container_ptr;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = int*;

  VertexProperties data;
  edge_container_ptr out_edges;

  adjlistBC_vertex_stored_type() : data(), out_edges() {}
  explicit adjlistBC_vertex_stored_type(const VertexProperties& aData)
      : data(aData), out_edges() {}
  explicit adjlistBC_vertex_stored_type(VertexProperties&& aData)
      : data(std::move(aData)), out_edges() {}
  void swap(self& rhs) {
    using std::swap;
    swap(data, rhs.data);
    swap(out_edges, rhs.out_edges);
  }
  self& operator=(self rhs) {
    this->swap(rhs);
    return *this;
  }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
void swap(adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                       VertexProperties, EdgeProperties>& lhs,
          adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                       VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

/*************************************************************************
 *        out-edges container pointers (or not): create and destroy
 * **********************************************************************/

template <typename OutEdgeCont>
struct adjlistBC_out_edges_factory {
  template <typename VertexValue>
  static void create_out_edges(VertexValue& /*unused*/) {}
  template <typename VertexValue>
  static void destroy_out_edges(VertexValue& /*unused*/) {}
  template <typename VertexContainer>
  static void destroy_all_out_edges(VertexContainer& /*unused*/) {}
};

template <typename OutEdgeCont>
struct adjlistBC_out_edges_factory<OutEdgeCont*> {
  template <typename VertexValue>
  static void create_out_edges(VertexValue& vp) {
    vp.out_edges = new OutEdgeCont();
  }
  template <typename VertexValue>
  static void destroy_out_edges(VertexValue& vp) {
    delete vp.out_edges;
  }
  template <typename VertexValue>
  static void destroy_all_out_edges(
      ::boost::container::vector<VertexValue>& vcont) {
    using VertexContainer = ::boost::container::vector<VertexValue>;
    using Iter = typename VertexContainer::iterator;
    for (Iter vi = vcont.begin(); vi != vcont.end(); ++vi) {
      delete vi->out_edges;
    }
  }
  template <typename VertexValue>
  static void destroy_all_out_edges(BC_pooled_vector<VertexValue>& vcont) {
    using VertexContainer =
        typename BC_pooled_vector<VertexValue>::container_type;
    using Iter = typename VertexContainer::iterator;
    for (Iter vi = vcont.m_data.begin(); vi != vcont.m_data.end(); ++vi) {
      if (vi->which() == 0) {
        delete get<VertexValue>(*vi).out_edges;
      }
    }
  }
};

template <typename OutEdgeCont>
struct adjlistBC_out_edges_range {
  using Iter = typename OutEdgeCont::iterator;
  static Iter begin(OutEdgeCont& econt) { return econt.begin(); }
  static Iter end(OutEdgeCont& econt) { return econt.end(); }
  template <typename Desc>
  static Iter from_desc(OutEdgeCont& econt, Desc d) {
    return BC_desc_to_iterator(econt, d);
  }
};

template <typename OutEdgeCont>
struct adjlistBC_out_edges_range<OutEdgeCont*> {
  using Iter = typename OutEdgeCont::iterator;
  static Iter begin(OutEdgeCont* econt) { return econt->begin(); }
  static Iter end(OutEdgeCont* econt) { return econt->end(); }
  template <typename Desc>
  static Iter from_desc(OutEdgeCont* econt, Desc d) {
    return BC_desc_to_iterator(*econt, d);
  }
};

/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void adjlistBC_erase_in_edge(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                 VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {}

// for vector of edge-desc (in-edges)
template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void adjlistBC_erase_in_edge(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {
  vp.in_edges.erase(std::find(vp.in_edges.begin(), vp.in_edges.end(), e));
}

// for OutEdgeListS = listS, setS, ...
template <typename Container, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void adjlistBC_erase_edge(Container& cont, EdgeDesc e, VertexCont& /*unused*/,
                          VertexDesc /*unused*/) {
  cont.erase(e.edge_id);
}
template <typename Container, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void adjlistBC_erase_edge(Container* cont, EdgeDesc e, VertexCont& vcont,
                          VertexDesc v) {
  adjlistBC_erase_edge(*cont, e, vcont, v);
}

// for OutEdgeListS = vecS
template <typename ValueType, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void adjlistBC_erase_edge(::boost::container::vector<ValueType>& cont,
                          EdgeDesc e, VertexCont& vertex_cont, VertexDesc v) {
  using Container = ::boost::container::vector<ValueType>;
  using Iter = typename Container::iterator;
  using std::swap;

  Iter it = BC_desc_to_iterator(cont, e.edge_id);
  Iter it_last = cont.end();
  --it_last;
  if (it != it_last) {
    swap(*it, *it_last);
    // If this graph has in-edge references, then they must be updated.
    adjlistBC_update_in_edge_id(
        BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), v,
        it_last - cont.begin(), it - cont.begin());
  }
  cont.erase(it_last, cont.end());
}

// for OutEdgeListS = poolS
template <typename ValueType, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void adjlistBC_erase_edge(BC_pooled_vector<ValueType>& cont, EdgeDesc e,
                          VertexCont& /*unused*/, VertexDesc /*unused*/) {
  cont.m_data[e.edge_id] = cont.m_first_hole;
  cont.m_first_hole = BC_hole_desc(e.edge_id);
  --(cont.m_num_elements);
}

/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void adjlistBC_add_in_edge(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                 VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void adjlistBC_add_in_edge(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {
  vp.in_edges.push_back(e);
}

// for OutEdgeListS = listS, multisetS, ...
template <typename Container, typename EdgeProperties, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_add_edge(Container& cont, const EdgeProperties& ep, VertexDesc v) {
  using ValueType = typename Container::value_type;
  return std::pair<typename Container::iterator, bool>(
      cont.insert(cont.end(), ValueType(v, ep)), true);
}

// for OutEdgeListS = setS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<typename ::boost::container::set<ValueType>::iterator, bool>
adjlistBC_add_edge(::boost::container::set<ValueType>& cont,
                   const EdgeProperties& ep, VertexDesc v) {
  return cont.insert(ValueType(v, ep));
}

// for OutEdgeListS = unordered_setS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<typename ::boost::unordered_set<ValueType>::iterator, bool>
adjlistBC_add_edge(::boost::unordered_set<ValueType>& cont,
                   const EdgeProperties& ep, VertexDesc v) {
  return cont.insert(ValueType(v, ep));
}

// for OutEdgeListS = vecS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_add_edge(
    ::boost::container::vector<ValueType>& cont, const EdgeProperties& ep,
    VertexDesc v) {
  cont.push_back(ValueType(v, ep));
  return std::pair<std::size_t, bool>(cont.size() - 1, true);
}

// for OutEdgeListS = poolS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_add_edge(
    BC_pooled_vector<ValueType>& cont, const EdgeProperties& ep, VertexDesc v) {
  using Iter = typename BC_pooled_vector<ValueType>::iterator;

  if (cont.m_first_hole == BC_hole_desc()) {
    Iter it = cont.m_data.insert(cont.m_data.end(), ValueType(v, ep));
    ++(cont.m_num_elements);
    return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
  }
  Iter it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(v, ep);
  ++(cont.m_num_elements);
  return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
}

template <typename Container, typename EdgeProperties, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_add_edge(Container* cont, const EdgeProperties& ep, VertexDesc v) {
  return adjlistBC_add_edge(*cont, ep, v);
}

// for OutEdgeListS = listS, multisetS, ...
template <typename Container, typename EdgeProperties, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_add_edge(Container& cont, EdgeProperties&& ep, VertexDesc v) {
  using ValueType = typename Container::value_type;
  return std::pair<typename BC_select_descriptor<Container>::type, bool>(
      cont.insert(cont.end(), ValueType(v, std::move(ep))), true);
}

// for OutEdgeListS = setS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<typename ::boost::container::set<ValueType>::iterator, bool>
adjlistBC_add_edge(::boost::container::set<ValueType>& cont,
                   EdgeProperties&& ep, VertexDesc v) {
  return cont.insert(ValueType(v, std::move(ep)));
}

// for OutEdgeListS = unordered_setS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<typename ::boost::unordered_set<ValueType>::iterator, bool>
adjlistBC_add_edge(::boost::unordered_set<ValueType>& cont, EdgeProperties&& ep,
                   VertexDesc v) {
  return cont.insert(ValueType(v, std::move(ep)));
}

// for OutEdgeListS = vecS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_add_edge(
    ::boost::container::vector<ValueType>& cont, EdgeProperties&& ep,
    VertexDesc v) {
  cont.push_back(ValueType(v, std::move(ep)));
  return std::pair<std::size_t, bool>(cont.size() - 1, true);
}

// for OutEdgeListS = poolS
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_add_edge(
    BC_pooled_vector<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
  if (cont.m_first_hole == BC_hole_desc()) {
    auto it =
        cont.m_data.insert(cont.m_data.end(), ValueType(v, std::move(ep)));
    ++(cont.m_num_elements);
    return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
  }
  auto it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(v, std::move(ep));
  ++(cont.m_num_elements);
  return std::pair<std::size_t, bool>(it - cont.m_data.begin(), true);
}

template <typename Container, typename EdgeProperties, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_add_edge(Container* cont, EdgeProperties&& ep, VertexDesc v) {
  return adjlistBC_add_edge(*cont, std::move(ep), v);
}

/*************************************************************************
 *        helper functions for finding an edgee between two vertices
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC    (multi)setBC      unordered_(multi)setBC
// any-dir.:    O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)

// for OutEdgeListS = listS
template <typename ValueType, typename VertexDesc>
std::pair<typename ::boost::container::list<ValueType>::iterator, bool>
adjlistBC_find_edge_to(::boost::container::list<ValueType>& cont,
                       VertexDesc v) {
  using Iter = typename ::boost::container::list<ValueType>::iterator;
  using ResultType = std::pair<Iter, bool>;
  Iter it = cont.begin();
  for (; it != cont.end(); ++it) {
    if (it->target == v) {
      return ResultType(it, true);
    }
  }
  return ResultType(it, false);
}

// for OutEdgeListS = setS, multisetS, ...
template <typename Container, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_find_edge_to(Container& cont, VertexDesc v) {
  using ValueType = typename Container::value_type;
  using Iter = typename Container::iterator;
  using ResultType = std::pair<Iter, bool>;
  auto ep = ValueType(v);
  Iter it = cont.find(ep);
  return ResultType(it, (it != cont.end()));
}

// for OutEdgeListS = vecS
template <typename ValueType, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_find_edge_to(
    ::boost::container::vector<ValueType>& cont, VertexDesc v) {
  using Container = ::boost::container::vector<ValueType>;
  using Iter = typename Container::iterator;
  using ResultType = std::pair<std::size_t, bool>;

  for (Iter it = cont.begin(); it != cont.end(); ++it) {
    if (it->target == v) {
      return ResultType(it - cont.begin(), true);
    }
  }
  return ResultType(cont.size(), false);
}

// for OutEdgeListS = poolS
template <typename ValueType, typename VertexDesc>
std::pair<std::size_t, bool> adjlistBC_find_edge_to(
    BC_pooled_vector<ValueType>& cont, VertexDesc v) {
  using Container = typename BC_pooled_vector<ValueType>::container_type;
  using Iter = typename Container::iterator;
  using ResultType = std::pair<std::size_t, bool>;

  for (Iter it = cont.m_data.begin(); it != cont.m_data.end(); ++it) {
    if ((it->which() == 0) && (get<ValueType>(*it).target == v)) {
      return ResultType(it - cont.m_data.begin(), true);
    }
  }
  return ResultType(cont.m_data.size(), false);
}

template <typename Container, typename VertexDesc>
std::pair<typename BC_select_descriptor<Container>::type, bool>
adjlistBC_find_edge_to(Container* cont, VertexDesc v) {
  return adjlistBC_find_edge_to(*cont, v);
}

/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC    (multi)setBC      unordered_(multi)setBC
// directedS:   O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)    O(log(E/V))       O(1)

// for OutEdgeListS = listS
template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlistBC_erase_edges_to(::boost::container::list<ValueType>& cont,
                              VertexCont& /*unused*/, VertexDesc /*unused*/,
                              VertexDesc v, std::size_t& e_count) {
  using Iter = typename ::boost::container::list<ValueType>::iterator;
  for (Iter it = cont.begin(); it != cont.end();) {
    if (it->target == v) {
      it = cont.erase(it);
      --e_count;
    } else {
      ++it;
    }
  }
}

// for OutEdgeListS = setS, multisetS, ...
template <typename Container, typename VertexCont, typename VertexDesc>
void adjlistBC_erase_edges_to(Container& cont, VertexCont& /*unused*/,
                              VertexDesc /*unused*/, VertexDesc v,
                              std::size_t& e_count) {
  using ValueType = typename Container::value_type;
  auto ep = ValueType(v);
  e_count -= cont.erase(ep);
}

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlistBC_update_in_edge_id(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                 VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlistBC_update_in_edge_id(
    adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {
  using VertexValue =
      adjlistBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>;
  using InEdgeCont = typename VertexValue::in_edge_container;
  using InEdgeIter = typename InEdgeCont::iterator;

  for (InEdgeIter ei = vp.in_edges.begin(); ei != vp.in_edges.end(); ++ei) {
    if ((ei->source == v) && (ei->edge_id == old_id)) {
      ei->edge_id = new_id;
      break;
    }
  }
}

// for OutEdgeListS = vecS
template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlistBC_erase_edges_to(::boost::container::vector<ValueType>& cont,
                              VertexCont& vertex_cont, VertexDesc u,
                              VertexDesc v, std::size_t& e_count) {
  using Container = ::boost::container::vector<ValueType>;
  using Iter = typename Container::iterator;
  using std::swap;

  Iter it_last = cont.end();
  for (Iter it = cont.begin(); it != it_last;) {
    if (BC_get_value(*it).target == v) {
      --it_last;
      if (it != it_last) {
        swap(*it, *it_last);
        // If this graph has in-edge references, then they must be updated.
        adjlistBC_update_in_edge_id(
            BC_get_value(*BC_desc_to_iterator(vertex_cont, it->target)), u,
            it_last - cont.begin(), it - cont.begin());
      }
      --e_count;
    } else {
      ++it;
    }
  }
  cont.erase(it_last, cont.end());
}

// for OutEdgeListS = poolS
template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlistBC_erase_edges_to(BC_pooled_vector<ValueType>& cont,
                              VertexCont& /*unused*/, VertexDesc /*unused*/,
                              VertexDesc v, std::size_t& e_count) {
  using Container = typename BC_pooled_vector<ValueType>::container_type;
  using Iter = typename Container::iterator;

  for (Iter it = cont.m_data.begin(); it != cont.m_data.end(); ++it) {
    if ((it->which() == 0) && (BC_get_value(*it).target == v)) {
      *it = cont.m_first_hole;
      cont.m_first_hole = BC_hole_desc(it - cont.m_data.begin());
      --(cont.m_num_elements);
      --e_count;
    }
  }
}

template <typename Container, typename VertexCont, typename VertexDesc>
void adjlistBC_erase_edges_to(Container* cont, VertexCont& vcont, VertexDesc u,
                              VertexDesc v, std::size_t& e_count) {
  adjlistBC_erase_edges_to(*cont, vcont, u, v, e_count);
}

/*************************************************************************
 *        functions for clearing the edges of a vertex
 * **********************************************************************/
// NOTE: The function 'adjlistBC_clear_vertex' works for all graph types.

// NOTE: Time complexities:
//                 vecBC         poolBC         listBC
// directedS:   O(V*(E/V))     O(V*(E/V))     O(V*(E/V))
// bidir:       O((E/V)^2)     O((E/V)^2)     O((E/V)^2)

template <typename DirectedS, typename VertexCont, typename VertexValue,
          typename VertexDesc>
typename enable_if<is_same<DirectedS, directedS>, void>::type
adjlistBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v,
                       std::size_t& e_count) {
  using VertexIter = typename VertexCont::iterator;

  // first, just clear the out-going edges. No need to synchronize in-edges (there are none).
  e_count -= BC_get_size(vp.out_edges);
  BC_clear_all(vp.out_edges);
  VertexIter vi = BC_desc_to_iterator(cont, v);

  // now, the stupid part... we have to traverse all other vertices (and their edges)
  // to look for in-edges that lead to "v", and erase them.
  for (VertexIter ui = cont.begin(); ui != cont.end(); ++ui) {
    if ((ui == vi) || !BC_is_elem_valid(*ui)) {
      continue;
    }
    VertexValue& up = BC_get_value(*ui);
    adjlistBC_erase_edges_to(up.out_edges, cont, BC_iterator_to_desc(cont, ui),
                             v, e_count);
  }
}

template <typename DirectedS, typename VertexCont, typename VertexValue,
          typename VertexDesc>
typename disable_if<is_same<DirectedS, directedS>, void>::type
adjlistBC_clear_vertex(VertexCont& cont, VertexValue& vp, VertexDesc v,
                       std::size_t& e_count) {
  using OutEdgeCont = typename VertexValue::edge_container;
  using OutEdgeIter = typename OutEdgeCont::iterator;
  using InEdgeCont = typename VertexValue::in_edge_container;
  using InEdgeIter = typename InEdgeCont::iterator;

  // first, remove the in-edge references from the adjacent vertices of v.
  for (OutEdgeIter ei = BC_get_begin_iter(vp.out_edges);
       ei != BC_get_end_iter(vp.out_edges); ++ei) {
    if (!BC_is_elem_valid(*ei)) {
      continue;
    }
    VertexValue& wp =
        BC_get_value(*BC_desc_to_iterator(cont, BC_get_value(*ei).target));
    for (InEdgeIter iei = wp.in_edges.begin(); iei != wp.in_edges.end();
         ++iei) {
      if ((iei->source == v) &&
          (BC_desc_to_iterator(vp.out_edges, iei->edge_id) == ei)) {
        wp.in_edges.erase(iei);
        break;
      }
    }
  }

  // then, clear the out-going edges.
  e_count -= BC_get_size(vp.out_edges);
  BC_clear_all(vp.out_edges);

  // finally, remove the required out-edges of the "parent" vertices of v.
  for (InEdgeIter iei = vp.in_edges.begin(); iei != vp.in_edges.end(); ++iei) {
    VertexDesc u = iei->source;
    VertexValue& up = BC_get_value(*BC_desc_to_iterator(cont, iei->source));
    adjlistBC_erase_edges_to(up.out_edges, cont, u, v, e_count);
  }
  vp.in_edges.clear();
}

/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
// NOTE: The function 'adjlistBC_erase_vertex' works for all graph types.
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

template <typename Container, typename VertexDesc>
void adjlistBC_erase_vertex(Container& cont, VertexDesc v) {
  using VertexValue = typename Container::value_type;
  using OutEdgeFactory =
      adjlistBC_out_edges_factory<typename VertexValue::edge_container_ptr>;
  OutEdgeFactory::destroy_out_edges(*v);
  cont.erase(v);
}

template <typename Container>
void adjlistBC_update_out_edges_impl(Container& cont, std::size_t old_v_id,
                                     std::size_t new_v_id) {
  using OEIter = typename Container::iterator;
  for (OEIter ei = cont.begin(); ei != cont.end(); ++ei) {
    if (BC_is_elem_valid(*ei) && (BC_get_value(*ei).target == old_v_id)) {
      BC_get_value(*ei).target = new_v_id;
    }
  }
}

template <typename Container>
void adjlistBC_update_assoc_out_edges_impl(Container& cont,
                                           std::size_t old_v_id,
                                           std::size_t new_v_id) {
  using OEIter = typename Container::iterator;
  using ValueType = typename Container::value_type;
  using std::inserter;
  auto v_test = ValueType(old_v_id);
  std::pair<OEIter, OEIter> eq_rg = cont.equal_range(v_test);
  std::vector<ValueType> v_temp;
  std::copy(eq_rg.first, eq_rg.second, std::back_inserter(v_temp));
  for (auto it = v_temp.begin(); it != v_temp.end(); ++it) {
    it->target = new_v_id;
  }
  cont.erase(eq_rg.first, eq_rg.second);
  std::move(v_temp.begin(), v_temp.end(), inserter(cont, cont.end()));
}

template <typename ValueType>
void adjlistBC_update_out_edges_impl(::boost::container::set<ValueType>& cont,
                                     std::size_t old_v_id,
                                     std::size_t new_v_id) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::container::multiset<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(::boost::unordered_set<ValueType>& cont,
                                     std::size_t old_v_id,
                                     std::size_t new_v_id) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::unordered_multiset<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}

template <typename Container>
void adjlistBC_update_out_edges_impl(Container* cont, std::size_t old_v_id,
                                     std::size_t new_v_id) {
  adjlistBC_update_out_edges_impl(*cont, old_v_id, new_v_id);
}

template <typename Container>
void adjlistBC_update_out_edges_impl(Container& cont, std::size_t old_v_id,
                                     std::size_t new_v_id,
                                     typename Container::iterator ei) {
  BC_get_value(*ei).target = new_v_id;
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::container::set<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id,
    typename ::boost::container::set<ValueType>::iterator /*unused*/) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::container::multiset<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id,
    typename ::boost::container::multiset<ValueType>::iterator /*unused*/) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::unordered_set<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id,
    typename ::boost::unordered_set<ValueType>::iterator /*unused*/) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlistBC_update_out_edges_impl(
    ::boost::unordered_multiset<ValueType>& cont, std::size_t old_v_id,
    std::size_t new_v_id,
    typename ::boost::unordered_multiset<ValueType>::iterator /*unused*/) {
  adjlistBC_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename Container>
void adjlistBC_update_out_edges_impl(Container* cont, std::size_t old_v_id,
                                     std::size_t new_v_id,
                                     typename Container::iterator ei) {
  adjlistBC_update_out_edges_impl(*cont, old_v_id, new_v_id, ei);
}

// O(E)
template <typename DirectedS, typename ValueType>
typename enable_if<is_same<DirectedS, directedS>, void>::type
adjlistBC_update_out_edges(::boost::container::vector<ValueType>& cont,
                           std::size_t old_v_id, std::size_t new_v_id) {
  using VContainer = ::boost::container::vector<ValueType>;
  using VIter = typename VContainer::iterator;

  for (VIter vi = cont.begin(); vi != cont.end(); ++vi) {
    adjlistBC_update_out_edges_impl(vi->out_edges, old_v_id, new_v_id);
  }
}

// O((E/V)^2)
template <typename DirectedS, typename ValueType>
typename disable_if<is_same<DirectedS, directedS>, void>::type
adjlistBC_update_out_edges(::boost::container::vector<ValueType>& cont,
                           std::size_t old_v_id, std::size_t new_v_id) {
  using OutEdgeCont = typename ValueType::edge_container;
  using OEIter = typename OutEdgeCont::iterator;
  using InEdgeCont = typename ValueType::in_edge_container;
  using InEdgeIter = typename InEdgeCont::iterator;

  // first, update in-edge vertices
  for (InEdgeIter iei = cont[new_v_id].in_edges.begin();
       iei != cont[new_v_id].in_edges.end(); ++iei) {
    ValueType& up = cont[iei->source];
    adjlistBC_update_out_edges_impl(
        up.out_edges, old_v_id, new_v_id,
        BC_desc_to_iterator(up.out_edges, iei->edge_id));
  }

  // second, update out-edge vertices
  for (OEIter ei = BC_get_begin_iter(cont[new_v_id].out_edges);
       ei != BC_get_end_iter(cont[new_v_id].out_edges); ++ei) {
    if (!BC_is_elem_valid(*ei)) {
      continue;
    }
    ValueType& wp = cont[BC_get_value(*ei).target];
    for (InEdgeIter iei = wp.in_edges.begin(); iei != wp.in_edges.end();
         ++iei) {
      if ((iei->source == old_v_id) &&
          (ei == BC_desc_to_iterator(cont[new_v_id].out_edges, iei->edge_id))) {
        iei->source = new_v_id;
        break;
      }
    }
  }
}

template <typename ValueType>
void adjlistBC_erase_vertex(::boost::container::vector<ValueType>& cont,
                            std::size_t v) {
  using Container = ::boost::container::vector<ValueType>;
  using Iter = typename Container::iterator;
  using OutEdgeFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  using DirectedS = typename ValueType::directed_tag;
  using std::swap;

  Iter it = BC_desc_to_iterator(cont, v);
  OutEdgeFactory::destroy_out_edges(*it);
  Iter it_last = cont.end();
  --it_last;
  if (it == it_last) {
    cont.erase(it_last);
    return;
  }
  swap(*it, *it_last);
  std::size_t old_id = it_last - cont.begin();
  std::size_t new_id = it - cont.begin();
  cont.erase(it_last);
  adjlistBC_update_out_edges<DirectedS>(cont, old_id, new_id);
}

template <typename ValueType>
void adjlistBC_erase_vertex(BC_pooled_vector<ValueType>& cont, std::size_t v) {
  using OutEdgeFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  // the first_hole will become v and v will be referring to first_hole (i.e., like a simple linked-list):
  OutEdgeFactory::destroy_out_edges(get<ValueType>(cont.m_data[v]));
  cont.m_data[v] = cont.m_first_hole;
  cont.m_first_hole = BC_hole_desc(v);
  --(cont.m_num_elements);
}

/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

template <typename Container, typename VertexProperties>
typename BC_select_descriptor<Container>::type adjlistBC_add_vertex(
    Container& cont, const VertexProperties& vp) {
  using ValueType = typename Container::value_type;
  using Iter = typename Container::iterator;
  Iter it = cont.insert(cont.end(), ValueType(vp));
  using OEFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  OEFactory::create_out_edges(*it);
  return it;
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlistBC_add_vertex(::boost::container::vector<ValueType>& cont,
                                 const VertexProperties& vp) {
  auto it = cont.insert(cont.end(), ValueType(vp));
  using OEFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  OEFactory::create_out_edges(*it);
  return it - cont.begin();
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlistBC_add_vertex(BC_pooled_vector<ValueType>& cont,
                                 const VertexProperties& vp) {
  using OEFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;

  if (cont.m_first_hole == BC_hole_desc()) {
    auto it = cont.m_data.insert(cont.m_data.end(), ValueType(vp));
    ++(cont.m_num_elements);
    OEFactory::create_out_edges(get<ValueType>(*it));
    return it - cont.m_data.begin();
  }
  auto it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(vp);
  ++(cont.m_num_elements);
  OEFactory::create_out_edges(get<ValueType>(*it));
  return it - cont.m_data.begin();
}

template <typename Container, typename VertexProperties>
typename BC_select_descriptor<Container>::type adjlistBC_add_vertex(
    Container* cont, const VertexProperties& vp) {
  return adjlistBC_add_vertex(*cont, vp);
}

template <typename Container, typename VertexProperties>
typename BC_select_descriptor<Container>::type adjlistBC_add_vertex(
    Container& cont, VertexProperties&& vp) {
  using ValueType = typename Container::value_type;
  return cont.insert(cont.end(), ValueType(std::move(vp)));
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlistBC_add_vertex(::boost::container::vector<ValueType>& cont,
                                 VertexProperties&& vp) {
  using Iter = typename ::boost::container::vector<ValueType>::iterator;
  Iter it = cont.insert(cont.end(), ValueType(std::move(vp)));
  using OEFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  OEFactory::create_out_edges(*it);
  return it - cont.begin();
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlistBC_add_vertex(BC_pooled_vector<ValueType>& cont,
                                 VertexProperties&& vp) {
  using OEFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;

  if (cont.m_first_hole == BC_hole_desc()) {
    auto it = cont.m_data.insert(cont.m_data.end(), ValueType(std::move(vp)));
    ++(cont.m_num_elements);
    OEFactory::create_out_edges(get<ValueType>(*it));
    return it - cont.m_data.begin();
  }
  auto it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = get<BC_hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(std::move(vp));
  ++(cont.m_num_elements);
  OEFactory::create_out_edges(get<ValueType>(*it));
  return it - cont.m_data.begin();
}

template <typename Container, typename VertexProperties>
typename BC_select_descriptor<Container>::type adjlistBC_add_vertex(
    Container* cont, VertexProperties&& vp) {
  return adjlistBC_add_vertex(*cont, std::move(vp));
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct adjlistBC_vertex_container {

  using VConfig = adjlistBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                          VertexProperties, EdgeProperties>;
  using vertex_container = typename VConfig::container;
  using vertices_size_type = typename vertex_container::size_type;
  using vertex_descriptor = typename VConfig::descriptor;
  using vertex_stored_type = typename VConfig::stored_type;
  using vertex_value_type = typename VConfig::value_type;

  using EConfig = adjlistBC_edge_config<VertexListS, OutEdgeListS, DirectedS,
                                        VertexProperties, EdgeProperties>;
  using edge_container = typename EConfig::container;
  using edges_size_type = typename edge_container::size_type;
  using edge_descriptor = typename EConfig::descriptor;
  using edge_stored_type = typename EConfig::stored_type;
  using edge_value_type = typename EConfig::value_type;

  mutable vertex_container m_vertices;
  std::size_t m_num_edges{0};

  adjlistBC_vertex_container() : m_vertices() {}

  ~adjlistBC_vertex_container() { clear(); }

  adjlistBC_vertex_container(const adjlistBC_vertex_container&) = delete;
  adjlistBC_vertex_container& operator=(const adjlistBC_vertex_container&) =
      delete;

  void swap(adjlistBC_vertex_container& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    // swap(m_vertices, rhs.m_vertices);
    swap(m_num_edges, rhs.m_num_edges);
  }

  adjlistBC_vertex_container(adjlistBC_vertex_container&& rhs) noexcept
      : m_vertices() {
    swap(rhs);
  }
  adjlistBC_vertex_container& operator=(
      adjlistBC_vertex_container&& rhs) noexcept {
    swap(rhs);
    return *this;
  }

  std::size_t size() const { return m_vertices.size(); }
  std::size_t capacity() const { return m_vertices.capacity(); }

  std::size_t num_edges() const { return m_num_edges; }

  vertex_stored_type& get_stored_vertex(vertex_descriptor v) const {
    return BC_get_value(*BC_desc_to_iterator(m_vertices, v));
  }

  const edge_stored_type& get_stored_edge(const edge_descriptor& e) const {
    return BC_get_value(
        *BC_desc_to_iterator(get_stored_vertex(e.source).out_edges, e.edge_id));
  }

  std::size_t get_out_degree(vertex_descriptor v) const {
    return BC_get_size(get_stored_vertex(v).out_edges);
  }

  std::size_t get_in_degree(vertex_descriptor v) const {
    return BC_get_size(get_stored_vertex(v).in_edges);
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  vertex_descriptor add_vertex(const VertexProperties& vp) {
    return adjlistBC_add_vertex(m_vertices, vp);
  }
  vertex_descriptor add_vertex(VertexProperties&& vp) {
    return adjlistBC_add_vertex(m_vertices, std::move(vp));
  }

  // NOTE: this operation only invalidates existing vertex-iterators,
  // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vecBC).
  // NOTE: This WORKS for ALL vertex container types.
  void clear_vertex(vertex_descriptor v) {
    adjlistBC_clear_vertex<DirectedS>(m_vertices, get_stored_vertex(v), v,
                                      m_num_edges);
  }

  // NOTE: this operation only invalidates existing vertex-descriptors,
  // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vecBC).
  // NOTE: This WORKS for ALL vertex container types.
  void remove_vertex(vertex_descriptor v) {
    clear_vertex(v);
    adjlistBC_erase_vertex(m_vertices, v);
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u,
                                            vertex_descriptor v,
                                            const EdgeProperties& ep) {
    using RawEDesc = typename edge_descriptor::edge_id_type;

    std::pair<RawEDesc, bool> raw_result =
        adjlistBC_add_edge(get_stored_vertex(u).out_edges, ep, v);

    if (raw_result.second) {
      ++m_num_edges;
      adjlistBC_add_in_edge(get_stored_vertex(v),
                            edge_descriptor(u, raw_result.first));
      return std::pair<edge_descriptor, bool>(
          edge_descriptor(u, raw_result.first), true);
    }
    return std::pair<edge_descriptor, bool>(edge_descriptor(), false);
  }
  std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u,
                                            vertex_descriptor v,
                                            EdgeProperties&& ep) {
    using RawEDesc = typename edge_descriptor::edge_id_type;

    std::pair<RawEDesc, bool> raw_result =
        adjlistBC_add_edge(get_stored_vertex(u).out_edges, std::move(ep), v);

    if (raw_result.second) {
      ++m_num_edges;
      adjlistBC_add_in_edge(get_stored_vertex(v),
                            edge_descriptor(u, raw_result.first));
      return std::pair<edge_descriptor, bool>(
          edge_descriptor(u, raw_result.first), true);
    }
    return std::pair<edge_descriptor, bool>(edge_descriptor(), false);
  }

  // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
  // NOTE: This WORKS for ALL vertex container types.
  void remove_edge(const edge_descriptor& e) {
    adjlistBC_erase_in_edge(get_stored_vertex(get_stored_edge(e).target), e);
    adjlistBC_erase_edge(get_stored_vertex(e.source).out_edges, e, m_vertices,
                         e.source);
    --m_num_edges;
  }

  // NOTE: This WORKS for ALL vertex container types.
  void clear() {
    using OEFactory = adjlistBC_out_edges_factory<
        typename vertex_value_type::edge_container_ptr>;
    OEFactory::destroy_all_out_edges(m_vertices);
    m_vertices.clear();
    m_num_edges = 0;
  }

  // NOTE: This WORKS for ALL vertex container types.
  using VIterSelect =
      adjlistBC_select_vertex_iterator<VertexListS, vertex_container>;
  using vertex_iterator = typename VIterSelect::type;

  std::pair<vertex_iterator, vertex_iterator> vertices() const {
    return VIterSelect::create_range(m_vertices);
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  using OEIterSelect =
      adjlistBC_select_out_edge_iterator<OutEdgeListS, edge_container,
                                         edge_descriptor>;
  using out_edge_iterator = typename OEIterSelect::type;

  std::pair<out_edge_iterator, out_edge_iterator> out_edges(
      vertex_descriptor u) const {
    return OEIterSelect::create_range(u, get_stored_vertex(u).out_edges);
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  std::pair<edge_descriptor, bool> get_edge(vertex_descriptor u,
                                            vertex_descriptor v) const {
    using RawEDesc = typename edge_descriptor::edge_id_type;

    std::pair<RawEDesc, bool> raw_result =
        adjlistBC_find_edge_to(get_stored_vertex(u).out_edges, v);

    if (raw_result.second) {
      return std::pair<edge_descriptor, bool>(
          edge_descriptor(u, raw_result.first), true);
    }
    return std::pair<edge_descriptor, bool>(edge_descriptor(), false);
  }

  // NOTE: This WORKS for ALL vertex container types.
  using in_edge_iterator = typename vertex_stored_type::in_edge_iterator;

  std::pair<in_edge_iterator, in_edge_iterator> in_edges(
      vertex_descriptor v) const {
    return std::pair<in_edge_iterator, in_edge_iterator>(
        get_stored_vertex(v).in_edges.begin(),
        get_stored_vertex(v).in_edges.end());
  }

  // NOTE: This WORKS for ALL vertex container types.
  using edge_iterator =
      adjlistBC_edge_iterator<vertex_container, edge_descriptor>;

  std::pair<edge_iterator, edge_iterator> edges() const {
    return std::pair<edge_iterator, edge_iterator>(
        edge_iterator::begin(m_vertices), edge_iterator::end(m_vertices));
  }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
void swap(adjlistBC_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& lhs,
          adjlistBC_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

}  // namespace boost::graph::detail

#endif
