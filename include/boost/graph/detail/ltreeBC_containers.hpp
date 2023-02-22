// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file ltreeBC_containers.hpp
 * 
 * This library implements adjacency-list data structures that underpin the adjacency_list_v2 template. 
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date June 2013
 */

#ifndef BOOST_LTREEBC_CONTAINERS_HPP
#define BOOST_LTREEBC_CONTAINERS_HPP

#include <boost/config.hpp>

#include <boost/graph/detail/adjlistBC_containers.hpp>
#include <boost/graph/detail/adjlistBC_iterators.hpp>
#include <boost/graph/detail/boost_container_generators.hpp>

#include <boost/graph/graph_selectors.hpp>  // for directedS, undirectedS, bidirectionalS.

#include <iterator>
#include <queue>
#include <stack>
#include <utility>
#include <vector>

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
struct ltreeBC_vertex_stored_type;  // forward declaration.

// this is fine because boost-containers allow incomplete types:
template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct ltreeBC_vertex_config {
  using stored_type =
      ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>;
  using container = typename BC_container_gen<VertexListS, stored_type>::type;
  using value_type = typename container::value_type;
  using descriptor = typename BC_select_descriptor<container>::type;

  static descriptor null_vertex() { return BC_null_desc<descriptor>::value(); }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct ltreeBC_edge_stored_type {
  using self = ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                        VertexProperties, EdgeProperties>;
  using VConfig = ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                        VertexProperties, EdgeProperties>;
  using vertex_descriptor = typename VConfig::descriptor;

  vertex_descriptor target;
  mutable EdgeProperties data;

  explicit ltreeBC_edge_stored_type(vertex_descriptor aTarget)
      : target(aTarget), data() {}
  ltreeBC_edge_stored_type(vertex_descriptor aTarget,
                           const EdgeProperties& aData)
      : target(aTarget), data(aData) {}
  ltreeBC_edge_stored_type()
      : ltreeBC_edge_stored_type(VConfig::null_vertex()) {}
  ltreeBC_edge_stored_type(vertex_descriptor aTarget, EdgeProperties&& aData)
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
void swap(ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>& lhs,
          ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>& rhs) {
  using std::swap;
  swap(lhs.target, rhs.target);
  swap(lhs.data, rhs.data);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
std::size_t hash_value(
    const ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>& ep) {
  return BC_desc_get_hash(ep.target);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct ltreeBC_edge_config {
  using vertex_descriptor =
      typename ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties,
                                     EdgeProperties>::descriptor;

  using stored_type =
      ltreeBC_edge_stored_type<VertexListS, OutEdgeListS, DirectedS,
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
struct ltreeBC_vertex_stored_type {
  using self = ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                          VertexProperties, EdgeProperties>;
  using directed_tag = DirectedS;
  using Config = ltreeBC_edge_config<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_container_ptr = typename Config::container_ptr;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = edge_descriptor*;

  VertexProperties data;
  edge_container_ptr out_edges;
  edge_descriptor in_edge;

  ltreeBC_vertex_stored_type() : data(), out_edges(), in_edge() {}
  explicit ltreeBC_vertex_stored_type(const VertexProperties& aData)
      : data(aData), out_edges(), in_edge() {}
  explicit ltreeBC_vertex_stored_type(VertexProperties&& aData)
      : data(std::move(aData)), out_edges(), in_edge() {}
  void swap(self& rhs) {
    using std::swap;
    swap(data, rhs.data);
    swap(out_edges, rhs.out_edges);
    swap(in_edge, rhs.in_edge);
  }
  self& operator=(self rhs) {
    this->swap(rhs);
    return *this;
  }
};

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties>
struct ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                  VertexProperties, EdgeProperties> {
  using self = ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                                          VertexProperties, EdgeProperties>;
  using directed_tag = directedS;
  using Config = ltreeBC_edge_config<VertexListS, OutEdgeListS, directedS,
                                     VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_container_ptr = typename Config::container_ptr;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = int*;

  VertexProperties data;
  edge_container_ptr out_edges;

  ltreeBC_vertex_stored_type() : data(), out_edges() {}
  explicit ltreeBC_vertex_stored_type(const VertexProperties& aData)
      : data(aData), out_edges() {}
  explicit ltreeBC_vertex_stored_type(VertexProperties&& aData)
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
void swap(ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& lhs,
          ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                     VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

//NOTE: ltreeBC_out_edges_factory == adjlistBC_out_edges_factory
//NOTE: ltreeBC_add_edge = adjlistBC_add_edge
//NOTE: ltreeBC_find_edge_to = adjlistBC_find_edge_to
//NOTE: ltreeBC_add_vertex = adjlistBC_add_vertex

/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC    (multi)setBC      unordered_(multi)setBC
// directedS:   O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)    O(log(E/V))       O(1)

//NOTE: ltreeBC_erase_edges_to = adjlistBC_erase_edges_to
// As long as these special versions of update_in_edge_id are used (by ADL):

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlistBC_update_in_edge_id(
    ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                               VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlistBC_update_in_edge_id(
    ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                               VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {
  using VertexValue =
      ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                                 VertexProperties, EdgeProperties>;
  using EdgeDesc = typename VertexValue::edge_descriptor;

  if (vp.in_edge == EdgeDesc::null_value()) {
    return;
  }
  if ((vp.in_edge.source == v) && (vp.in_edge.edge_id == old_id)) {
    vp.in_edge.edge_id = new_id;
  }
}

/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))

// for OutEdgeListS = listS, setS, ...
template <typename Container, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void ltreeBC_erase_edge(Container& cont, EdgeDesc e, VertexCont& vcont,
                        VertexDesc v) {
  adjlistBC_erase_edge(cont, e, vcont, v);
}
template <typename Container, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void ltreeBC_erase_edge(Container* cont, EdgeDesc e, VertexCont& vcont,
                        VertexDesc v) {
  adjlistBC_erase_edge(*cont, e, vcont, v);
}

// for OutEdgeListS = vecS
template <typename ValueType, typename EdgeDesc, typename VertexCont,
          typename VertexDesc>
void ltreeBC_erase_edge(::boost::container::vector<ValueType>& cont, EdgeDesc e,
                        VertexCont& vertex_cont, VertexDesc v) {
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

/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

template <typename VertexListS, typename OutEdgeListS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void ltreeBC_add_in_edge(
    ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, directedS,
                               VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties, typename EdgeDesc>
void ltreeBC_add_in_edge(
    ltreeBC_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS,
                               VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {
  vp.in_edge = e;
}

/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
// NOTE: The function 'ltreeBC_erase_vertex' works for all graph types.
// NOTE: Time complexities:
//               vecBC      poolBC      listBC
// directedS:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

template <typename Container, typename VDescContainer>
void ltreeBC_erase_vertices(Container& cont, VDescContainer& v_list) {
  for (typename VDescContainer::iterator it = v_list.begin();
       it != v_list.end(); ++it) {
    adjlistBC_erase_vertex(cont, *it);
  }
}

//NOTE: ltreeBC_update_out_edges_impl = adjlistBC_update_out_edges_impl

// O(E)
template <typename DirectedS, typename ValueType>
typename enable_if<is_same<DirectedS, directedS>, void>::type
ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont,
                         std::size_t old_v_id, std::size_t new_v_id) {
  adjlistBC_update_out_edges<DirectedS>(cont, old_v_id, new_v_id);
}

// O(E/V)
template <typename DirectedS, typename ValueType>
typename disable_if<is_same<DirectedS, directedS>, void>::type
ltreeBC_update_out_edges(::boost::container::vector<ValueType>& cont,
                         std::size_t old_v_id, std::size_t new_v_id) {
  using OutEdgeCont = typename ValueType::edge_container;
  using OEIter = typename OutEdgeCont::iterator;
  using EdgeDesc = typename ValueType::edge_descriptor;
  // first, update in-edge vertices
  if (cont[new_v_id].in_edge != EdgeDesc::null_value()) {
    EdgeDesc& e_in = cont[new_v_id].in_edge;
    ValueType& up = cont[e_in.source];
    // std::cout << "address 3: " << cont[new_v_id].data << " " << up.data << " " << e_in.edge_id << " " << up.out_edges.size() << std::endl;
    adjlistBC_update_out_edges_impl(
        up.out_edges, old_v_id, new_v_id,
        BC_desc_to_iterator(up.out_edges, e_in.edge_id));
  }
  // second, update out-edge vertices
  for (OEIter ei = BC_get_begin_iter(cont[new_v_id].out_edges);
       ei != BC_get_end_iter(cont[new_v_id].out_edges); ++ei) {
    if (!BC_is_elem_valid(*ei)) {
      continue;
    }
    ValueType& wp = cont[BC_get_value(*ei).target];
    EdgeDesc& e_in = wp.in_edge;
    if (e_in != EdgeDesc::null_value()) {
      if ((e_in.source == old_v_id) &&
          (ei == BC_desc_to_iterator(cont[new_v_id].out_edges, e_in.edge_id))) {
        e_in.source = new_v_id;
        break;
      }
    }
  }
}

// O(1)
template <typename DirectedS, typename ValueType>
typename enable_if<is_same<DirectedS, directedS>, void>::type
ltreeBC_invalidate_in_edges(::boost::container::vector<ValueType>& /*unused*/,
                            std::size_t /*unused*/) {}

// O(1)
template <typename DirectedS, typename ValueType>
typename disable_if<is_same<DirectedS, directedS>, void>::type
ltreeBC_invalidate_in_edges(::boost::container::vector<ValueType>& cont,
                            std::size_t u) {
  using OutEdgeCont = typename ValueType::edge_container;
  using OEIter = typename OutEdgeCont::iterator;
  using EdgeDesc = typename ValueType::edge_descriptor;

  for (OEIter ei = BC_get_begin_iter(cont[u].out_edges);
       ei != BC_get_end_iter(cont[u].out_edges); ++ei) {
    if (!BC_is_elem_valid(*ei)) {
      continue;
    }
    cont[BC_get_value(*ei).target].in_edge = EdgeDesc::null_value();
  }
}

template <typename ValueType>
void ltreeBC_erase_vertices(::boost::container::vector<ValueType>& cont,
                            std::vector<std::size_t>& v_list) {
  using Container = ::boost::container::vector<ValueType>;
  using Iter = typename Container::iterator;
  using OutEdgeFactory =
      adjlistBC_out_edges_factory<typename ValueType::edge_container_ptr>;
  using DirectedS = typename ValueType::directed_tag;
  using std::swap;

  for (auto v_it = v_list.begin(); v_it != v_list.end(); ++v_it) {
    // NOTE: First, invalidate the in-edge references of the children, then destroy out-edge list.
    ltreeBC_invalidate_in_edges<DirectedS>(cont, *v_it);
    Iter it = BC_desc_to_iterator(cont, *v_it);
    OutEdgeFactory::destroy_out_edges(*it);
    Iter it_last = cont.end();
    --it_last;
    if (it == it_last) {
      cont.erase(it_last);
      continue;
    }
    swap(*it, *it_last);
    std::size_t old_id = it_last - cont.begin();
    std::size_t new_id = it - cont.begin();
    cont.erase(it_last);
    ltreeBC_update_out_edges<DirectedS>(cont, old_id, new_id);
    for (auto v_it2 = v_it + 1; v_it2 != v_list.end(); ++v_it2) {
      if (*v_it2 == old_id) {
        *v_it2 = new_id;
      }
    }
  }
}

template <typename ValueType>
void ltreeBC_erase_vertices(BC_pooled_vector<ValueType>& cont,
                            std::vector<std::size_t>& v_list) {
  for (unsigned long& it : v_list) {
    adjlistBC_erase_vertex(cont, it);
  }
}

// O(E)
template <typename DirectedS, typename Container, typename VertexValue,
          typename Vertex>
typename enable_if<is_same<DirectedS, directedS>,
                   VertexValue>::type::edge_descriptor
ltreeBC_get_in_edge(Container& vcont, VertexValue& vp, Vertex v) {
  using VIter = typename Container::iterator;
  using EdgeCont = typename VertexValue::edge_container;
  using OEIter = typename EdgeCont::iterator;
  using EdgeDesc = typename VertexValue::edge_descriptor;

  for (VIter it = vcont.begin(); it != vcont.end(); ++it) {
    if (BC_is_elem_valid(*it)) {
      VertexValue& up = BC_get_value(*it);
      for (OEIter ei = BC_get_begin_iter(up.out_edges);
           ei != BC_get_end_iter(up.out_edges); ++ei) {
        if (BC_is_elem_valid(*ei) && (BC_get_value(*ei).target == v)) {
          return EdgeDesc(BC_iterator_to_desc(vcont, it),
                          BC_iterator_to_desc(up.out_edges, ei));
        }
      }
    }
  }

  return EdgeDesc::null_value();
}

// O(1)
template <typename DirectedS, typename Container, typename VertexValue,
          typename Vertex>
typename disable_if<is_same<DirectedS, directedS>,
                    VertexValue>::type::edge_descriptor
ltreeBC_get_in_edge(Container& vcont, VertexValue& vp, Vertex v) {
  return vp.in_edge;
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
struct ltreeBC_vertex_container {

  using VConfig = ltreeBC_vertex_config<VertexListS, OutEdgeListS, DirectedS,
                                        VertexProperties, EdgeProperties>;
  using vertex_container = typename VConfig::container;
  using vertices_size_type = typename vertex_container::size_type;
  using vertex_descriptor = typename VConfig::descriptor;
  using vertex_stored_type = typename VConfig::stored_type;
  using vertex_value_type = typename VConfig::value_type;

  using EConfig = ltreeBC_edge_config<VertexListS, OutEdgeListS, DirectedS,
                                      VertexProperties, EdgeProperties>;
  using edge_container = typename EConfig::container;
  using edges_size_type = typename edge_container::size_type;
  using edge_descriptor = typename EConfig::descriptor;
  using edge_stored_type = typename EConfig::stored_type;
  using edge_value_type = typename EConfig::value_type;

  mutable vertex_container m_vertices;
  vertex_descriptor m_root;

  ltreeBC_vertex_container()
      : m_vertices(), m_root(BC_null_desc<vertex_descriptor>::value()) {}
  ~ltreeBC_vertex_container() { clear(); }

  ltreeBC_vertex_container(const ltreeBC_vertex_container&) = delete;
  ltreeBC_vertex_container& operator=(const ltreeBC_vertex_container&) = delete;

  void swap(ltreeBC_vertex_container& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    // swap(m_vertices, rhs.m_vertices);
    swap(m_root, rhs.m_root);
  }

  ltreeBC_vertex_container(ltreeBC_vertex_container&& rhs) noexcept
      : m_vertices(), m_root(BC_null_desc<vertex_descriptor>::value()) {
    swap(rhs);
  }
  ltreeBC_vertex_container& operator=(ltreeBC_vertex_container&& rhs) noexcept {
    swap(rhs);
    return *this;
  }

  /********************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE SAME AS ADJ-LIST-BC - NOTE NOTE NOTE *
 * ******************************************************************************/

  std::size_t size() const { return m_vertices.size(); }
  std::size_t capacity() const { return m_vertices.capacity(); }

  std::size_t num_edges() const {
    return (m_vertices.size() == 0 ? 0 : m_vertices.size() - 1);
  }

  vertex_stored_type& get_stored_vertex(vertex_descriptor v) const {
    return BC_get_value(*BC_desc_to_iterator(m_vertices, v));
  }

  const edge_stored_type& get_stored_edge(edge_descriptor e) const {
    return BC_get_value(
        *BC_desc_to_iterator(get_stored_vertex(e.source).out_edges, e.edge_id));
  }

  std::size_t get_out_degree(vertex_descriptor v) const {
    return BC_get_size(get_stored_vertex(v).out_edges);
  }

  std::size_t get_in_degree(vertex_descriptor v) const {
    if ((v != BC_null_desc<vertex_descriptor>::value()) && (v != m_root) &&
        (BC_is_elem_valid(*BC_desc_to_iterator(m_vertices, v)))) {
      return 1;
    }
    return 0;
  }

  // NOTE: This WORKS for ALL vertex container types.
  void clear() {
    using OEFactory = adjlistBC_out_edges_factory<
        typename vertex_value_type::edge_container_ptr>;
    OEFactory::destroy_all_out_edges(m_vertices);
    m_vertices.clear();
    m_root = BC_null_desc<vertex_descriptor>::value();
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
  using edge_iterator =
      adjlistBC_edge_iterator<vertex_container, edge_descriptor>;

  std::pair<edge_iterator, edge_iterator> edges() const {
    return std::pair<edge_iterator, edge_iterator>(
        edge_iterator::begin(m_vertices), edge_iterator::end(m_vertices));
  }

  /***************************************************************************************
 * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE DIFFERENT FROM ADJ-LIST-BC - NOTE NOTE NOTE *
 * *************************************************************************************/

  // NOTE: This WORKS for ALL vertex container types.
  using in_edge_iterator = typename vertex_stored_type::in_edge_iterator;

  std::pair<in_edge_iterator, in_edge_iterator> in_edges(
      vertex_descriptor v) const {
    if (get_stored_vertex(v).in_edge == EConfig::null_edge()) {
      return std::pair<in_edge_iterator, in_edge_iterator>(NULL, NULL);
    }
    return std::pair<in_edge_iterator, in_edge_iterator>(
        &(get_stored_vertex(v).in_edge), &(get_stored_vertex(v).in_edge) + 1);
  }

  // NOTE: This WORKS for ALL vertex container types.
  vertex_descriptor get_parent(vertex_descriptor v) const {
    return ltreeBC_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v)
        .source;
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  std::size_t get_depth(const vertex_descriptor& u) const {
    using OEIter = typename edge_container::iterator;
    using TaskType = std::pair<std::size_t, vertex_descriptor>;

    std::size_t max_depth = 0;
    std::stack<TaskType> tasks;
    tasks.push(TaskType(0, u));
    while (!tasks.empty()) {
      TaskType cur = tasks.top();
      tasks.pop();
      ++(cur.first);
      if (cur.first > max_depth) {
        max_depth = cur.first;
      }
      vertex_stored_type& vp = get_stored_vertex(cur.second);
      for (OEIter ei = BC_get_begin_iter(vp.out_edges);
           ei != BC_get_end_iter(vp.out_edges); ++ei) {
        tasks.push(TaskType(cur.first, BC_get_value(*ei).target));
      }
    }
    return max_depth;
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
  // NOTE: This WORKS for ALL edge container types.
  template <typename VP>
  vertex_descriptor add_new_vertex(VP&& vp) {
    return adjlistBC_add_vertex(m_vertices, std::forward<VP>(vp));
  }
  template <typename VP>
  void add_root_vertex(VP&& vp) {
    m_root = add_new_vertex(std::forward<VP>(vp));
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EP>
  std::pair<edge_descriptor, bool> add_new_edge(vertex_descriptor u,
                                                vertex_descriptor v, EP&& ep) {
    using RawEDesc = typename edge_descriptor::edge_id_type;
    std::pair<RawEDesc, bool> raw_result = adjlistBC_add_edge(
        get_stored_vertex(u).out_edges, std::forward<EP>(ep), v);
    if (raw_result.second) {
      ltreeBC_add_in_edge(get_stored_vertex(v),
                          edge_descriptor(u, raw_result.first));
      return std::pair<edge_descriptor, bool>(
          edge_descriptor(u, raw_result.first), true);
    }
    return std::pair<edge_descriptor, bool>(edge_descriptor(), false);
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename VP, typename EP>
  std::pair<vertex_descriptor, edge_descriptor> add_child(vertex_descriptor v,
                                                          VP&& vp, EP&& ep) {
    vertex_descriptor new_node = add_new_vertex(std::forward<VP>(vp));
    std::pair<edge_descriptor, bool> new_edge =
        add_new_edge(v, new_node, std::forward<EP>(ep));
    return std::make_pair(new_node, new_edge.first);
  }

  template <typename Vertex_OIter, typename Edge_OIter>
  std::pair<Vertex_OIter, Edge_OIter> clear_children_impl(vertex_descriptor v,
                                                          Vertex_OIter vit_out,
                                                          Edge_OIter eit_out) {
    using OEIter = typename edge_container::iterator;

    std::vector<vertex_descriptor> death_row;
    std::queue<vertex_descriptor> bft_queue;
    bft_queue.push(v);
    // Put all children on death-row:
    while (!bft_queue.empty()) {
      vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
      bft_queue.pop();
      for (OEIter ei = BC_get_begin_iter(v_value.out_edges);
           ei != BC_get_end_iter(v_value.out_edges); ++ei) {
        if (!BC_is_elem_valid(*ei)) {
          continue;
        }
        death_row.push_back(BC_get_value(*ei).target);
        bft_queue.push(BC_get_value(*ei).target);
        vertex_stored_type& u_value =
            get_stored_vertex(BC_get_value(*ei).target);
        *(vit_out++) = std::move(u_value.data);
        *(eit_out++) = std::move(BC_get_value(*ei).data);
      }
    }

    // Check if we removed the root:
    if (v == m_root) {  // v is the root vertex.
      clear();
    } else {
      // remove the out-edges.
      BC_clear_all(get_stored_vertex(v).out_edges);
      // Execute the death-row vertices!
      ltreeBC_erase_vertices(m_vertices, death_row);
    }

    return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
  }

  template <typename OutputIter>
  OutputIter clear_children_impl(vertex_descriptor v, OutputIter it_out) {
    return clear_children_impl(v, it_out, ignore_output_iter()).first;
  }

  void clear_children_impl(vertex_descriptor v) {
    clear_children_impl(v, ignore_output_iter(), ignore_output_iter());
  }

  template <typename Vertex_OIter, typename Edge_OIter>
  std::pair<Vertex_OIter, Edge_OIter> remove_branch_impl(vertex_descriptor v,
                                                         Vertex_OIter vit_out,
                                                         Edge_OIter eit_out) {
    using OEIter = typename edge_container::iterator;

    edge_descriptor p_edge =
        ltreeBC_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v);

    if (p_edge.source == VConfig::null_vertex()) {
      *(eit_out++) = edge_value_type();
    } else {
      *(eit_out++) = std::move(get_stored_edge(p_edge).data);
    }

    std::vector<vertex_descriptor> death_row;
    death_row.push_back(v);
    std::queue<vertex_descriptor> bft_queue;
    bft_queue.push(v);
    // Put all children on death-row:
    while (!bft_queue.empty()) {
      vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
      *(vit_out++) = std::move(v_value.data);
      bft_queue.pop();
      for (OEIter ei = BC_get_begin_iter(v_value.out_edges);
           ei != BC_get_end_iter(v_value.out_edges); ++ei) {
        if (!BC_is_elem_valid(*ei)) {
          continue;
        }
        death_row.push_back(BC_get_value(*ei).target);
        bft_queue.push(BC_get_value(*ei).target);
        *(eit_out++) = std::move(BC_get_value(*ei).data);
      }
    }

    // Check if we removed the root:
    if (p_edge.source == VConfig::null_vertex()) {
      // v must be the root vertex.
      clear();
      return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
    }  // remove the edge.
    ltreeBC_erase_edge(get_stored_vertex(p_edge.source).out_edges, p_edge,
                       m_vertices, p_edge.source);

    // Execute the death-row vertices!
    ltreeBC_erase_vertices(m_vertices, death_row);

    return std::pair<Vertex_OIter, Edge_OIter>(vit_out, eit_out);
  }

  template <typename OutputIter>
  OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
    return remove_branch_impl(v, it_out, ignore_output_iter()).first;
  }

  void remove_branch_impl(vertex_descriptor v) {
    remove_branch_impl(v, ignore_output_iter(), ignore_output_iter());
  }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS,
          typename VertexProperties, typename EdgeProperties>
void swap(ltreeBC_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>& lhs,
          ltreeBC_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                   VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

}  // namespace boost::graph::detail

#endif
