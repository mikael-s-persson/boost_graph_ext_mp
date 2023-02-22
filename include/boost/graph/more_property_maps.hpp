// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file more_property_maps.hpp
 *
 * This library provides additional property-maps that are especially handy when writing custom graph classes.
 *
 * \author Sven Mikael Persson
 * \date May 2013
 */

#ifndef BOOST_MORE_PROPERTY_MAPS_HPP
#define BOOST_MORE_PROPERTY_MAPS_HPP

#include <boost/graph/properties.hpp>
#include <boost/property_map/property_map.hpp>

namespace boost {

/****************************************************************
 *       Whole-bundle property-map
 * **************************************************************/

/**
 * This property-map delivers the entire vertex or edge bundle associated to a
 * vertex / edge descriptor of a graph.
 * \note This property-map relies of the existence of an operator[] for the graph and 
 * for the relevant descriptor type. This operator[] should deliver the vertex/edge bundle.
 * \tparam Graph The graph type.
 * \tparam PropertyMapTag The property tag type that identifies if it is a vertex_bundle_t or edge_bundle_t.
 */
template <typename Graph, typename PropertyMapTag>
struct whole_bundle_property_map
    : public put_get_helper<
          typename mpl::if_<is_same<PropertyMapTag, vertex_bundle_t>,
                            typename Graph::vertex_bundled,
                            typename Graph::edge_bundled>::type&,
          whole_bundle_property_map<Graph, PropertyMapTag>> {
 private:
  Graph* pg;

 public:
  using is_vertex_bundle = is_same<PropertyMapTag, vertex_bundle_t>;
  using is_const_graph = is_const<Graph>;
  using value_type =
      typename mpl::if_<is_vertex_bundle, typename Graph::vertex_bundled,
                        typename Graph::edge_bundled>::type;
  using reference =
      typename mpl::if_<is_const_graph, const value_type&, value_type&>::type;
  using key_type =
      typename mpl::if_<is_vertex_bundle,
                        typename graph_traits<Graph>::vertex_descriptor,
                        typename graph_traits<Graph>::edge_descriptor>::type;
  using category = typename mpl::if_<is_const_graph, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit whole_bundle_property_map(Graph* aPG) : pg(aPG) {}
  whole_bundle_property_map() : whole_bundle_property_map(nullptr) {}
  reference operator[](key_type k) const { return (*pg)[k]; }
};

/****************************************************************
 *       Property-graph property-map
 * **************************************************************/

/**
 * This property-map uses a graph's "get" function to obtain the 
 * property value associated to a given descriptor.
 * \tparam T The value-type of the property.
 * \tparam Graph The graph type.
 * \tparam PropertyMapTag The tag associated to the property.
 */
template <typename T, typename Graph, typename PropertyMapTag>
struct tagged_from_bundle_property_map
    : public put_get_helper<
          T&, tagged_from_bundle_property_map<T, Graph, PropertyMapTag>> {
 private:
  Graph* pg;
  PropertyMapTag tag;

 public:
  using is_vertex_prop = is_same<typename property_kind<PropertyMapTag>::type,
                                 vertex_property_tag>;
  using is_const_graph = is_const<Graph>;
  using value_type = T;
  using reference = T&;
  using key_type =
      typename mpl::if_<is_vertex_prop,
                        typename graph_traits<Graph>::vertex_descriptor,
                        typename graph_traits<Graph>::edge_descriptor>::type;
  using category = typename mpl::if_<is_const<T>, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit tagged_from_bundle_property_map(
      Graph* aPG, PropertyMapTag aTag = PropertyMapTag())
      : pg(aPG), tag(aTag) {}
  tagged_from_bundle_property_map()
      : tagged_from_bundle_property_map(nullptr) {}
  reference operator[](key_type k) const {
    return get_property_value((*pg)[k], tag);
  }
};

/****************************************************************
 *       Property-graph property-map
 * **************************************************************/

/**
 * This property-map uses a graph's "get" function to obtain the 
 * property value associated to a given descriptor.
 * \tparam T The value-type of the property.
 * \tparam Graph The graph type.
 * \tparam PropertyMapTag The tag associated to the property.
 */
template <typename T, typename Graph, typename PropertyMapTag>
struct propgraph_property_map
    : public put_get_helper<T&,
                            propgraph_property_map<T, Graph, PropertyMapTag>> {
 private:
  Graph* pg;
  PropertyMapTag tag;

 public:
  using is_vertex_prop = is_same<typename property_kind<PropertyMapTag>::type,
                                 vertex_property_tag>;
  using is_const_graph = is_const<Graph>;
  using value_type = T;
  using reference = T&;
  using key_type =
      typename mpl::if_<is_vertex_prop,
                        typename graph_traits<Graph>::vertex_descriptor,
                        typename graph_traits<Graph>::edge_descriptor>::type;
  using category = typename mpl::if_<is_const<T>, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit propgraph_property_map(Graph* aPG,
                                  PropertyMapTag aTag = PropertyMapTag())
      : pg(aPG), tag(aTag) {}
  propgraph_property_map() : propgraph_property_map(nullptr) {}
  reference operator[](key_type k) const { return get(tag, *pg, k); }
};

/****************************************************************
 *       Bundle-data-member property-map
 * **************************************************************/

/**
 * This property-map delivers a data-member of the vertex or edge bundle associated to a
 * vertex / edge descriptor of a graph. This is similar to the property-map obtained by 
 * calling get(&SomeBundle::SomeMember, my_graph), and can be used to implement such a 
 * functionality for custom graph classes.
 * \note This property-map relies of the existence of an operator[] for the graph and 
 * for the relevant descriptor type. This operator[] should deliver the vertex/edge bundle.
 * \tparam T The value-type of the property (i.e., data-member of the bundle).
 * \tparam Graph The graph type.
 * \tparam PropertyMapTag The property tag type that identifies if it is a vertex_bundle_t or edge_bundle_t.
 */
template <typename T, typename Graph, typename PropertyMapTag>
class bundle_member_property_map
    : public put_get_helper<
          T&, bundle_member_property_map<T, Graph, PropertyMapTag>> {
 public:
  using self = bundle_member_property_map<T, Graph, PropertyMapTag>;
  using is_vertex_bundle = is_same<PropertyMapTag, vertex_bundle_t>;
  using is_const_graph = is_const<Graph>;
  using bundle_type =
      typename mpl::if_<is_vertex_bundle, typename Graph::vertex_bundled,
                        typename Graph::edge_bundled>::type;
  using member_ptr_type = T bundle_type::*;

 private:
  Graph* pg;
  member_ptr_type mem_ptr;

 public:
  using value_type = T;
  using reference = T&;
  using key_type =
      typename mpl::if_<is_vertex_bundle,
                        typename graph_traits<Graph>::vertex_descriptor,
                        typename graph_traits<Graph>::edge_descriptor>::type;
  using category = typename mpl::if_<is_const_graph, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit bundle_member_property_map(Graph* aPG,
                                      member_ptr_type aMemPtr = nullptr)
      : pg(aPG), mem_ptr(aMemPtr) {}
  bundle_member_property_map() : bundle_member_property_map(nullptr) {}
  reference operator[](key_type p) const { return (*pg)[p].*mem_ptr; }
};

/****************************************************************
 *       Sub-object put-get helper
 * **************************************************************/

/**
 * This is a CRTP-style base-class that can be used to imbue a property-map class 
 * with put and get functions. The property-map class is only required to provide 
 * an operator[].
 * \note This helper class applies to property-maps whose value is a subobject of the key object.
 */
template <typename Reference, typename LvaluePropertyMap>
struct subobject_put_get_helper {};

template <typename PropertyMap, typename Reference, typename K>
const typename property_traits<PropertyMap>::value_type& get(
    const subobject_put_get_helper<Reference, PropertyMap>& pa, const K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
}

template <typename PropertyMap, typename Reference, typename K>
Reference get(const subobject_put_get_helper<Reference, PropertyMap>& pa,
              K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
}

template <typename PropertyMap, typename Reference, typename K, typename V>
void put(const subobject_put_get_helper<Reference, PropertyMap>& pa, K& k,
         V&& v) {
  static_cast<const PropertyMap&>(pa)[k] = std::forward<V>(v);
}

/****************************************************************
 *       Self property-map
 * **************************************************************/

/**
 * This property-map is an "identity" map that preserves the reference 
 * semantics. In other words, it takes a key by reference and delivers 
 * it back unchanged (by reference).
 * \tparam T The value-type of the property (and the key).
 */
template <typename T>
struct self_property_map
    : public subobject_put_get_helper<T&, self_property_map<T>> {
  using self = self_property_map;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using key_type = T;
  using category = typename mpl::if_<is_const<T>, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  self_property_map() = default;
  ;
  reference operator[](reference p) const { return p; }
  const_reference operator[](const_reference p) const { return p; }
};

/****************************************************************
 *       Data-member property-map
 * **************************************************************/

/**
 * This property-map class can be used to map an object (by reference) to one 
 * of its data-members (by reference).
 * \tparam T The value type of the property (the data-member).
 * \tparam PropertyType The type of the object from which the data-member is picked.
 */
template <typename T, typename PropertyType>
class data_member_property_map
    : public subobject_put_get_helper<
          T&, data_member_property_map<T, PropertyType>> {
 public:
  using member_ptr_type = T PropertyType::*;
  using self = data_member_property_map<T, PropertyType>;

 private:
  member_ptr_type mem_ptr;

 public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using key_type = PropertyType;
  using category = typename mpl::if_<is_const<T>, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit data_member_property_map(member_ptr_type aMemPtr)
      : mem_ptr(aMemPtr) {}
  data_member_property_map() : data_member_property_map(nullptr) {}
  reference operator[](key_type& p) const { return p.*mem_ptr; }
  const_reference operator[](const key_type& p) const { return p.*mem_ptr; }
};

/**
 * This property-map class can be used to map an object (by reference) to one 
 * of its data-members (by reference).
 * \note This specialization applies to const objects (and const data-member).
 * \tparam T The value type of the property (the data-member).
 * \tparam PropertyType The type of the object from which the data-member is picked.
 */
template <typename T, typename PropertyType>
class data_member_property_map<const T, const PropertyType>
    : public subobject_put_get_helper<
          const T&, data_member_property_map<const T, const PropertyType>> {
 public:
  using value_type = T;
  using reference = const T&;
  using const_reference = const T&;
  using key_type = const PropertyType;

  using member_ptr_type = T key_type::*;
  using self = data_member_property_map<const T, const PropertyType>;

 private:
  member_ptr_type mem_ptr;

 public:
  using category = typename mpl::if_<is_const<T>, readable_property_map_tag,
                                     lvalue_property_map_tag>::type;

  explicit data_member_property_map(member_ptr_type aMemPtr = NULL)
      : mem_ptr(aMemPtr) {}
  data_member_property_map() : data_member_property_map(nullptr) {}
  reference operator[](key_type& p) const { return p.*mem_ptr; }
};

/****************************************************************
 *       Composite property-map
 * **************************************************************/

/**
 * This property-map allows for the composition of two property-maps.
 * The resulting property-map has the value-type of the output-map and 
 * the key-type of the input-map.
 * \tparam OutputMap The output map (or outer-map) whose key-type is the value-type of the input map.
 * \tparam InputMap The input map (or inner-map) whose value-type is the key-type of the output map.
 */
template <typename OutputMap, typename InputMap>
class composite_property_map {
 public:  // private:  would be private is friends were more portable.
  OutputMap prop_out;
  InputMap prop_in;

  using value_type = typename property_traits<OutputMap>::value_type;
  using key_type =
      typename remove_const<typename property_traits<InputMap>::key_type>::type;
  using category = typename property_traits<OutputMap>::category;
  using reference = typename property_traits<OutputMap>::reference;
  using const_reference = const reference;

  explicit composite_property_map(OutputMap aPropOut,
                                  InputMap aPropIn = InputMap())
      : prop_out(aPropOut), prop_in(aPropIn) {}
  composite_property_map() : composite_property_map(OutputMap()) {}

  reference operator[](const key_type& k) const { return prop_out[prop_in[k]]; }
  reference operator[](key_type& k) const { return prop_out[prop_in[k]]; }
};

template <typename OutputMap, typename InputMap>
typename composite_property_map<OutputMap, InputMap>::value_type get(
    const composite_property_map<OutputMap, InputMap>& m,
    const typename composite_property_map<OutputMap, InputMap>::key_type& p) {
  return m.prop_out[m.prop_in[p]];
}

template <typename OutputMap, typename InputMap, typename V>
void put(
    const composite_property_map<OutputMap, InputMap>& m,
    const typename composite_property_map<OutputMap, InputMap>::key_type& p,
    V&& value) {
  put(m.prop_out, m.prop_in[p], std::forward<V>(value));
}

template <typename OutputMap, typename InputMap, typename V>
void put(const composite_property_map<OutputMap, InputMap>& m,
         typename composite_property_map<OutputMap, InputMap>::key_type& p,
         V&& value) {
  put(m.prop_out, m.prop_in[p], std::forward<V>(value));
}

/**
 * This function template can be used to construct a property-map that takes a 
 * vertex-descriptor into a graph and maps it to a data member of the vertex-bundle.
 * This property-map is constructed from bundle-to-data-member map and a graph.
 * \tparam BundleMemberMap The property-map type that can map a bundle to one of its data-members.
 * \tparam Graph The graph type to which the vertex-bundle belongs.
 * \param bundle_prop The bundle-to-data-member property-map.
 * \param g The graph.
 * \return A property-map that can map a vertex-descriptor to a vertex-bundle data-member.
 */
template <typename BundleMemberMap, typename Graph>
composite_property_map<BundleMemberMap,
                       whole_bundle_property_map<Graph, vertex_bundle_t>>
bundle_prop_to_vertex_prop(BundleMemberMap bundle_prop, Graph& g) {
  using ResultType =
      composite_property_map<BundleMemberMap,
                             whole_bundle_property_map<Graph, vertex_bundle_t>>;
  return ResultType(bundle_prop,
                    whole_bundle_property_map<Graph, vertex_bundle_t>(&g));
}

/**
 * This function template can be used to construct a property-map that takes a 
 * edge-descriptor into a graph and maps it to a data member of the edge-bundle.
 * This property-map is constructed from bundle-to-data-member map and a graph.
 * \tparam BundleMemberMap The property-map type that can map a bundle to one of its data-members.
 * \tparam Graph The graph type to which the edge-bundle belongs.
 * \param bundle_prop The bundle-to-data-member property-map.
 * \param g The graph.
 * \return A property-map that can map a edge-descriptor to a edge-bundle data-member.
 */
template <typename BundleMemberMap, typename Graph>
composite_property_map<BundleMemberMap,
                       whole_bundle_property_map<Graph, edge_bundle_t>>
bundle_prop_to_edge_prop(BundleMemberMap bundle_prop, Graph& g) {
  using ResultType =
      composite_property_map<BundleMemberMap,
                             whole_bundle_property_map<Graph, edge_bundle_t>>;
  return ResultType(bundle_prop,
                    whole_bundle_property_map<Graph, edge_bundle_t>(&g));
}

}  // namespace boost

#endif
