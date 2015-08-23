// Copyright Sven Mikael Persson 2013.
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

/**
 * \file bst_inorder_iterator.hpp
 * 
 * This library implements a simple binary search tree in-order iterator for a BGL tree.
 * 
 * \author Sven Mikael Persson <mikael.s.persson@gmail.com>
 * \date August 2013
 */

#ifndef BOOST_BST_INORDER_ITERATOR_HPP
#define BOOST_BST_INORDER_ITERATOR_HPP

#include <boost/bind.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/property_map/property_map.hpp>

#include <boost/graph/tree_concepts.hpp>
#include <boost/graph/tree_adaptor.hpp>

namespace boost {

namespace graph {
  
  
namespace detail {
  
  
  enum bst_traversal_status {
    OnLeftBranch,
    OnMiddleBranch,
    OnRightBranch
  };
  
  
  template <typename TreeType, typename VertexType>
  VertexType bst_go_down_left(const TreeType& aTree, VertexType aStart) {
    // look down-left until we reach the leaf:
    while(out_degree(aStart, aTree))
      aStart = *(child_vertices(aStart, aTree).first);
    return aStart;
  };
  
  template <typename TreeType, typename VertexType>
  VertexType bst_go_down_right(const TreeType& aTree, VertexType aStart) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    // look down-right until we reach the leaf:
    while(out_degree(aStart, aTree) > 1) {
      std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(aStart, aTree);
      ++(cur_children.first);
      aStart = *(cur_children.first);
    };
    return aStart;
  };
  
  template <typename TreeType, typename VertexType>
  void bst_move_up_to_next(const TreeType& aTree, VertexType& aU, bst_traversal_status& aStatus) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    typedef typename graph_traits<TreeType>::in_edge_iterator in_edge_iter;
    aStatus = OnRightBranch;
    while(true) {
      in_edge_iter ei, ei_end;
      tie(ei, ei_end) = in_edges(aU, aTree);
      if(ei == ei_end) // at the root, go to the end (root, rightbranch)
        return;
      VertexType v = source(*ei, aTree);
      child_vertex_iter vil, vi_end;
      tie(vil, vi_end) = child_vertices(v, aTree);
      if(*vil == aU) { // u is the left child of v.
        aU = v;
        aStatus = OnMiddleBranch;
        return;
      }
      // u must be the right child of v. keep going.
      aU = v;
    };
  };
  
  template <typename TreeType, typename VertexType>
  void bst_move_up_to_prev(const TreeType& aTree, VertexType& aU, bst_traversal_status& aStatus) {
    typedef typename tree_traits<TreeType>::child_vertex_iterator child_vertex_iter;
    typedef typename graph_traits<TreeType>::in_edge_iterator in_edge_iter;
    child_vertex_iter vil, vi_end;
    aStatus = OnLeftBranch;
    VertexType u = aU;
    while(true) {
      in_edge_iter ei, ei_end;
      tie(ei, ei_end) = in_edges(u, aTree);
      if(ei == ei_end) // at the root, so, aU must be the beginning node.
        return;
      VertexType v = source(*ei, aTree);
      child_vertex_iter vil, vi_end;
      tie(vil, vi_end) = child_vertices(v, aTree);
      if(*vil == u) { // u is the left child of v.
        u = v;
        continue;
      };
      // u must be the right child of v. keep going.
      aU = v;
      aStatus = OnMiddleBranch;
      return;
    };
  };
  
};  // detail


}; // graph



template <typename CompleteBinaryTree, typename ValueType>
class bst_inorder_iterator {
  public:
    typedef bst_inorder_iterator<CompleteBinaryTree, ValueType> self;
    typedef CompleteBinaryTree tree_type;
    
    typedef std::ptrdiff_t difference_type;
    typedef ValueType value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef std::bidirectional_iterator_tag iterator_category;
    
  private:
    typedef typename graph_traits<tree_type>::vertex_descriptor vertex_type;
    typedef typename tree_traits<tree_type>::child_vertex_iterator child_vertex_iter;
    typedef typename graph_traits<tree_type>::in_edge_iterator in_edge_iter;
    
    
    tree_type* m_tree;
    vertex_type m_u;
    graph::detail::bst_traversal_status m_status;
    
    bst_inorder_iterator(tree_type* aTree, vertex_type aU, graph::detail::bst_traversal_status aStatus) : m_tree(aTree), m_u(aU), m_status(aStatus) { };
    
  public:
    
    vertex_type base() const { return m_u; };
    
    bst_inorder_iterator(tree_type* aTree, vertex_type aU) : m_tree(aTree), m_u(aU), m_status(graph::detail::OnRightBranch) {
      // must figure out what the case is.
      if((!m_tree) || (m_u == graph_traits<tree_type>::null_vertex()))
        return;
      if(m_u == get_root_vertex(*m_tree)) {
        m_status = graph::detail::OnMiddleBranch;
        return;
      };
      // first check if there are any children:
      if(out_degree(m_u, *m_tree)) {
        m_status = graph::detail::OnMiddleBranch; // not on leaf.
        return;
      };
      // then, check if m_u is a left or right child of its parent:
      in_edge_iter ei, ei_end;
      tie(ei, ei_end) = in_edges(m_u, *m_tree);
      vertex_type v = source(*ei, *m_tree);
      child_vertex_iter vil, vi_end;
      tie(vil, vi_end) = child_vertices(v, *m_tree);
      if(*vil == m_u) // u is the left child of v.
        m_status = graph::detail::OnLeftBranch;
      else // u must be the right child of v.
        m_status = graph::detail::OnRightBranch;
    };
    
    static self begin(tree_type* aTree) {
      if((aTree) && (num_vertices(*aTree)))
        return self(aTree, graph::detail::bst_go_down_left(*aTree, get_root_vertex(*aTree)), graph::detail::OnLeftBranch);
      else
        return self(NULL, graph_traits<tree_type>::null_vertex(), graph::detail::OnRightBranch);
    };
    
    static self end(tree_type* aTree) {
      if((aTree) && (num_vertices(*aTree)))
        return self(aTree, get_root_vertex(*aTree), graph::detail::OnRightBranch);
      else
        return self(NULL, graph_traits<tree_type>::null_vertex(), graph::detail::OnRightBranch);
    };
    
    friend bool operator==( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree == rhs.m_tree) && (lhs.m_u == rhs.m_u) && (lhs.m_status == rhs.m_status));
    };
    friend bool operator!=( const self& lhs, const self& rhs) { 
      return ((lhs.m_tree != rhs.m_tree) || (lhs.m_u != rhs.m_u) || (lhs.m_status != rhs.m_status));
    };
    
    self& operator++() { 
      if(!m_tree)
        return *this;
      switch(m_status) {
        case graph::detail::OnLeftBranch:
          // on a left-leaf, must move up to the parent as a middle-point
          {
            std::pair<in_edge_iter, in_edge_iter> eis = in_edges(m_u, *m_tree);
            if(eis.first == eis.second) { // at the root, go to the end (root, rightbranch).
              m_status = graph::detail::OnRightBranch;
              return *this;
            };
            m_u = source(*(eis.first), *m_tree);
            m_status = graph::detail::OnMiddleBranch;
          };
          break;
        case graph::detail::OnMiddleBranch:
          // on a middle-point, must move down to the right once, and then left to the bottom.
          {
            if( out_degree(m_u, *m_tree) > 1 ) {
              // go to the right child.
              std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
              ++(cur_children.first);
              m_u = graph::detail::bst_go_down_left(*m_tree, *(cur_children.first));
              if(m_u == *(cur_children.first))
                m_status = graph::detail::OnRightBranch;
              else
                m_status = graph::detail::OnLeftBranch;
              break;
            };
          };
          // this means that we must move up to the next value (no right child here).
        case graph::detail::OnRightBranch:
          graph::detail::bst_move_up_to_next(*m_tree, m_u, m_status);
          break;
      };
      return *this;
    };
    
    self& operator--() { 
      if(!m_tree)
        return *this;
      switch(m_status) {
        case graph::detail::OnRightBranch:
          if(m_u == get_root_vertex(*m_tree)) {
            // go to the left child, and down the right:
            m_u = graph::detail::bst_go_down_right(*m_tree, m_u);
            if(m_u == get_root_vertex(*m_tree))
              m_status = graph::detail::OnMiddleBranch;
            else {
              if(out_degree(m_u, *m_tree) == 0) 
                m_status = graph::detail::OnRightBranch;  // on the leaf.
              else
                m_status = graph::detail::OnMiddleBranch; // not on leaf.
            };
            break;
          };
          // this means that we are either on a right-leaf or on a mis-labeled middle-node, 
          // in either case, try the middle-branch case:
        case graph::detail::OnMiddleBranch:
          // on a middle-point or right-point, must move down to the left once (if exists), and then right to the bottom.
          {
            if(out_degree(m_u, *m_tree)) {
              // go to the left child, and down the right:
              std::pair<child_vertex_iter, child_vertex_iter> cur_children = child_vertices(m_u, *m_tree);
              m_u = graph::detail::bst_go_down_right(*m_tree, *(cur_children.first));
              if(m_u == *(cur_children.first)) {
                if(out_degree(m_u, *m_tree) == 0)
                  m_status = graph::detail::OnLeftBranch;
                else
                  m_status = graph::detail::OnMiddleBranch;
              } else {
                if(out_degree(m_u, *m_tree) == 0) 
                  m_status = graph::detail::OnRightBranch;  // on the leaf.
                else
                  m_status = graph::detail::OnMiddleBranch; // not on leaf.
              };
              break;
            };
          };
          // this means that we must move up to the previous value (no left child here).
        case graph::detail::OnLeftBranch:
          graph::detail::bst_move_up_to_prev(*m_tree, m_u, m_status);
          break;
      };
      return *this;
    };
    
    self operator++(int) { self result(*this); return ++result; };
    self operator--(int) { self result(*this); return --result; };
    
    reference operator*() { return (*m_tree)[m_u]; };
    pointer operator->() { return &(*m_tree)[m_u]; };
};




};


#endif









