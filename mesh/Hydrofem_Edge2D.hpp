// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Edge2D_HPP__
#define __Hydrofem_Edge2D_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Edge.hpp"

namespace hydrofem
{

/** \brief Edges in 2D */
class Edge2D
  :
  public Edge
{
public:

  //! \brief Ctor
  Edge2D() : Edge(-1,-1) { m_elems.resize(2,-1); }

  //! \brief Ctor
  Edge2D(const int node1, const int node2) : Edge(node1,node2)
  { m_elems.resize(2,-1); }

  //! \brief Copy Ctor
  Edge2D(const Edge2D& edge) : Edge(edge.m_nodes.at(0),edge.m_nodes.at(1))
  {
    m_elems.resize(2);
    m_elems.at(0) = edge.m_elems.at(0);
    m_elems.at(1) = edge.m_elems.at(1);
    m_is_boundary = edge.m_is_boundary;
    length        = edge.length;
  }
  
  //! \brief Dtor
  virtual ~Edge2D() {}
  
  //! \brief Copy assignment
  Edge2D& operator=(const Edge2D& edge)
  {
    m_nodes       = edge.m_nodes;
    m_elems       = edge.m_elems;
    m_is_boundary = edge.m_is_boundary;
    length        = edge.length;
    return *this;
  }
  
  //! \brief equality check
  bool operator==(const Edge2D& edge)
  {
    
    return (m_nodes.at(0)==edge.m_nodes.at(0)) &&
           (m_nodes.at(1)==edge.m_nodes.at(1)) &&
           (m_elems.at(0)==edge.m_elems.at(0)) &&
           (m_elems.at(1)==edge.m_elems.at(1)) &&
           (m_is_boundary==edge.m_is_boundary) &&
           (length==edge.length);
    
  }

  /////////////////////////////////////////////////////////////////////////////////////////
  //!   \brief function getSignofEdgeInElement
  //!   function return the sign of edge in the given element, it is assumed that
  //!   the direction of the edge is from element with smaller index into the
  //!   element with larger index, thus is function called with element with
  //!   smaller index the result is 1.0, otherwise -1.0
  //!   input parameters : elem_ind - index of element
  //!   Remark : the functin terminates in case the edge does not belong to element
  /////////////////////////////////////////////////////////////////////////////////////////
  double getSignofEdgeInElement(const int elem_ind) const;
  
};
/// end class Edge2D


class SearchEdgeStruct
{
public:
  
  SearchEdgeStruct() {}
  SearchEdgeStruct(const int node1, const int node2, const int edge_ind);
  bool isTheSameEdge(const int node1, const int node2) const;
  int m_edge_ind;
  int m_nodes[2];
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Edge2D_HPP__ */
