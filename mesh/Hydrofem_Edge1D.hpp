// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Edge1D_HPP__
#define __Hydrofem_Edge1D_HPP__

#include "Hydrofem_Edge.hpp"

namespace hydrofem
{

/** \brief An edge in 1D */
class Edge1D
  :
  public Edge
{
public:
  
  //! \brief Ctor
  Edge1D() : Edge(-1,-1) { m_elems.resize(2,-1); }
  
  //! \brief Ctor
  Edge1D(const int node1, const int node2) : Edge(node1,node2)
  { m_elems.resize(2,-1); }

  //! \brief Copy Ctor
  Edge1D(const Edge1D& edge) : Edge(edge.m_nodes.at(0),edge.m_nodes.at(1))
  {
    m_elems.resize(2);
    m_elems.at(0) = edge.m_elems.at(0);
    m_elems.at(1) = edge.m_elems.at(1);
    m_is_boundary = edge.m_is_boundary;
    length        = edge.length;
  }
  
  //! \brief Dtor
  virtual ~Edge1D() {}
  
  //! \brief Copy assignment
  Edge1D& operator=(const Edge1D& edge)
  {
    m_nodes       = edge.m_nodes;
    m_elems       = edge.m_elems;
    m_is_boundary = edge.m_is_boundary;
    length        = edge.length;
    return *this;
  }
  
  //! \brief equality check
  bool operator==(const Edge1D& edge)
  {
    
    return (m_nodes.at(0)==edge.m_nodes.at(0)) &&
           (m_nodes.at(1)==edge.m_nodes.at(1)) &&
           (m_elems.at(0)==edge.m_elems.at(0)) &&
           (m_elems.at(1)==edge.m_elems.at(1)) &&
           (m_is_boundary==edge.m_is_boundary) &&
           (length==edge.length);
    
  }
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Edge1D_HPP__ */

