// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Element2D_HPP__
#define __Hydrofem_Element2D_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_MeshTypes.hpp"

#include "Hydrofem_Element.hpp"

namespace hydrofem
{

//! \brief Elements in 2D
class Element2D
  :
  public Element
{
public:
  
  //! \brief Ctor
  Element2D() : Element()
  {
    m_element_type = ObjType::typeGenericElement;
    area = -1.0;
    orientation = 0;
  }
  
  //! \brief Copy Ctor
  Element2D(const Element2D& elem) : Element()
  {
    m_element_type = elem.getElementType();
    m_nodes = elem.m_nodes;
    m_edges = elem.m_edges;
    area = elem.getArea();
    orientation = elem.orientation;
  }
  
  //! \brief Dtor
  virtual ~Element2D();
  
  /// Methods
  //@{ 
  
  /////////////////////////////////////////////////////////////////////////////////////////
  //!  \brief  function GetIndexOfTheEdge
  //!          function the index of the edge edge_ind in the local array of elmement edges
  //!          input parameters : edge_ind - index of edge ( global )
  //!    Remark : the functin terminates in case the edge does not belong to element
  /////////////////////////////////////////////////////////////////////////////////////////
  int getIndexOfTheEdge(const int edge_ind) const override;
  
  //@}
  
};
/// end class Element2D

//! \brief A triangular element
class Triangle
  :
  public Element2D
{
public:
  
  Triangle(const int node1, const int node2, const int node3) : Element2D()
  {
    m_element_type = ObjType::typeTriangle;
    m_nodes.resize(3);
    m_nodes[0] = node1;
    m_nodes[1] = node2;
    m_nodes[2] = node3;
    m_edges.resize(3,-1);
  }
  
};
/// end class Triangle

//! \brief A quadrilateral element
class Quadrilateral
  :
  public Element2D
{
public:
  
  Quadrilateral(const int node1, const int node2, const int node3, const int node4)
  {
    m_element_type = ObjType::typeQuadrilateral;
    m_nodes.resize(4);
    m_nodes[0] = node1;
    m_nodes[1] = node2;
    m_nodes[2] = node3;
    m_nodes[3] = node4;
    m_edges.resize(4,-1);
  }
  
};
/// end class Quadrilateral

}
// end namespace hydrofem

#endif /** __Hydrofem_Element2D_HPP__ */
