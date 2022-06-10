// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Element1D_HPP__
#define __Hydrofem_Element1D_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Element.hpp"
#include "Hydrofem_MeshTypes.hpp"

namespace hydrofem
{

/** \brief A 1D finite element */
class Element1D
  :
  public Element
{
  
public:
  
  //! \brief Default Ctor
  Element1D()
  {
    m_element_type = ObjType::typeLine;
    m_nodes.resize(2);
    m_nodes[0] = -1;
    m_nodes[1] = -1;
    m_edges.resize(1,-1);
    area = -1.0;
    orientation = 0;
  }
             
  //! Ctor from given node indicies
  Element1D(const int node1,
            const int node2)
  {
    m_element_type = ObjType::typeLine;
    m_nodes.resize(2);
    m_nodes[0] = node1;
    m_nodes[1] = node2;
    m_edges.resize(1,-1);
    area = -1.0;
    orientation = 0;
  }
  
  //! Ctor from given node indicies and points
  Element1D(const int node1,
            const int node2,
            const SPoint& point1,
            const SPoint& point2)
  {
    m_element_type = ObjType::typeLine;
    m_nodes.resize(2);
    m_nodes[0] = node1;
    m_nodes[1] = node2;
    m_edges = m_nodes;
    area = std::fabs(point1.x()-point2.x());
    orientation = 0;
  }
  
  //! Copy Ctor
  Element1D(const Element1D& elem)
  {
    m_element_type = elem.getElementType();
    m_nodes = elem.m_nodes;
    m_edges = elem.m_edges;
    area = elem.getArea();
    orientation = elem.orientation;
  }
  
  virtual ~Element1D() {}  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Element1D_HPP__ */
