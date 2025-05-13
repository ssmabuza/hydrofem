// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Element_HPP__
#define __Hydrofem_Element_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_MeshTypes.hpp"

namespace hydrofem
{

/** \brief Element base class */
class Element
{
public:
  
  //! \brief Ctor
  Element() {}
  
  //! \brief copy Ctor
  Element(const Element& element)
  {
    m_element_type = element.m_element_type;
    m_nodes = element.m_nodes;
    m_edges = element.m_edges;
    area = element.area;
    orientation = element.orientation;
    m_domain_ind = element.m_domain_ind;
  }
  
  //! \brief Dtor
  virtual ~Element() {}
  
  //! \brief gets the type of element built
  virtual inline ObjType getElementType() const 
  { return m_element_type; }

  //! \brief
  virtual inline const std::vector<int>& getNodes() const 
  { return m_nodes; }
  
  //! \brief
  virtual inline std::vector<int>& setNodes()
  { return m_nodes; }

  //! \brief
  virtual inline const std::vector<int>& getEdges() const 
  { return m_edges; }
  
  //! \brief
  virtual inline std::vector<int>& setEdges()
  { return m_edges; }
  
  //! \brief
  virtual inline double getArea() const
  { return area; }
  
  //! \brief
  virtual inline double& setArea()
  { return area; }
    
  //! \brief
  virtual inline int getOrientation() const
  { return orientation; }
    
  //! \brief
  virtual inline int& setOrientation()
  { return orientation; }

  //! \brief gets the local edge index given a global index
  virtual inline int getIndexOfTheEdge(const int /*edge_ind*/) const
  { return 0; }
  
  //! \brief
  virtual inline bool operator==(const Element& element)
  {
    if (m_element_type != element.m_element_type)
      return false;
    if (orientation != element.orientation)
      return false;
    if (area != element.area)
      return false;
    if (m_nodes.size() != element.m_nodes.size())
      return false;
    for (std::size_t i (0); i < m_nodes.size(); ++i)
      if (m_nodes[i] != element.m_nodes[i])
        return false;
    if (m_edges.size() != element.m_edges.size())
      return false;
    for (std::size_t i (0); i < m_nodes.size(); ++i)
      if (m_edges[i] != element.m_edges[i])
        return false;
    return true;
  }

  //! \brief element type
  ObjType m_element_type;
  //! \brief vertex indexes
  std::vector<int> m_nodes;
  //! \brief edges 
  std::vector<int> m_edges;
  //! \brief element area 
  double area;
  //! \brief nodes oriented either clockwise (orientation = -1) or counterclockwise (orientation = +1)
  int orientation;
  //! \brief Domain index 0 by default
  int m_domain_ind = 0;
    
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Element_HPP__ */
