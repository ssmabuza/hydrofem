// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Edge2D.hpp"

namespace hydrofem
{

double Edge2D::getSignofEdgeInElement(const int elem_ind) const
{
  if (elem_ind != m_elems.at(0) && elem_ind != m_elems.at(1))
  {
    std::cout << "Error in Edge2D::getSignofFaceInElement : incorrect elem_ind : " << elem_ind << std::endl;
    exit(1);
  }
  
  if (m_elems.at(1) == -1)
    return 1.0;
  
  if (elem_ind == m_elems.at(0))
  {
    if (m_elems.at(0) < m_elems.at(1))
      return 1.0;
    else
      return -1.0;
  } else {
    if (m_elems.at(1) < m_elems.at(0))
      return 1.0;
    else
      return -1.0;
  }
}
/// end double Edge2D::getSignofEdgeInElement

SearchEdgeStruct::SearchEdgeStruct(const int node1, const int node2, const int edge_ind) : m_edge_ind(edge_ind)
{
  m_nodes[0] = node1;
  m_nodes[1] = node2;
}
/// end constructor SearchEdgeStruct::SearchEdgeStruct

bool SearchEdgeStruct::isTheSameEdge(const int node1, const int node2) const
{
  if ((node1 == m_nodes[0] && node2 == m_nodes[1]) || (node1 == m_nodes[1] && node2 == m_nodes[0]))
    return true;
  return false;
}
/// end bool SearchEdgeStruct::isTheSameEdge

}
// end namespace hydrofem

