// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Edge_HPP__
#define __Hydrofem_Edge_HPP__

#include "Hydrofem.hpp"

namespace hydrofem
{

/**
 * @brief An edge connecting two nodes (in 1D edge==element)
 */
class Edge
{
public:
  
  /** @brief Ctor */
  Edge(int node1, int node2) : m_nodes(std::vector<int>(2))
  {
    m_nodes[0] = node1;
    m_nodes[1] = node2;
    m_is_boundary = false;
    length = -1.0;
  }
  
  /** @brief Copy Ctor */
  Edge(const Edge& edge)
    :
    m_nodes(edge.m_nodes),
    m_elems(edge.m_elems),
    m_is_boundary(edge.m_is_boundary),
    length(edge.length) {}

  /** @brief Dtor */
  virtual ~Edge() {}
  
  /** @brief nodes in global mesh connected by the edge */
  std::vector<int> m_nodes;
  
  /** @brief elements connected to edge */
  std::vector<int> m_elems;

  /** @brief boundary flag */
  bool m_is_boundary;
  
  /** @brief length of edge */
  double length;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Edge_HPP__ */

