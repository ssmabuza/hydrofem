// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_BC_Scalar.hpp"

namespace hydrofem
{

void BC_Scalar::initialize()
{
  /**
   * Nothing to do since by default, we have no boundary calculations.
   * 
   */
}

void BC_Scalar::initializeBoundaryPointsToDirichletEverywhere()
{
  for (int elem_ind(0); elem_ind < m_mesh->numOfElements(); ++elem_ind)
  {
    for (int ledge_ind(0); ledge_ind < int(m_mesh->getElement(elem_ind).m_edges.size()); ++ledge_ind)
    {
      const int gedge_ind = m_mesh->getElement(elem_ind).m_edges[ledge_ind];
      if (m_mesh->getEdge(gedge_ind).m_is_boundary)
      {
        const int &i = m_mesh->getEdge(gedge_ind).m_nodes[0];
        const int &j = m_mesh->getEdge(gedge_ind).m_nodes[1];
        if (m_bc_info_points.find(i)==m_bc_info_points.end())
        {
          m_bc_info_points[i] = std::make_shared<BCInfo>();
          m_bc_info_points[i]->m_boundary_condition_type = typeDirichlet;
        }
        if (m_bc_info_points.find(j)==m_bc_info_points.end())
        {
          m_bc_info_points[j] = std::make_shared<BCInfo>();
          m_bc_info_points[j]->m_boundary_condition_type = typeDirichlet;
        }
      }
    }
  }

}

}
// end namespace hydrofem
