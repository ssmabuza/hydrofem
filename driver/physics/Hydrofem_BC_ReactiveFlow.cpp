// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_BC_ReactiveFlow.hpp"

namespace hydrofem
{

void BC_ReactiveFlow::initialize()
{
  assert(m_mesh);
  assert(m_velocity);
  
  for (int elem_ind(0); elem_ind < m_mesh->numOfElements(); ++elem_ind)
  {
    for (int ledge_ind(0); ledge_ind < int(m_mesh->getElement(elem_ind).m_edges.size()); ++ledge_ind)
    {
      const int gedge_ind = m_mesh->getElement(elem_ind).m_edges[ledge_ind];
      if (m_mesh->getEdge(gedge_ind).m_is_boundary)
      {
        m_bc_info_edges[gedge_ind] = std::make_unique<BCInfo>();
        
        const auto gedge_vertices = m_mesh->getEdgeVertices(gedge_ind);
        const auto edge_centre = 0.5*(gedge_vertices[0]+gedge_vertices[1]);
        const auto normall = m_mesh->evalEdgeNormal(elem_ind,ledge_ind);
        const double v_dot_n = (*m_velocity)(edge_centre)*normall;
        
        if (v_dot_n > 0)
          m_bc_info_edges[gedge_ind]->m_boundary_condition_type = typeGammaPlus;
        else if (v_dot_n < 0)
          m_bc_info_edges[gedge_ind]->m_boundary_condition_type = typeGammaMinus;
        else if (edge_centre.x() > 0.01)
          m_bc_info_edges[gedge_ind]->m_boundary_condition_type = typeSigma;
        else  
          m_bc_info_edges[gedge_ind]->m_boundary_condition_type = typeGammaZero;
      }
    }
  }
  
  // update index in system for each side/boundary type
  int index = 0;
  for (std::unordered_map<int,std::unique_ptr<BCInfo>>::iterator it = m_bc_info_edges.begin(); it != m_bc_info_edges.end(); ++it)
  {
    if (it->second->m_boundary_condition_type == typeGammaMinus)
    {
      it->second->m_index_in_system = index;
      ++index;
    }
  }
  
  index = 0;
  for (std::unordered_map<int,std::unique_ptr<BCInfo>>::iterator it = m_bc_info_edges.begin(); it != m_bc_info_edges.end(); ++it)
  {
    if (it->second->m_boundary_condition_type == typeGammaPlus)
    {
      it->second->m_index_in_system = index;
      ++index;
    }
  }
  
  index = 0;
  for (std::unordered_map<int,std::unique_ptr<BCInfo>>::iterator it = m_bc_info_edges.begin(); it != m_bc_info_edges.end(); ++it)
  {
    if (it->second->m_boundary_condition_type == typeGammaZero)
    {
      it->second->m_index_in_system = index;
      ++index;
    }
  }

    index = 0;
  for (std::unordered_map<int,std::unique_ptr<BCInfo>>::iterator it = m_bc_info_edges.begin(); it != m_bc_info_edges.end(); ++it)
  {
    if (it->second->m_boundary_condition_type == typeSigma)
    {
      it->second->m_index_in_system = index;
      ++index;
    }
  }
  

}

}
// end namespace hydrofem
