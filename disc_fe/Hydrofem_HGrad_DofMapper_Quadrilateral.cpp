// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_HGrad_DofMapper_Quadrilateral.hpp"


namespace hydrofem
{

int HGrad_DofMapper_Quadrilateral::
eval_global(const int i_elem, const int i, const int j) const
{
  // first get the global index
  const int i_elem_glob = m_mesh->globalElemInd(i_elem);
  assert( (i >= 0) && (i <= m_p) && (j >= 0) && (j <= m_p) && (i_elem_glob < m_nelems) && (i_elem_glob >= 0));
  int index=-1;
  const std::vector<int>& pointInd = m_mesh->getElement(i_elem).m_nodes;
  const std::vector<int>& edgeInd = m_mesh->getElement(i_elem).m_edges;
  if ( (i == 0) && (j == 0) )
  {
    // on point 1:
    index = pointInd.at(0);
  } else if ( (i == 0) && (j == m_p) ) {
    // on point 2:
    index = pointInd.at(3);
  } else if ( (i == m_p) && (j == 0) ) {
    // on point 3:
    index = pointInd.at(1);
  } else if ( (i == m_p) && (j == m_p) ) {
    // on point 4:
    index = pointInd.at(2);
  } else if ( (i > 0) && (i < m_p) && ((j == 0) || (j == m_p)) ) {
    
    // on edge 1 or 3:
    int g_edge_ind = -1;
    if (j == m_p)
      // edge 3
      g_edge_ind = edgeInd.at(2);
    else if (j == 0)
      // edge 1
      g_edge_ind = edgeInd.at(0);
    
    const Edge& edge = m_mesh->getEdge(g_edge_ind);
    if (edge.m_nodes.at(0) < edge.m_nodes.at(1))
    {
      index = nnodes() + g_edge_ind * (m_p - 1) + i;
    } else {
      index = nnodes() + g_edge_ind * (m_p - 1) + ( (m_p - 1) - i );
    }
    
  } else if ( ((i == 0) || (i == m_p)) && (j > 0) && (j < m_p) ) {
    
    // on edge 2 or 4:
    int g_edge_ind = -1;
    if (i == 0)
      // edge 4
      g_edge_ind = edgeInd.at(3);
    else if (i == m_p)
      // edge 2
      g_edge_ind = edgeInd.at(1);
    assert(g_edge_ind>=0);
    
    const Edge& edge = m_mesh->getEdge(g_edge_ind);
    if (edge.m_nodes.at(0) < edge.m_nodes.at(1))
    {
      index = nnodes() + g_edge_ind * (m_p - 1) + j;
    } else {
      index = nnodes() + g_edge_ind * (m_p - 1) + ( (m_p - 1) - j );
    }
    
  } else {
    
    // element interior
    index = nnodes() + nedges()*(m_p - 1) + (i - 1)*(m_p - 1) + j;
    
  }
  
  assert(index >= 0);
  return index;
}

int HGrad_DofMapper_Quadrilateral::local(const int i, const int j) const
{
  assert((i >= 0) && (i <= m_p) && (j >= 0) && (j <= m_p));
  return i*(m_p+1)+j;
}

}
// end namespace hydrofem
