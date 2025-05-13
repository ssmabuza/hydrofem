// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_HGrad_DofMapper_Triangle.hpp"

namespace hydrofem
{

int HGrad_DofMapper_Triangle::local(const int i, const int j, const int /*k*/) const
{
  return i * (int(2 * m_p + 3) - i) / int(2) + j;
}
  
int HGrad_DofMapper_Triangle::eval_global(const int i_elem, const int i, const int j, const int k) const
{
  const int i_elem_glob = m_mesh->globalElemInd(i_elem);
  assert(i_elem_glob < m_mesh->numOfElements());
  assert((0 <= i) && (0 <= j) && (0 <= k) && (m_p == i + j + k));
  
  //! get the global element index from this proc's local element index
  const int elemInd = i_elem_glob;
  //! get the global point indexes on the local element
  const auto & pointInd = m_mesh->getElement(i_elem).m_nodes;
  //! get the global edge indexes on each element given the local element index
  const auto & edgeInd = m_mesh->getElement(i_elem).m_edges;

  // degree of freedom corresponds to a node value
  if (i == m_p)
  {
    return int(pointInd.at(0));
  }
  else if (j == m_p)
  {
    return int(pointInd.at(1));
  }
  else if (k == m_p)
  {
    return int(pointInd.at(2));
  }
  // degree of freedom lies on an edge value
  else if (i == 0)
  {
    if (pointInd.at(1) < pointInd.at(2))
    {
      return m_npoints + edgeInd.at(0) * (m_p - 1) + j - 1;
    }
    else
    {
      return m_npoints + edgeInd.at(0) * (m_p - 1) + k - 1;
    }
  }
  else if (j == 0)
  {
    if (pointInd.at(0) < pointInd.at(2))
    {
      return m_npoints + edgeInd.at(1) * (m_p - 1) + i - 1;
    }
    else
    {
      return m_npoints + edgeInd.at(1) * (m_p - 1) + k - 1;
    }
  }
  else if (k == 0)
  {
    if (pointInd.at(0) < pointInd.at(1))
    {
      return m_npoints + edgeInd.at(2) * (m_p - 1) + i - 1;
    }
    else
    {
      return m_npoints + edgeInd.at(2) * (m_p - 1) + j - 1;
    }
  }
  // inner degree of freedom
  else
  {
    return m_npoints + m_nedges * (m_p - 1)
          + elemInd * ((m_p - 2) * (m_p - 1)) / 2
          + (i + j - 1) * (i + j - 2) / 2 + j - 1;
  }
}

}
// end namespace hydrofem
