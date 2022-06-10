// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_HGrad_DofMapper_Line.hpp"

namespace hydrofem
{

int HGrad_DofMapper_Line::eval_global(const int i_elem, const int k) const
{
  int i_elem_glob = m_mesh->globalElemInd(i_elem);
  assert((k >= 0) && (k <= m_p) && (i_elem_glob >= 0) && (i_elem_glob < m_nelems));
  return i_elem_glob * m_p  + k;
}
  
int HGrad_DofMapper_Line::local(const int k) const
{
  assert((k >= 0) && (k <= m_p));
  return k;
}

}
// end namespace hydrofem
