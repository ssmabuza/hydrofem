// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Element2D.hpp"

namespace hydrofem
{

Element2D::~Element2D() {}
/// end destructor Element2D::~Element2D()

int Element2D::
getIndexOfTheEdge(const int edge_ind) const
{
  for (std::vector<int>::size_type i = 0; i < m_edges.size(); ++i)
    if (m_edges[i] == edge_ind)
      return i;

  std::cout << "In Element2D::getIndexOfTheEdge : can not find required edge_ind: " << edge_ind << std::endl;
  exit(1);
}
/// end int Element2D::getIndexOfTheEdge

}
// end namespace hydrofem
