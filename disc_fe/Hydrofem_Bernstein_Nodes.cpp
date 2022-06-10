// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Bernstein_Nodes.hpp"

#include "Hydrofem_ReferenceQuadrilateral.hpp"

#include "Hydrofem_HGrad_DofMapper_Line.hpp"
#include "Hydrofem_HGrad_DofMapper_Triangle.hpp"
#include "Hydrofem_HGrad_DofMapper_Quadrilateral.hpp"

namespace hydrofem
{

std::vector<SPoint>
Bernstein_Nodes::getNodes(const int elem_ind) const
{
  std::vector<SPoint> points;
  const auto nodes = m_dofmapper->mesh()->getElementVertices(elem_ind);
  const int m_p = m_dofmapper->p();
  if (m_dofmapper->mesh()->meshType()==MeshType::typeLineMesh)
  {
    const double hx = m_dofmapper->mesh()->getElementArea(elem_ind);
    for (int i (0); i <= m_p; ++i)
    {
      SPoint x(nodes[0].x() + i*hx);
      points.push_back(x);
    }
  }
  else if (m_dofmapper->mesh()->meshType()==MeshType::typeTriangularMesh)
  {
    for (int i(0); i <= int(m_p); ++i)
    {
      for (int k(m_p - i), j(0); j <= int(m_p - i); ++j, --k)
      {
        SPoint x((nodes[0].x() * double(i) + nodes[1].x() * double(j) + nodes[2].x() * double(k)) / double(m_p),
                 (nodes[0].y() * double(i) + nodes[1].y() * double(j) + nodes[2].y() * double(k)) / double(m_p));
        points.push_back(x);
      }
    }
  }
  else if (m_dofmapper->mesh()->meshType()==MeshType::typeQuadrilateralMesh)
  {
    double hx = 1.0/m_p, hy = 1.0/m_p;
    for (int i(0); i <= int(m_p); ++i)
    {
      for (int j(0); j <= int(m_p); ++j)
      {
        const auto ref_point = SPoint(i*hx,j*hy);
        points.push_back(referenceQuadToPhysical(nodes[0],nodes[1],nodes[2],nodes[3],ref_point));
      }
    }
  }  
  return points;
}

}
// end namespace hydrofem
