// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEQuadrature.hpp"
#include "Hydrofem_Quadrature_Point.hpp"

namespace hydrofem
{

Quadrature_Point::Quadrature_Point(const SPoint& point)
{
  assert(int(point.size()) == 1);
  buildQuadrature(point);
}

void Quadrature_Point::buildQuadrature(const SPoint& point)
{
  using namespace FEQuadrature;

  auto& _q_points = get_q_pointsView();
  auto& _q_weights = get_q_weightsView();
  
  _q_points.clear();
  _q_weights.clear();
  
  _q_points.push_back(point);
  _q_weights.push_back(1.0);
}

}
// end namespace hydrofem
