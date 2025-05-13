// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Quadrature_Quad.hpp"
#include "Hydrofem_FEQuadrature_Quadrilateral.hpp"

namespace hydrofem
{

using namespace FEQuadrature_Quadrilateral;

Quadrature_Quad::Quadrature_Quad(const int order, const std::vector<SPoint>& quad)
{
  buildQuadrature(order,quad);
}

void Quadrature_Quad::buildQuadrature(const int order, const std::vector<SPoint>& quad)
{
  auto& _q_points = get_q_pointsView();
  auto& _q_weights = get_q_weightsView();
  // reference quadrature points (x_i,y_j)
  std::vector<SPoint> ref_qp;
  // quadrature weights w_iw_j
  std::vector<double> qweights;
  // get weights and nodes
  GetQuadrature(ref_qp,qweights,order);
  // map the reference quadrature points to current element
  MapPhysicalRectangleToReferenceRectangle(_q_points,ref_qp,quad.at(0),quad.at(1),quad.at(2),quad.at(3));
  // get the weights
  _q_weights.resize(_q_points.size());
  for (int qpInd = 0; qpInd < int(_q_points.size()); ++qpInd)
  {
    const double detJ = JacobianDeterminant(ref_qp[qpInd],quad.at(0),quad.at(1),quad.at(2),quad.at(3));
    _q_weights[qpInd] = qweights[qpInd]*detJ;
  }
}

}
// end namespace hydrofem
