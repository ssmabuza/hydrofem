// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Quadrature_Tri.hpp"

#include "Hydrofem_Triangle_Dunavant.hpp"
#include "Hydrofem_FEQuadrature_Dunavant.hpp"

namespace hydrofem
{

Quadrature_Tri::Quadrature_Tri(const int order, const std::vector<SPoint>& tri)
{
  buildQuadrature(order,tri);
}

void Quadrature_Tri::buildQuadrature(const int order, const std::vector<SPoint>& tri)
{
  auto& _q_points = get_q_pointsView();
  auto& _q_weights = get_q_weightsView();
  
  using namespace TriangleDunavant;
  // quadrature rule
  auto rule = dunavant_degree(order);
  // number of quadrature points
  auto order_ = dunavant_order_num(rule);
  // quadrature weights
  double* qw = new double[order_];
  // quadrature points reference triangle
  double* qp_ref = new double[2*order_];
  // compute the weights and quadrature points
  dunavant_rule(rule,order_,qp_ref,qw);
  // current element coords
  std::vector<double> tri_ = {tri[0].x(), tri[0].y(), tri[1].x(), tri[1].y(), tri[2].x(), tri[2].y()};
  // physical quadrature points
  double* qp = new double[2*order_];
  // get the quadrature points at the current triangle   
  reference_to_physical_t3(tri_.data(),order_,qp_ref,qp);
  // get the triangle area
  double area = triangle_area(tri_.data());

  for (int i = 0; i < order_; ++i)
    qw[i] = qw[i]*area;

  _q_points.clear();
  _q_weights.clear();
  for (int i = 0; i < order_; ++i)
  {
    _q_points.push_back(SPoint(qp[2*i+0],qp[2*i+1]));
    _q_weights.push_back(qw[i]);
  }

  delete [] qp;
  delete [] qw;
  delete [] qp_ref;
}

}
// end namespace hydrofem
