// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Triangle_Dunavant.hpp"
#include "Hydrofem_FEQuadrature_Dunavant.hpp"

namespace hydrofem
{

namespace FEQuadrature_Dunavant
{

double CalculateIntegralOverTriangle(const int Order,
                                     const SPoint& point1, 
                                     const SPoint& point2, 
                                     const SPoint& point3, 
                                     const std::function<double(double,double)>& fun)
{
  using namespace TriangleDunavant;
  // quadrature rule
  auto rule = dunavant_degree(Order);
  // number of quadrature points
  auto order = dunavant_order_num(rule);
  // quadrature weights
  double* qw = new double[order];
  // quadrature points reference triangle
  double* qp_ref = new double[2*order];
  // compute the weights and quadrature points
  dunavant_rule(rule,order,qp_ref,qw);
  // current element coords
  double tri[2*3] = {point1.x(), point1.y(), point2.x(), point2.y(), point3.x(), point3.y()};
  // physical quadrature points
  double* qp = new double[2*order];
  // get the quadrature points at the current triangle   
  reference_to_physical_t3(tri,order,qp_ref,qp);
  // get the triangle area
  double area = triangle_area(tri);
  // integrate the function
  double val = 0.0;
  for (int i = 0; i < order; ++i)
    val += qw[i]*fun(qp[2*i+0],qp[2*i+1]);
  val *= area;
  
  delete [] qw; qw = nullptr;
  delete [] qp_ref; qp_ref = nullptr;
  delete [] qp; qp = nullptr;
  
  return val;
}

double CalculateIntegralOverTriangle(const int Order,
                                     const SPoint& point1,
                                     const SPoint& point2,
                                     const SPoint& point3,
                                     const std::vector<std::function<double(double,double)>>& funs,
                                     const std::vector<double>& multipliers)
{
  using namespace TriangleDunavant;
  // quadrature rule
  auto rule = dunavant_degree(Order);
  // number of quadrature points
  auto order = dunavant_order_num(rule);
  // quadrature weights
  double* qw = new double[order];
  // quadrature points reference triangle
  double* qp_ref = new double[2*order];
  // compute the weights and quadrature points
  dunavant_rule(rule,order,qp_ref,qw);
  // current element coords
  double tri[2*3] = {point1.x(), point1.y(), point2.x(), point2.y(), point3.x(), point3.y()};
  // physical quadrature points
  double* qp = new double[2*order];
  // get the quadrature points at the current triangle   
  reference_to_physical_t3(tri,order,qp_ref,qp);
  // get the triangle area
  double area = triangle_area(tri);
  // integrate the function
  double val = 0.0;
  for (int i = 0; i < order; ++i)
  {
    double prod = 1.0;
    for (int fInd = 0; fInd < int(funs.size()); ++fInd)
      prod *= funs[fInd](qp[2*i+0],qp[2*i+1]);
    for (int mInd = 0; mInd < int(multipliers.size()); ++mInd)
      prod *= multipliers[mInd];
    val += qw[i]*prod;
  }
  val *= area;
  
  delete [] qw; qw = 0;
  delete [] qp_ref; qp_ref = 0;
  delete [] qp; qp = 0;
  
  return val;  
}

}
// end namespace FEQuadrature_Dunavant

}
// end namespace hydrofem
