// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEQuadrature_Quadrilateral_HPP__
#define __Hydrofem_FEQuadrature_Quadrilateral_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

namespace FEQuadrature_Quadrilateral
{

// Isoparametric mapping of the points in the reference element to the physical element
// Reference element is : [-1,1]x[-1,1]
void MapPhysicalRectangleToReferenceRectangle(std::vector<SPoint>& phy_points,
                                              const std::vector<SPoint>& ref_points,
                                              const SPoint& phy_point1,
                                              const SPoint& phy_point2,
                                              const SPoint& phy_point3,
                                              const SPoint& phy_point4);

// Computes the quadrilateral area assuming points are order in counterclockwise form
double ComputeElementArea(const SPoint& point1,
                          const SPoint& point2,
                          const SPoint& point3,
                          const SPoint& point4);
                          
// Computes the Jacobian determinant at a point on the reference element
double JacobianDeterminant(const SPoint& ref_point,
                           const SPoint& point1,
                           const SPoint& point2,
                           const SPoint& point3,
                           const SPoint& point4);

// Rules for the quadrature weights and nodes in a reference element
void GetQuadrature(std::vector<SPoint>& qp,                                      
                   std::vector<double>& qweights,
                   const int rule_1d);
                                      
// Computes the integral over a quadrilateral shape
// Uses same number of nodes in the x and y direction
double CalculateIntegralOverRectangle(const int Rule, 
                                      const SPoint& point1,
                                      const SPoint& point2,
                                      const SPoint& point3,
                                      const SPoint& point4,
                                      const std::function<double(double,double)>& fun);

// Computes the integral over a quadrilateral shape
// Uses same number of quadrature nodes and weights in the x and y direction
template <int Rule>
inline double IntegrateOverQuadrilateral(const SPoint& point1,
                                         const SPoint& point2,
                                         const SPoint& point3,
                                         const SPoint& point4,
                                         const std::function<double(double,double)>& fun)
{
  return CalculateIntegralOverRectangle(Rule,
                                        point1,
                                        point2,
                                        point3,
                                        point4,
                                        fun);
}
 
template <int Rule>
inline double IntegrateOverQuadrilateral(const std::vector<SPoint>& points,
                                         const std::function<double(double,double)>& fun)
{
  assert(int(points.size())==4);
  return CalculateIntegralOverRectangle(Rule, 
                                        points[0],
                                        points[1],
                                        points[2],
                                        points[3],
                                        fun);
}

}
// end namespace FEQuadrature_Quadrilateral

}
// end namespace hydrofem

#endif /** __Hydrofem_FEQuadrature_Quadrilateral_HPP__ */
