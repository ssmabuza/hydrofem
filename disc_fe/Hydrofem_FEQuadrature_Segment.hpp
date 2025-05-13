// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEQuadrature_Segment_HPP__
#define __Hydrofem_FEQuadrature_Segment_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEQuadrature.hpp"

namespace hydrofem
{

namespace FEQuadrature_Segment
{
  
template<int Order>
inline 
double IntegrateOverEdge(const SPoint& point1,
                         const SPoint& point2,
                         const std::function<double(double,double)>& fun)
{
  return FEQuadrature::CalculateIntegralOverEdge(Order, point1, point2, fun);
}

template<int Order>
inline 
double IntegrateOverEdge(const std::vector<SPoint>& points,
                         const std::function<double(double,double)>& fun)
{
  assert(points.size() == static_cast<std::size_t>(2));
  return FEQuadrature::CalculateIntegralOverEdge(Order, points[0], points[1], fun);
}

template<int Order>
inline 
double IntegrateOverEdge(const SPoint& point1,
                         const SPoint& point2,
                         const std::vector<std::function<double(double,double)>>& funs,
                         const std::vector<double>& multipliers)
{
  return FEQuadrature::CalculateIntegralOverEdge(Order, point1, point2, funs, multipliers);
}

template<int Order>
inline 
double IntegrateOverEdge(const std::vector<SPoint>& points, 
                         const std::vector<std::function<double(double,double)>>& funs,
                         const std::vector<double>& multipliers)
{
  assert(points.size() == static_cast<std::size_t>(2));
  return FEQuadrature::CalculateIntegralOverEdge(Order, points[0], points[1], funs, multipliers);
}

}
// end namespace FEQuadrature_Segment

}
// end namespace hydrofem

#endif /** __Hydrofem_FEQuadrature_Segment_HPP__ */
