// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEQuadrature_Dunavant_HPP__
#define __Hydrofem_FEQuadrature_Dunavant_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

namespace FEQuadrature_Dunavant
{

double CalculateIntegralOverTriangle(const int Order,
                                     const SPoint& point1, 
                                     const SPoint& point2, 
                                     const SPoint& point3, 
                                     const std::function<double(double,double)>& fun);

double CalculateIntegralOverTriangle(const int Order,
                                     const SPoint& point1,
                                     const SPoint& point2,
                                     const SPoint& point3,
                                     const std::vector<std::function<double(double,double)>>& funs,
                                     const std::vector<double>& multipliers);

template<int Order>
inline double IntegrateOverTriangle(const SPoint& point1,
                                    const SPoint& point2,
                                    const SPoint& point3,
                                    const std::function<double(double,double)>& fun)
{
  return CalculateIntegralOverTriangle(Order,
                                       point1, 
                                       point2, 
                                       point3, 
                                       fun);
}
  
template<int Order>
inline double IntegrateOverTriangle(const SPoint& point1, 
                                    const SPoint& point2, 
                                    const SPoint& point3, 
                                    const std::vector<std::function<double(double,double)>>& funs,
                                    const std::vector<double>& multipliers)
{
  return CalculateIntegralOverTriangle(Order,
                                       point1, 
                                       point2, 
                                       point3, 
                                       funs,
                                       multipliers);
}

template<int Order>
inline double IntegrateOverTriangle(const std::vector<SPoint>& points,
                                    const std::function<double(double,double)>& fun)
{
  assert(int(points.size())==3);
  return CalculateIntegralOverTriangle(Order,
                                       points[0], 
                                       points[1], 
                                       points[2], 
                                       fun);
}
  
template<int Order>
inline double IntegrateOverTriangle(const std::vector<SPoint>& points,
                                    const std::vector<std::function<double(double,double)>>& funs,
                                    const std::vector<double>& multipliers)
{
  assert(int(points.size())==3);
  return CalculateIntegralOverTriangle(Order,
                                       points[0], 
                                       points[1], 
                                       points[2], 
                                       funs,
                                       multipliers);
}

}
// end namespace FEQuadrature_Dunavant

}
// end namespace hydrofem

#endif /** __Hydrofem_FEQuadrature_Dunavant_HPP__ */
