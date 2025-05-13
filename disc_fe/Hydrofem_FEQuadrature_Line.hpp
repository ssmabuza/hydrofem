// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEQuadrature_Line_HPP__
#define __Hydrofem_FEQuadrature_Line_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEQuadrature.hpp"

namespace hydrofem
{

namespace FEQuadrature_Line
{
  
double CalculateIntegralOverLine(const int Rule, 
                                 const SPoint& point1, 
                                 const SPoint& point2, 
                                 const std::function<double(double)>& fun);  
  
template <int order>
inline double IntegrateOverLine(const SPoint& point1,
                                const SPoint& point2,
                                const std::function<double(double)>& fun)
{
  return CalculateIntegralOverLine(order,point1,point2,fun);
}

template <int order>
inline double IntegrateOverLine(const std::vector<SPoint>& points,
                                const std::function<double(double)>& fun)
{
  return CalculateIntegralOverLine(order,points.at(0),points.at(1),fun);
}

}
// end namespace FEQuadrature_Line

}
// end namespace hydrofem

#endif /** __Hydrofem_FEQuadrature_Line_HPP__ */
