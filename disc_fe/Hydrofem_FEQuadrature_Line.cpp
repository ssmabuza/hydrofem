// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEQuadrature_Line.hpp"

namespace hydrofem
{

using namespace FEQuadrature;
  
namespace FEQuadrature_Line
{
  
double CalculateIntegralOverLine(const int Rule, 
                                 const SPoint& point1, 
                                 const SPoint& point2, 
                                 const std::function<double(double)>& fun)
{
  const int rule_ = Rule;
  const double hx = std::fabs(0.5*(point2.x()-point1.x()));
  const double cx = 0.5*(point2.x()+point1.x());
  double res = 0.0;
  if (rule_ == 1)
    res = 2.0*hx*fun(cx);
  else if (rule_ == 2)
    for (int k = 0; k < 2; ++k)
      res += hx*weightsGauss1D3[k]*fun(cx + coordsGauss1D3[k]*hx);
  else if (rule_ == 3)
    for (int k = 0; k < 3; ++k)
      res += hx*weightsGauss1D5[k]*fun(cx + coordsGauss1D5[k]*hx);
  else if (rule_ == 4)
    for (int k = 0; k < 4; ++k)
      res += hx*weightsGauss1D7[k]*fun(cx + coordsGauss1D7[k]*hx);
  else if (rule_ == 5)
    for (int k = 0; k < 5; ++k)
      res += hx*weightsGauss1D9[k]*fun(cx + coordsGauss1D9[k]*hx);
  else if (rule_ == 6)
    for (int k = 0; k < 6; ++k)
      res += hx*weightsGauss1D11[k]*fun(cx + coordsGauss1D11[k]*hx);
  else if (rule_ == 7)
    for (int k = 0; k < 7; ++k)
      res += hx*weightsGauss1D13[k]*fun(cx + coordsGauss1D13[k]*hx);
  else if (rule_ == 8)
    for (int k = 0; k < 8; ++k)
      res += hx*weightsGauss1D15[k]*fun(cx + coordsGauss1D15[k]*hx);
  else if (rule_ == 9)
    for (int k = 0; k < 9; ++k)
      res += hx*weightsGauss1D17[k]*fun(cx + coordsGauss1D17[k]*hx);
  else if (rule_ == 10)
    for (int k = 0; k < 10; ++k)
      res += hx*weightsGauss1D19[k]*fun(cx + coordsGauss1D19[k]*hx);
  else if (rule_ == 11)
    for (int k = 0; k < 11; ++k)
      res += hx*weightsGauss1D21[k]*fun(cx + coordsGauss1D21[k]*hx);
  else {
    std::cerr << "Invaling Gauss quadrature rule with " << rule_ << " points, chosen.\n";
    std::cerr << "Only rules up to " << 11 << " points allowed." << std::endl;
    exit(1);
  }
  
  return res;
}
/// end double CalculateIntegralOverEdge

}
// end namespace FEQuadrature_Line

}
// end namespace hydrofem
