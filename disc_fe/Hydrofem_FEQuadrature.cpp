// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEQuadrature.hpp"
#include "Hydrofem_ElementShapeTools.hpp"

namespace hydrofem
{

namespace FEQuadrature
{

double CalculateIntegralOverEdge1(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun)
{
  const SPoint midpoint = 0.5*(point1 + point2);
  const SPoint AB = point2 - point1;
  const double length = AB.norm();
  return length*fun(midpoint.x(), midpoint.y());
}
/// end double CalculateIntegralOverEdge1

double CalculateIntegralOverEdge1(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers)
{
  const SPoint midpoint = 0.5*(point1 + point2);
  const SPoint AB = point2 - point1;
  const double length = AB.norm();
  double prod = 1.0;
  for (int i = 0; i < int(funs.size()); ++i)
    prod *= funs[i](midpoint.x(),midpoint.y());
  for (int i = 0; i < int(multipliers.size()); ++i)
    prod *= multipliers[i];
  return length*prod;
}
/// end double CalculateIntegralOverEdge1

double CalculateIntegralOverEdge3(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun)
{
  const double hx = 0.5*(point2.x()-point1.x());
  const double hy = 0.5*(point2.y()-point1.y());
  const double cx = 0.5*(point2.x()+point1.x());
  const double cy = 0.5*(point2.y()+point1.y());
  double x, y;
  double res = 0.0;
  for (int k = 0; k < 2; ++k)
  {
    x = cx + coordsGauss1D3[k]*hx;
    y = cy + coordsGauss1D3[k]*hy;
    res += weightsGauss1D3[k]*fun(x,y);
  }
  res *= std::sqrt(hx*hx+hy*hy);
  return res;
}
/// end double CalculateIntegralOverEdge3

double CalculateIntegralOverEdge3(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers)
{
  const double hx = 0.5*(point2.x()-point1.x());
  const double hy = 0.5*(point2.y()-point1.y());
  const double cx = 0.5*(point2.x()+point1.x());
  const double cy = 0.5*(point2.y()+point1.y());
  double x, y;
  double res = 0.0;
  for (int k = 0; k < 2; ++k)
  {
    x = cx + coordsGauss1D3[k]*hx;
    y = cy + coordsGauss1D3[k]*hy;
    double prod = 1.0;
    for (int i = 0; i < int(funs.size()); ++i)
      prod *= funs[i](x,y);
    for (int i = 0; i < int(multipliers.size()); ++i)
      prod *= multipliers[i];
    res += weightsGauss1D3[k]*prod;
  }
  res *= std::sqrt(hx*hx+hy*hy);
  return res;
}
/// end double CalculateIntegralOverEdge3

double CalculateIntegralOverEdge5(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun)
{
  const double hx = 0.5*(point2.x()-point1.x());
  const double hy = 0.5*(point2.y()-point1.y());
  const double cx = 0.5*(point2.x()+point1.x());
  const double cy = 0.5*(point2.y()+point1.y());
  double x, y;
  double res = 0.0;
  for (int k = 0; k < 3; ++k)
  {
    x = cx + coordsGauss1D5[k]*hx;
    y = cy + coordsGauss1D5[k]*hy;
    res += weightsGauss1D5[k]*fun(x,y);
  }
  res *= std::sqrt(hx*hx+hy*hy);
  return res;
}
/// end double CalculateIntegralOverEdge5

double CalculateIntegralOverEdge5(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers)
{
  const double hx = 0.5*(point2.x()-point1.x());
  const double hy = 0.5*(point2.y()-point1.y());
  const double cx = 0.5*(point2.x()+point1.x());
  const double cy = 0.5*(point2.y()+point1.y());
  double x, y;
  double res = 0.0;
  for (int k = 0; k < 3; ++k)
  {
    x = cx + coordsGauss1D5[k]*hx;
    y = cy + coordsGauss1D5[k]*hy;
    double prod = 1.0;
    for (int i = 0; i < int(funs.size()); ++i)
      prod *= funs[i](x,y);
    for (int i = 0; i < int(multipliers.size()); ++i)
      prod *= multipliers[i];
    res += weightsGauss1D5[k]*prod;
  }
  res *= sqrt(hx*hx+hy*hy);
  return res;
}
/// end double CalculateIntegralOverEdge5

double CalculateIntegralOverEdge(const int Rule, 
                                 const SPoint& point1, 
                                 const SPoint& point2, 
                                 const std::function<double(double,double)>& fun)
{
  const auto rule_ = Rule;
  if (rule_ == 1)
    return CalculateIntegralOverEdge1(point1,point2,fun);
  else if (rule_ == 2)
    return CalculateIntegralOverEdge3(point1,point2,fun);
  else if (rule_ == 3)
    return CalculateIntegralOverEdge5(point1,point2,fun);
  else if (rule_ == 4)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 4; ++k)
    {
      x = cx + coordsGauss1D7[k]*hx;
      y = cy + coordsGauss1D7[k]*hy;
      res += weightsGauss1D7[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 5)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 5; ++k)
    {
      x = cx + coordsGauss1D9[k]*hx;
      y = cy + coordsGauss1D9[k]*hy;
      res += weightsGauss1D9[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 6)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 6; ++k)
    {
      x = cx + coordsGauss1D11[k]*hx;
      y = cy + coordsGauss1D11[k]*hy;
      res += weightsGauss1D11[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 7)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 7; ++k)
    {
      x = cx + coordsGauss1D13[k]*hx;
      y = cy + coordsGauss1D13[k]*hy;
      res += weightsGauss1D13[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 8)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 8; ++k)
    {
      x = cx + coordsGauss1D15[k]*hx;
      y = cy + coordsGauss1D15[k]*hy;
      res += weightsGauss1D15[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 9)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 9; ++k)
    {
      x = cx + coordsGauss1D17[k]*hx;
      y = cy + coordsGauss1D17[k]*hy;
      res += weightsGauss1D17[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 10)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 10; ++k)
    {
      x = cx + coordsGauss1D19[k]*hx;
      y = cy + coordsGauss1D19[k]*hy;
      res += weightsGauss1D19[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 11)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 11; ++k)
    {
      x = cx + coordsGauss1D21[k]*hx;
      y = cy + coordsGauss1D21[k]*hy;
      res += weightsGauss1D21[k]*fun(x,y);
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  } else {
    std::string err_msg = "Invaling Gauss quadrature rule with " + std::to_string(rule_) +
                          " points, chosen.\n Only rules up to " + std::to_string(11) + " points allowed.";
    throw std::logic_error(err_msg);
  }
  
}
/// end double CalculateIntegralOverEdge

double CalculateIntegralOverEdge(const int Rule,
                                 const SPoint& point1,
                                 const SPoint& point2,
                                 const std::vector<std::function<double(double,double)>>& funs,
                                 const std::vector<double>& multipliers)
{
  
  auto rule_ = Rule;
  if (rule_ == 1)
    return CalculateIntegralOverEdge1(point1,point2,funs,multipliers);
  else if (rule_ == 2)
    return CalculateIntegralOverEdge3(point1,point2,funs,multipliers);
  else if (rule_ == 3)
    return CalculateIntegralOverEdge5(point1,point2,funs,multipliers);
  else if (rule_ == 4)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 4; ++k)
    {
      x = cx + coordsGauss1D7[k]*hx;
      y = cy + coordsGauss1D7[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D7[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 5)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 4; ++k)
    {
      x = cx + coordsGauss1D9[k]*hx;
      y = cy + coordsGauss1D9[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D9[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 6)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 6; ++k)
    {
      x = cx + coordsGauss1D11[k]*hx;
      y = cy + coordsGauss1D11[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D11[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 7)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 7; ++k)
    {
      x = cx + coordsGauss1D13[k]*hx;
      y = cy + coordsGauss1D13[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D13[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 8)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 8; ++k)
    {
      x = cx + coordsGauss1D15[k]*hx;
      y = cy + coordsGauss1D15[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D15[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 9)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 9; ++k)
    {
      x = cx + coordsGauss1D17[k]*hx;
      y = cy + coordsGauss1D17[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D17[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 10)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 10; ++k)
    {
      x = cx + coordsGauss1D19[k]*hx;
      y = cy + coordsGauss1D19[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D19[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  }
  else if (rule_ == 11)
  {
    const double hx = 0.5*(point2.x()-point1.x());
    const double hy = 0.5*(point2.y()-point1.y());
    const double cx = 0.5*(point2.x()+point1.x());
    const double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    double res = 0.0;
    for (int k = 0; k < 11; ++k)
    {
      x = cx + coordsGauss1D21[k]*hx;
      y = cy + coordsGauss1D21[k]*hy;
      double prod = 1.0;
      for (int i = 0; i < int(funs.size()); ++i)
        prod *= funs[i](x,y);
      for (int i = 0; i < int(multipliers.size()); ++i)
        prod *= multipliers[i];
      res += weightsGauss1D21[k]*prod;
    }
    res *= sqrt(hx*hx+hy*hy);
    return res;    
  } else {
    std::string err_msg = "Invaling Gauss quadrature rule with " + std::to_string(rule_) +
                          " points, chosen.\n Only rules up to " + std::to_string(11) + " points allowed.";
    throw std::logic_error(err_msg);
  }
}
/// end double CalculateIntegralOverEdge

double CalculateFirstMoment5(const SPoint& point1,
                             const SPoint& point2, 
                             const std::function<double(double,double)>& fun)
{
  double res = 0.0;
  for (int k = 0; k < 3; ++k)
  {
    const double a = 0.5*(1.0 - coordsGauss1D5[k]);
    const double b = 0.5*(1.0 + coordsGauss1D5[k]);
    const double x = a*point1.x() + b*point2.x();
    const double y = a*point1.y() + b*point2.y();

    res += weightsGauss1D5[k]*coordsGauss1D5[k]*fun(x,y);
  }
  res *= 3.0;
  return res;
}
/// end double CalculateFirstMoment5

double CalculateFirstMoment5(const SPoint& point1, 
                            const SPoint& point2, 
                            const std::vector<std::function<double(double,double)>>& funs,
                            const std::vector<double>& multipliers)
{
  double res = 0.0;
  for (int k = 0; k < 3; ++k)
  {
    const double a = 0.5*(1.0 - coordsGauss1D5[k]);
    const double b = 0.5*(1.0 + coordsGauss1D5[k]);
    const double x = a*point1.x() + b*point2.x();
    const double y = a*point1.y() + b*point2.y();

    double prod = 1.0;
    for (int i = 0; i < int(funs.size()); ++i)
      prod *= funs[i](x,y);
    for (int i = 0; i < int(multipliers.size()); ++i)
      prod *= multipliers[i];
    res += weightsGauss1D5[k]*prod;
  }
  res *= 3.0;
  return res;
}
/// end double CalculateFirstMoment5

double CalculateIntegralOverTriangle1(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun)
{
  const double area = computeTriangleArea(point1, point2, point3);
  const SPoint center = (point1 + point2 + point3)/3.0;
  return area*fun(center.x(), center.y());
}
/// end double CalculateIntegralOverTriangle1

double CalculateIntegralOverTriangle1(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers)
{
  const double area = computeTriangleArea(point1, point2, point3);
  const SPoint center = (point1 + point2 + point3)/3.0;
  
  double prod = 1.0;
  for (int i = 0; i < int(funs.size()); ++i)
    prod *= funs[i](center.x(), center.y());
  for (int i = 0; i < int(multipliers.size()); ++i)
    prod *= multipliers[i];
  
  return area*prod;
}
/// end double CalculateIntegralOverTriangle1

double CalculateIntegralOverTriangle2(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun)
{
  double res = 0.0;
  const double area = computeTriangleArea(point1, point2, point3);
  for (int k = 0; k < 3; ++k)
  {
    const SPoint point = coordsTri2[3*k + 0]*point1 + coordsTri2[3*k + 1]*point2 + coordsTri2[3*k + 2]*point3;
    res += weightsTri2[k]*fun(point.x(), point.y());
  }
  return res*area;
}
/// end double CalculateIntegralOverTriangle2

double CalculateIntegralOverTriangle2(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers)
{
  double res = 0.0;
  const double area = computeTriangleArea(point1, point2, point3);
  for (int k = 0; k < 3; ++k)
  {
    const SPoint point = coordsTri2[3*k + 0]*point1 + coordsTri2[3*k + 1]*point2 + coordsTri2[3*k + 2]*point3;

    double prod = 1.0;
    for (int i = 0; i < int(funs.size()); ++i)
      prod *= funs[i](point.x(), point.y());
    for (int i = 0; i < int(multipliers.size()); ++i)
      prod *= multipliers[i];
    
    res += weightsTri2[k]*prod;
  }
  return res*area;
}
/// end double CalculateIntegralOverTriangle2

double CalculateIntegralOverTriangle5(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun)
{
  double res = 0.0;
  const double area = computeTriangleArea(point1, point2, point3);
  for (int k = 0; k < 7; ++k)
  {
    const SPoint point = coordsTri5[3*k + 0]*point1 + coordsTri5[3*k + 1]*point2 + coordsTri5[3*k + 2]* point3;
    res += weightsTri5[k]*fun(point.x(), point.y());
  }
  return res*area;
}
/// end double CalculateIntegralOverTriangle5

double CalculateIntegralOverTriangle5(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers)
{
  double res = 0.0;
  const double area = computeTriangleArea(point1, point2, point3);
  for (int k = 0; k < 7; ++k)
  {
    const SPoint point = coordsTri5[3*k + 0]*point1 + coordsTri5[3*k + 1]*point2 + coordsTri5[3*k + 2]*point3;
    double prod = 1.0;
    for (int i = 0; i < int(funs.size()); ++i)
      prod *= funs[i](point.x(), point.y());
    for (int i = 0; i < int(multipliers.size()); ++i)
      prod *= multipliers[i];
    res += weightsTri5[k]*prod;
  }
  return res*area;
}
/// end double CalculateIntegralOverTriangle5

}
// end namespace FEQuadrature

}
// end namespace hydrofem
