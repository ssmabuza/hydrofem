// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEUtils.hpp"

namespace hydrofem
{

/// to_string with leading 0
std::string to_string_0(const int& i, const int length)
{
  const auto old_string = std::to_string(i);
  return std::string(length - old_string.length(), '0') + old_string;
}

namespace Intern
{
/// calculate factorial
int factorial(const int n)
{
  int ret = 1;
  for(int i(1); i <= n; ++i)
  {
    ret *= i;
  }

  return ret;
}

  
Binomial::Binomial(const int len_binomialMat) : _len_binomialMat(len_binomialMat)
{
  _binomialMat = new int *[_len_binomialMat];
  for (int i(0); i < _len_binomialMat; i++)
    _binomialMat[i]= new int [_len_binomialMat];

  for (int i(0); i < _len_binomialMat; i++)
  {
    for (int j(0); j < _len_binomialMat; j++)
    {
      _binomialMat[i][j]=0;
    }
  }

  for (int i(0); i < _len_binomialMat; i++)
    _binomialMat[i][0] += 1;

  for (int j(1); j < _len_binomialMat; j++)
    _binomialMat[0][j] += 1;

  for (int i(1); i < _len_binomialMat; i++)
  {
    for (int j(1); j < _len_binomialMat; j++)
    {
      _binomialMat[i][j] += _binomialMat[i][j-1] + _binomialMat[i-1][j];
    }
  }
}

Binomial::~Binomial()
{
  for (int i(0); i < _len_binomialMat; i++)
    delete[] _binomialMat[i];
  delete[] _binomialMat;
}

int Binomial::operator()(const int i, const int j) const
{
  assert((int(0) < i) && (int(0) < j) && (int(i + j) < _len_binomialMat));

  return int(_binomialMat[int(i) - int(1)][int(j) - int(1)]);
}

/// calculate barycentric coordinates
void BarycentricCoordinates(const SPoint& p,
                            const SPoint& a,
                            const SPoint& b,
                            const SPoint& c,
                            double& u, 
                            double& v, 
                            double& w)
{
  SPoint v0 = b - a;
  SPoint v1 = c - a;
  SPoint v2 = p - a;

  double d00 = v0 * v0;
  double d01 = v0 * v1;
  double d11 = v1 * v1;
  double d20 = v2 * v0;
  double d21 = v2 * v1;
  double denom = d00 * d11 - d01 * d01;
  v = (d11 * d20 - d01 * d21) / denom;
  w = (d00 * d21 - d01 * d20) / denom;
  u = double(1) - v - w;
}

void BarycentricCoordinates_x(const SPoint& /*p*/,
                              const SPoint& a,
                              const SPoint& b,
                              const SPoint& c,
                              double& u_x,
                              double& v_x,
                              double& w_x)
{
  double x1 = a.x();
  double x2 = b.x();
  double x3 = c.x();
  double y1 = a.y();
  double y2 = b.y();
  double y3 = c.y();
  
  double detT = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
  
  u_x = (y2-y3)/detT;
  v_x = (y3-y1)/detT;
  w_x = -u_x-v_x;
}

void BarycentricCoordinates_y(const SPoint& /*p*/,
                              const SPoint& a,
                              const SPoint& b,
                              const SPoint& c,
                              double& u_y,
                              double& v_y,
                              double& w_y)
{
  double x1 = a.x();
  double x2 = b.x();
  double x3 = c.x();
  double y1 = a.y();
  double y2 = b.y();
  double y3 = c.y();
  
  double detT = (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3);
  
  u_y = (x3-x2)/detT;
  v_y = (x1-x3)/detT;
  w_y = -u_y-v_y;
}

}
// end namespace Intern

}
// end namespace hydrofem
