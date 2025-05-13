// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEUtils_HPP__
#define __Hydrofem_FEUtils_HPP__

#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

/// to_string with leading 0
std::string to_string_0(const int& i, const int length = int(5));

// definition of pi
const double PI = double(4.0) * std::atan(1.0);

namespace Intern
{
/// calculate factorial
int factorial(const int n);

class Binomial
{
private:
  
  int _len_binomialMat;
  int ** _binomialMat;

public:
  
  Binomial(const int len_binomialMat);

  ~Binomial();

  int operator()(const int i, const int j) const;
};

// Binomial binomial(40);

/// calculate barycentric coordinates
void BarycentricCoordinates(const SPoint& p,
                            const SPoint& a,
                            const SPoint& b,
                            const SPoint& c,
                            double& u, 
                            double& v, 
                            double& w);

void BarycentricCoordinates_x(const SPoint& p,
                              const SPoint& a,
                              const SPoint& b,
                              const SPoint& c,
                              double& u_x,
                              double& v_x,
                              double& w_x);

void BarycentricCoordinates_y(const SPoint& p,
                              const SPoint& a,
                              const SPoint& b,
                              const SPoint& c,
                              double& u_y,
                              double& v_y,
                              double& w_y);


const Binomial binomial(40);

}
// end namespace Intern

}
// end namespace hydrofem

#endif /** __Hydrofem_FEUtils_HPP__ */
