// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Bernstein_Basis_Line.hpp"

namespace hydrofem
{

double Bernstein_Basis_Line::
eval(const double x, const int k) const
{
  assert((k <= m_p) && (k >= 0));
  return (Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k))) * std::pow(x,k) * std::pow(1.0-x,m_p-k);  
}
  
double Bernstein_Basis_Line::
eval(const double x,
     const SPoint& line1,
     const SPoint& line2,
     const int k) const
{
  return eval((x - line1.x())/(line2.x()-line1.x()),k);
}

double Bernstein_Basis_Line::
eval_x(const double x, const int k) const
{
  double value;
  if (k == 0)
    value = -m_p * (Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k))) * std::pow(1.0-x,m_p-1);
  else if (k == m_p)
    value = m_p*(Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k)))*std::pow(x,m_p-1);
  else 
    value = k*(Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k)))*std::pow(x,k-1)*std::pow(1.0-x,m_p-k)
            - (m_p-k)* (Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k))) *std::pow(x,k)*std::pow(1.0-x,m_p-k-1);
  return value;
}

double Bernstein_Basis_Line::
eval_x(const double x,
       const SPoint& line1,
       const SPoint& line2,
       const int k) const
{
  return eval_x((x - line1.x())/(line2.x()-line1.x()),k) * x/(line2.x()-line1.x());
}

}
// end namespace valiant
