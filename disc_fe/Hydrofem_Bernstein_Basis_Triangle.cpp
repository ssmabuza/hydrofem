// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_Bernstein_Basis_Triangle.hpp"

namespace hydrofem
{

/// evaluation on reference triangle
double Bernstein_Basis_Triangle::
eval(const double x,
     const double y,
     const int i,
     const int j,
     const int k) const
{
  return eval(x, y, 1.0 - x - y, i, j, k);
}

/// evaluation on arbitrary triangle using barycentric coordinates
double Bernstein_Basis_Triangle::
eval(const double u,
     const double v,
     const double w,
     const int i,
     const int j,
     const int k) const
{
  assert((0 <= i) && (0 <= j) && (0 <= k) && (int(p_) == i + j + k));
  assert(std::fabs(1.0 - u - v - w) < double(1.0e-8));
  return Intern::factorial(int(p_)) / (Intern::factorial(i) * Intern::factorial(j) * Intern::factorial(k)) 
          * std::pow(u, double(i)) * std::pow(v, double(j)) * std::pow(w, double(k));
}

/// evaluation on arbitrary triangle using barycentric coordinates
double Bernstein_Basis_Triangle::
eval(const double x,
     const double y,
     const SPoint& tri1,
     const SPoint& tri2,
     const SPoint& tri3,
     const int i,
     const int j,
     const int k) const
{
  double u, v, w;
  Intern::BarycentricCoordinates(SPoint(x,y), tri1, tri2, tri3, u, v, w);
  return eval(u,v,w,i,j,k);
}

double Bernstein_Basis_Triangle::
eval_u(const double u,
       const double v,
       const double w,
       const int i,
       const int j,
       const int k) const
{
  assert(i >= 0);
  if (i == 0)
    return 0.0;
  else
    return Intern::factorial(int(p_)) / (Intern::factorial(i-1) * Intern::factorial(j) * Intern::factorial(k)) 
            * std::pow(u, double(i-1)) * std::pow(v, double(j)) * std::pow(w, double(k));
}

double Bernstein_Basis_Triangle::
eval_v(const double u,
       const double v,
       const double w,
       const int i,
       const int j,
       const int k) const
{
  assert(j >= 0);
  if (j == 0)
    return 0.0;
  else
    return Intern::factorial(int(p_)) / (Intern::factorial(i) * Intern::factorial(j-1) * Intern::factorial(k)) 
            * std::pow(u, double(i)) * std::pow(v, double(j-1)) * std::pow(w, double(k));
}

double Bernstein_Basis_Triangle::
eval_w(const double u,
       const double v,
       const double w,
       const int i,
       const int j,
       const int k) const
{
  assert(k >= 0);
  if (k == 0)
    return 0.0;
  else
    return Intern::factorial(int(p_)) / (Intern::factorial(i) * Intern::factorial(j) * Intern::factorial(k-1)) 
            * std::pow(u, double(i)) * std::pow(v, double(j)) * std::pow(w, double(k-1));
}

double Bernstein_Basis_Triangle::
eval_x(const double x,
       const double y, 
       const SPoint& tri1,
       const SPoint& tri2,
       const SPoint& tri3, 
       const int i,
       const int j,
       const int k) const
{
  double u, v, w;
  double u_x, v_x, w_x;
  
  Intern::BarycentricCoordinates(SPoint(x,y),tri1,tri2,tri3,u,v,w);
  Intern::BarycentricCoordinates_x(SPoint(x,y),tri1,tri2,tri3,u_x,v_x,w_x);
  
  return u_x*eval_u(u,v,w,i,j,k) + 
         v_x*eval_v(u,v,w,i,j,k) + 
         w_x*eval_w(u,v,w,i,j,k);
}

double Bernstein_Basis_Triangle::
eval_y(const double x,
       const double y, 
       const SPoint& tri1,
       const SPoint& tri2,
       const SPoint& tri3, 
       const int i,
       const int j,
       const int k) const
{
  double u, v, w;
  double u_y, v_y, w_y;
  
  Intern::BarycentricCoordinates(SPoint(x,y),tri1,tri2,tri3,u,v,w);
  Intern::BarycentricCoordinates_y(SPoint(x,y),tri1,tri2,tri3,u_y,v_y,w_y);
  
  return u_y*eval_u(u,v,w,i,j,k) + 
         v_y*eval_v(u,v,w,i,j,k) + 
         w_y*eval_w(u,v,w,i,j,k);
}

SPoint Bernstein_Basis_Triangle::
grad(const double x,
     const double y, 
     const SPoint& tri1,
     const SPoint& tri2,
     const SPoint& tri3, 
     const int i,
     const int j,
     const int k) const
{
  double u, v, w;
  double u_x, v_x, w_x;
  double u_y, v_y, w_y;
  
  Intern::BarycentricCoordinates(SPoint(x,y),tri1,tri2,tri3,u,v,w);
  Intern::BarycentricCoordinates_x(SPoint(x,y),tri1,tri2,tri3,u_x,v_x,w_x);
  Intern::BarycentricCoordinates_y(SPoint(x,y),tri1,tri2,tri3,u_y,v_y,w_y);
  
  double B_u (eval_u(u,v,w,i,j,k));
  double B_v (eval_v(u,v,w,i,j,k));
  double B_w (eval_w(u,v,w,i,j,k));
  
  return SPoint(u_x*B_u + v_x*B_v + w_x*B_w,u_y*B_u + v_y*B_v + w_y*B_w);
}

}
// end namespace hydrofem
