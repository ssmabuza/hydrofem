// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_Bernstein_Basis_Line.hpp"
#include "Hydrofem_Bernstein_Basis_Quadrilateral.hpp"

namespace hydrofem
{

double Bernstein_Basis_Quadrilateral::eval(const double x,
                                           const double y,
                                           const int i,
                                           const int j) const
{
  return line_basis1->eval(x,i)*line_basis2->eval(y,j);
}
  
double
Bernstein_Basis_Quadrilateral::
eval(const double x,
     const double y,
     const SPoint& a,
     const SPoint& b,
     const SPoint& c,
     const SPoint& d,
     const int i,
     const int j) const
{
  const SPoint ref_pt = physicalQuadToReference(a,b,c,d,SPoint(x,y));
  double u, v;
  u = ref_pt.x();
  v = ref_pt.y();  
  return eval(u,v,i,j);
}

double
Bernstein_Basis_Quadrilateral::
eval_x(const double x,
       const double y,
       const int i,
       const int j) const
{
  return line_basis1->eval_x(x,i)*line_basis2->eval(y,j);
}

double
Bernstein_Basis_Quadrilateral::
eval_x(const double x,
       const double y,
       const SPoint& quad1,
       const SPoint& quad2,
       const SPoint& quad3,
       const SPoint& quad4,
       const int i,
       const int j) const
{
  double eps = 1.0e-9;
  double x_minus_eps = x - eps;
  double x_plus_eps = x + eps;
  return (eval(x_plus_eps,y,quad1,quad2,quad3,quad4,i,j)-
          eval(x_minus_eps,y,quad1,quad2,quad3,quad4,i,j))/(2*eps);
}

double
Bernstein_Basis_Quadrilateral::
eval_y(const double x,
       const double y,
       const int i,
       const int j) const
{
  return line_basis1->eval(x,i)*line_basis2->eval_x(y,j);
}

double
Bernstein_Basis_Quadrilateral::
eval_y(const double x,
       const double y,
       const SPoint& quad1,
       const SPoint& quad2,
       const SPoint& quad3,
       const SPoint& quad4,
       const int i,
       const int j) const
{
  double eps = 1.0e-9;
  double y_minus_eps = y - eps;
  double y_plus_eps = y + eps;
  return (eval(x,y_plus_eps,quad1,quad2,quad3,quad4,i,j) -
          eval(x,y_minus_eps,quad1,quad2,quad3,quad4,i,j))/(2*eps);
}

SPoint
Bernstein_Basis_Quadrilateral::
grad_eval(const double x,
     const double y,
     const SPoint& quad1,
     const SPoint& quad2,
     const SPoint& quad3,
     const SPoint& quad4,
     const int i,
     const int j) const
{
  return SPoint(eval_x(x,y,quad1,quad2,quad3,quad4,i,j),
                eval_y(x,y,quad1,quad2,quad3,quad4,i,j));
}

}
// end namespace hydrofem

