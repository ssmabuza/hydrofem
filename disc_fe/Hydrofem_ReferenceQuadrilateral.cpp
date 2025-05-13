// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_ElementShapeTools.hpp"
#include "Hydrofem_ReferenceQuadrilateral.hpp"

namespace hydrofem
{

double quadArea(const SPoint& point1,
                const SPoint& point2,
                const SPoint& point3,
                const SPoint& point4)
{
  return computeTriangleArea(point1,point2,point3)+computeTriangleArea(point1,point3,point4);
}
  
double getJacobianDeterminant(const SPoint& ref_point,
                              const SPoint& point1,
                              const SPoint& point2,
                              const SPoint& point3,
                              const SPoint& point4)
{
  double Xi  = 2.0 * ref_point.x() - 1.0;
  double Eta = 2.0 * ref_point.y() - 1.0;
  
  double dN1_dXi  = -0.25*(1.0 - Eta);
  double dN1_dEta = -0.25*(1.0 - Xi);
  double dN2_dXi  =  0.25*(1.0 - Eta);
  double dN2_dEta = -0.25*(1.0 + Xi);
  double dN3_dXi  =  0.25*(1.0 + Eta);
  double dN3_dEta =  0.25*(1.0 + Xi);
  double dN4_dXi  = -0.25*(1.0 + Eta);
  double dN4_dEta =  0.25*(1.0 - Xi);
  
  double dx_dXi  = point1.x()*dN1_dXi  + point2.x()*dN2_dXi  + point3.x()*dN3_dXi  + point4.x()*dN4_dXi;
  double dy_dXi  = point1.y()*dN1_dXi  + point2.y()*dN2_dXi  + point3.y()*dN3_dXi  + point4.y()*dN4_dXi;
  double dx_dEta = point1.x()*dN1_dEta + point2.x()*dN2_dEta + point3.x()*dN3_dEta + point4.x()*dN4_dEta;
  double dy_dEta = point1.y()*dN1_dEta + point2.y()*dN2_dEta + point3.y()*dN3_dEta + point4.y()*dN4_dEta;
  
  return dx_dXi*dy_dEta - dx_dEta*dy_dXi;
}
  
SPoint physicalQuadToReference(const SPoint& quad1,
                               const SPoint& quad2,
                               const SPoint& quad3,
                               const SPoint& quad4,
                               const SPoint& physicalPoint)
{
  double x1 = quad1.x();
  double x2 = quad2.x();
  double x3 = quad3.x();
  double x4 = quad4.x();
  double y1 = quad1.y();
  double y2 = quad2.y();
  double y3 = quad3.y();
  double y4 = quad4.y();
  
  double x = physicalPoint.x();
  double y = physicalPoint.y();
  
  //////////////////////////////////////////////////////////////////
  // Implementation from Section 23-13 IFEM notes //////////////////
  //////////////////////////////////////////////////////////////////
  
  double xb    = x1 - x2 + x3 - x4;
  double yb    = y1 - y2 + y3 - y4;
  double x_cx  = x1 + x2 - x3 - x4;
  double y_cx  = y1 + y2 - y3 - y4;
  double x_ce  = x1 - x2 - x3 + x4;
  double y_ce  = y1 - y2 - y3 + y4;
  double A     = 0.5*((x3 - x1)*(y4 - y2) - (x4 - x2)*(y3 - y1));
  double J1    = (x3 - x4)*(y1 - y2) - (x1 - x2)*(y3 - y4);
  double J2    = (x2 - x3)*(y1 - y4) - (x1 - x4)*(y2 - y3);
  double x0    = 0.25*(x1 + x2 + x3 + x4);
  double y0    = 0.25*(y1 + y2 + y3 + y4);
  double xP0   = x - x0;
  double yP0   = y - y0;
  double b_xi  = A - xP0*yb + yP0*xb;
  double b_eta = -A - xP0*yb + yP0*xb;
  double c_xi  = xP0*y_cx - yP0*x_cx;
  double c_eta = xP0*y_ce - yP0*x_ce;
  
  // transform to [-1,1] x [-1,1]
  double xi,eta;
  xi  = 2*c_xi/(-std::sqrt(b_xi*b_xi + 2*J1*c_xi)    - b_xi);
  eta = 2*c_eta/(std::sqrt(b_eta*b_eta + 2*J2*c_eta) - b_eta);
  
  // transform to [0,1] x [0,1] for Bezier unit cell
  double xi_hat = 0.5*(xi+1.0);
  double eta_hat = 0.5*(eta+1.0);
  
  return SPoint(xi_hat,eta_hat);
}
  
SPoint referenceQuadToPhysical(const SPoint& quad1,
                               const SPoint& quad2,
                               const SPoint& quad3,
                               const SPoint& quad4,
                               const SPoint& referencePoint)
{
  // first transform this from [0,1] x [0,1] to [-1,1] x [-1,1]
  double Xi = 2.0 * referencePoint.x() - 1.0;
  double Eta = 2.0 * referencePoint.y() - 1.0;
  
  // now start transforming to general quad
  
  double N1 = 0.25*(1.0-Xi)*(1.0-Eta);
  double N2 = 0.25*(1.0+Xi)*(1.0-Eta);
  double N3 = 0.25*(1.0+Xi)*(1.0+Eta);
  double N4 = 0.25*(1.0-Xi)*(1.0+Eta);
  
  double x1 = quad1.x();
  double x2 = quad2.x();
  double x3 = quad3.x();
  double x4 = quad4.x();
  double y1 = quad1.y();
  double y2 = quad2.y();
  double y3 = quad3.y();
  double y4 = quad4.y();
  
  // get the physical points
  return SPoint(x1*N1+x2*N2+x3*N3+x4*N4,
                y1*N1+y2*N2+y3*N3+y4*N4);
}
  
LMAT_<double> getJacobian(const SPoint& ref_point,
                          const SPoint& point1,
                          const SPoint& point2,
                          const SPoint& point3,
                          const SPoint& point4)
{
  
  double Xi  = 2.0 * ref_point.x() - 1.0;
  double Eta = 2.0 * ref_point.y() - 1.0;
  
  double dN1_dXi  = -0.25*(1.0 - Eta);
  double dN1_dEta = -0.25*(1.0 - Xi);
  double dN2_dXi  =  0.25*(1.0 - Eta);
  double dN2_dEta = -0.25*(1.0 + Xi);
  double dN3_dXi  =  0.25*(1.0 + Eta);
  double dN3_dEta =  0.25*(1.0 + Xi);
  double dN4_dXi  = -0.25*(1.0 + Eta);
  double dN4_dEta =  0.25*(1.0 - Xi);
  
  double dx_dXi  = point1.x()*dN1_dXi  + point2.x()*dN2_dXi  + point3.x()*dN3_dXi  + point4.x()*dN4_dXi;
  double dy_dXi  = point1.y()*dN1_dXi  + point2.y()*dN2_dXi  + point3.y()*dN3_dXi  + point4.y()*dN4_dXi;
  double dx_dEta = point1.x()*dN1_dEta + point2.x()*dN2_dEta + point3.x()*dN3_dEta + point4.x()*dN4_dEta;
  double dy_dEta = point1.y()*dN1_dEta + point2.y()*dN2_dEta + point3.y()*dN3_dEta + point4.y()*dN4_dEta;

  LMAT_<double> jac = createKArray<LMAT_<double>>(2,2);
  jac(0,0) = dx_dXi;
  jac(0,1) = dy_dXi;
  jac(1,0) = dx_dEta;
  jac(1,1) = dy_dEta;
  
  return jac;  
}

LMAT_<double> getJacobianInverse(const SPoint& quad1,
                                 const SPoint& quad2,
                                 const SPoint& quad3,
                                 const SPoint& quad4,
                                 const double x, 
                                 const double y)
{
  auto ref_point = physicalQuadToReference(quad1,quad2,quad3,quad4,SPoint(x,y));
  auto jac = getJacobian(ref_point,quad1,quad2,quad3,quad4);
  auto jac_inv = createKArray<LMAT_<double>>(2,2);
  auto jac_det = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
  jac_inv(0,0) =  jac(1,1)/jac_det;
  jac_inv(0,1) = -jac(0,1)/jac_det;
  jac_inv(1,0) = -jac(1,0)/jac_det;
  jac_inv(1,1) =  jac(0,0)/jac_det;
  return jac_inv;
}

double deta_dx(const SPoint& quad1,
               const SPoint& quad2,
               const SPoint& quad3,
               const SPoint& quad4,
               const double x, 
               const double y)
{
  auto ref_point = physicalQuadToReference(quad1,quad2,quad3,quad4,SPoint(x,y));
  auto jac = getJacobian(ref_point,quad1,quad2,quad3,quad4);
  auto jac_det = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
  return -jac(0,1)/jac_det;
}

double deta_dy(const SPoint& quad1,
               const SPoint& quad2,
               const SPoint& quad3,
               const SPoint& quad4,
               const double x, 
               const double y)
{
  auto ref_point = physicalQuadToReference(quad1,quad2,quad3,quad4,SPoint(x,y));
  auto jac = getJacobian(ref_point,quad1,quad2,quad3,quad4);
  auto jac_det = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
  return jac(0,0)/jac_det;
}

double dxi_dx(const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const double x, 
              const double y)
{
  auto ref_point = physicalQuadToReference(quad1,quad2,quad3,quad4,SPoint(x,y));
  auto jac = getJacobian(ref_point,quad1,quad2,quad3,quad4);
  auto jac_det = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
  return jac(1,1)/jac_det;
}

double dxi_dy(const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const double x, 
              const double y)
{
  auto ref_point = physicalQuadToReference(quad1,quad2,quad3,quad4,SPoint(x,y));
  auto jac = getJacobian(ref_point,quad1,quad2,quad3,quad4);
  auto jac_det = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
  return -jac(1,0)/jac_det;
}

}
// end namespace hydrofem

