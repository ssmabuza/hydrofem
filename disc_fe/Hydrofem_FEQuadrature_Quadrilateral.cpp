// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_ElementShapeTools.hpp"
#include "Hydrofem_FEQuadrature_Quadrilateral.hpp"

namespace hydrofem
{

namespace FEQuadrature_Quadrilateral
{

double CalculateIntegralOverRectangle(const int Rule,
                                      const SPoint& point1,
                                      const SPoint& point2,
                                      const SPoint& point3,
                                      const SPoint& point4,
                                      const std::function<double(double,double)>& fun)
{
//   const int num_qp(5*5);
  // reference quadrature points (x_i,y_j)
  std::vector<SPoint> ref_qp;
  // quadrature weights w_iw_j
  std::vector<double> qweights;
  // get weights and nodes
  GetQuadrature(ref_qp,qweights,Rule);
  // map the reference quadrature points to current element
  std::vector<SPoint> qp;
  MapPhysicalRectangleToReferenceRectangle(qp,ref_qp,point1,point2,point3,point4);
  // get the element area
  // integrate:
  double res = 0.0;
  for (int qpInd = 0; qpInd < int(qp.size()); ++qpInd)
  {
    auto detJ = JacobianDeterminant(ref_qp[qpInd],point1,point2,point3,point4);
    res += qweights[qpInd]*fun(qp[qpInd].x(),qp[qpInd].y())*detJ;
  }
  return res;
}

void MapPhysicalRectangleToReferenceRectangle(std::vector<SPoint>& phy_points,
                                              const std::vector<SPoint>& ref_points,
                                              const SPoint& phy_point1,
                                              const SPoint& phy_point2,
                                              const SPoint& phy_point3,
                                              const SPoint& phy_point4)
{
  phy_points.clear();
  phy_points.resize(ref_points.size(),SPoint(2));
  for (auto i = 0; i < int(ref_points.size()); ++i)
  {
    // shape functions
    //@{
    double Xi = ref_points[i].x();
    double Eta = ref_points[i].y();
    double N1 = 0.25*(1.0-Xi)*(1.0-Eta);
    double N2 = 0.25*(1.0+Xi)*(1.0-Eta);
    double N3 = 0.25*(1.0+Xi)*(1.0+Eta);
    double N4 = 0.25*(1.0-Xi)*(1.0+Eta);
    //@}
    
    // get the physical points
    phy_points[i].x() = phy_point1.x()*N1+phy_point2.x()*N2+phy_point3.x()*N3+phy_point4.x()*N4;
    phy_points[i].y() = phy_point1.y()*N1+phy_point2.y()*N2+phy_point3.y()*N3+phy_point4.y()*N4;
  }
}

double ComputeElementArea(const SPoint& point1,
                          const SPoint& point2,
                          const SPoint& point3,
                          const SPoint& point4)
{
  return computeTriangleArea(point1,point2,point3)+computeTriangleArea(point1,point3,point4);
}
                          
double JacobianDeterminant(const SPoint& ref_point,
                           const SPoint& point1,
                           const SPoint& point2,
                           const SPoint& point3,
                           const SPoint& point4)
{
  double Xi = ref_point.x();
  double Eta = ref_point.y();
  double dN1_dXi  = -0.25*(1.0-Eta);
  double dN1_dEta = -0.25*(1.0-Xi);
  double dN2_dXi  =  0.25*(1.0-Eta);
  double dN2_dEta = -0.25*(1.0+Xi);
  double dN3_dXi  =  0.25*(1.0+Eta);
  double dN3_dEta =  0.25*(1.0+Xi);
  double dN4_dXi  = -0.25*(1.0+Eta);
  double dN4_dEta =  0.25*(1.0-Xi);
  
  double dx_dXi  = point1.x()*dN1_dXi+point2.x()*dN2_dXi+point3.x()*dN3_dXi+point4.x()*dN4_dXi;
  double dy_dXi  = point1.y()*dN1_dXi+point2.y()*dN2_dXi+point3.y()*dN3_dXi+point4.y()*dN4_dXi;
  double dx_dEta = point1.x()*dN1_dEta+point2.x()*dN2_dEta+point3.x()*dN3_dEta+point4.x()*dN4_dEta;
  double dy_dEta = point1.y()*dN1_dEta+point2.y()*dN2_dEta+point3.y()*dN3_dEta+point4.y()*dN4_dEta;
  
  return dx_dXi*dy_dEta - dx_dEta*dy_dXi;
}

void GetQuadrature(std::vector<SPoint>& qp,                                      
                   std::vector<double>& qweights,
                   const int rule_1d)
{
  qp.clear();
  qp.resize(rule_1d*rule_1d,SPoint(2));
  qweights.clear();
  qweights.resize(rule_1d*rule_1d);
  
  // 1D Gaussian Quadrature Weights and Abscissae from: https://pomax.github.io/bezierinfo/legendre-gauss.html
  if (rule_1d == 1)
  {
    /// 1 point Gauss in 1D, exact on linears : 1 node in 2D
    qp[0] = SPoint(0.0,0.0);
    qweights[0] = 1.0000000000000;
    
  } else if (rule_1d == 2) {
    
    /// 2 point Gauss in 1D, exact on degree 3 polys : 4 nodes in 2D
    const double coordsGauss[2] = {-0.57735026918962576451, 0.57735026918962576451};
    const double weightsGauss[2] = {1.0000000000000000, 1.0000000000000000};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 3) {
  
    /// 3 point Gauss in 1D, exact on degree 5 polys : 9 nodes in 2D
    const double coordsGauss[3] = {-0.7745966692414837704, 0.0000000000000000000, 0.7745966692414837704};
    const double weightsGauss[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 4) {
    
    /// 4 point Gauss in 1D, exact on degree 7 polys : 16 nodes in 2D
    const double coordsGauss[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
    const double weightsGauss[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 5) {
    
    /// 5 point Gauss in 1D, exact on degree 9 polys : 25 nodes in 2D
    const double coordsGauss[5] = {-0.9061798459386640, -0.5384693101056831, 0.0000000000000000, 0.5384693101056831, 0.9061798459386640};
    const double weightsGauss[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 6) {
    
    /// 6 point Gauss in 1D, exact on degree 11 polys : 36 nodes in 2D
    const double coordsGauss[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
    const double weightsGauss[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 7) {
    
    /// 7 point Gauss in 1D, exact on degree 13 polys : 49 nodes in 2D
    const double coordsGauss[7] = {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0000000000000000, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585};
    const double weightsGauss[7] = {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 8) {
    
    /// 8 point Gauss in 1D, exact on degree 15 polys : 64 nodes in 2D
    const double coordsGauss[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
    const double weightsGauss[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 9) {
    
    /// 9 point Gauss in 1D, exact on degree 17 polys : 81 nodes in 2D
    const double coordsGauss[9] = {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0000000000000000, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
    const double weightsGauss[9] = {0.3123470770400029, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, 0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.3123470770400029};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 10) {
    
    /// 10 point Gauss in 1D, exact on degree 19 polys : 100 nodes in 2D
    const double coordsGauss[10] = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
    const double weightsGauss[10] = {0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else if (rule_1d == 11) {
    
    /// 11 point Gauss in 1D, exact on degree 21 polys : 121 nodes in 2D
    const double coordsGauss[11] = {-0.9782286581460570, -0.8870625997680953, -0.7301520055740494, -0.5190961292068118, -0.2695431559523450, 0.0000000000000000, 0.2695431559523450, 0.5190961292068118, 0.7301520055740494, 0.8870625997680953, 0.9782286581460570};
    const double weightsGauss[11] = {0.0556685671161737, 0.1255803694649046, 0.1862902109277343, 0.2331937645919905, 0.2628045445102467, 0.2729250867779006, 0.2628045445102467, 0.2331937645919905, 0.1862902109277343, 0.1255803694649046, 0.0556685671161737};

    for (auto i = 0; i < rule_1d; ++i)
      for (auto j = 0; j < rule_1d; ++j)
      {
        qp[i*rule_1d+j] = SPoint(coordsGauss[i],coordsGauss[j]);
        qweights[i*rule_1d+j] = weightsGauss[i]*weightsGauss[j];
      }

  } else {
    
    std::cerr << "Invalid quad rule specified: quad rule = " << rule_1d << '\n';
    std::cerr << "Only rules between 1 and 11 allowed" << std::endl;
    exit(1);
    
  }
}
                          
}
// end namespace FEQuadrature_Felippa

}
// end namespace hydrofem
