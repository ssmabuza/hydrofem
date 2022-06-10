// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_FEQuadrature_HPP__
#define __Hydrofem_FEQuadrature_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

namespace FEQuadrature
{

// 1D Gaussian Quadrature Weights and Abscissae from: https://pomax.github.io/bezierinfo/legendre-gauss.html
//@{

/// 2 point Gauss in 1D, exact on degree 3 polys
const double coordsGauss1D3[2] = {-0.57735026918962576451, 0.57735026918962576451};
const double weightsGauss1D3[2] = {1.0000000000000000, 1.0000000000000000};

/// 3 point Gauss in 1D, exact on degree 5 polys
const double coordsGauss1D5[3] = {-0.7745966692414837704, 0.0000000000000000000, 0.7745966692414837704};
const double weightsGauss1D5[3] = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};

/// 4 point Gauss in 1D, exact on degree 7 polys
const double coordsGauss1D7[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
const double weightsGauss1D7[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};

/// 5 point Gauss in 1D, exact on degree 9 polys
const double coordsGauss1D9[5] = {-0.9061798459386640, -0.5384693101056831, 0.0000000000000000, 0.5384693101056831, 0.9061798459386640};
const double weightsGauss1D9[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

/// 6 point Gauss in 1D, exact on degree 11 polys
const double coordsGauss1D11[6] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
const double weightsGauss1D11[6] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704};

/// 7 point Gauss in 1D, exact on degree 13 polys
const double coordsGauss1D13[7] = {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0000000000000000, 0.4058451513773972, 0.7415311855993945, 0.9491079123427585};
const double weightsGauss1D13[7] = {0.1294849661688697, 0.2797053914892766, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.2797053914892766, 0.1294849661688697};

/// 8 point Gauss in 1D, exact on degree 15 polys
const double coordsGauss1D15[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498, 0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
const double weightsGauss1D15[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620, 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};

/// 9 point Gauss in 1D, exact on degree 17 polys 
const double coordsGauss1D17[9] = {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0000000000000000, 0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261};
const double weightsGauss1D17[9] = {0.3123470770400029, 0.1806481606948574, 0.2606106964029354, 0.3123470770400029, 0.3302393550012598, 0.3123470770400029, 0.2606106964029354, 0.1806481606948574, 0.3123470770400029};

/// 10 point Gauss in 1D, exact on degree 19 polys
const double coordsGauss1D19[10] = {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717};
const double weightsGauss1D19[10] = {0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963, 0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820, 0.1494513491505806, 0.0666713443086881};

/// 11 point Gauss in 1D, exact on degree 21 polys
const double coordsGauss1D21[11] = {-0.9782286581460570, -0.8870625997680953, -0.7301520055740494, -0.5190961292068118, -0.2695431559523450, 0.0000000000000000, 0.2695431559523450, 0.5190961292068118, 0.7301520055740494, 0.8870625997680953, 0.9782286581460570};
const double weightsGauss1D21[11] = {0.0556685671161737, 0.1255803694649046, 0.1862902109277343, 0.2331937645919905, 0.2628045445102467, 0.2729250867779006, 0.2628045445102467, 0.2331937645919905, 0.1862902109277343, 0.1255803694649046, 0.0556685671161737};

//@}


// 2D Gaussian Quadrature Weights and Abscissae
//@{

/// Gauss quadrature in 2D triangles: exact on quadratic polynomials
const double oneThirdTri2 = 1.0/3.0;
const double coordsTri2[9] = {0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0};
const double weightsTri2[3] = {oneThirdTri2, oneThirdTri2, oneThirdTri2};

/// Gauss quadrature in 2D triangles: exact on polynomials of degree 5
const double beta0Tri5 = (6.0 + sqrt(15.0))/21.0;
const double beta1Tri5 = (6.0 - sqrt(15.0))/21.0;
const double alpha0Tri5 = 1.0 - 2.0*beta0Tri5;
const double alpha1Tri5 = 1.0 - 2.0*beta1Tri5;
const double oneThirdTri5 = 1.0/3.0;
const double coordsTri5[21] = {oneThirdTri5, oneThirdTri5, oneThirdTri5,
                               alpha0Tri5, beta0Tri5, beta0Tri5,
                               beta0Tri5, alpha0Tri5, beta0Tri5,
                               beta0Tri5, beta0Tri5, alpha0Tri5,
                               alpha1Tri5, beta1Tri5, beta1Tri5,
                               beta1Tri5, alpha1Tri5, beta1Tri5,
                               beta1Tri5, beta1Tri5, alpha1Tri5};
const double aux0Tri5 = (155.0 + sqrt(15.0))/1200.0;
const double aux1Tri5 = (155.0 - sqrt(15.0))/1200.0;
const double weightsTri5[7] = {0.225, aux0Tri5, aux0Tri5, aux0Tri5, aux1Tri5, aux1Tri5, aux1Tri5 };

/// Gauss quadrature in 2D quadrilaterals: exact on quadratic polynomials
const double betaTet2 = (5.0 - sqrt(5.0) )/20.0;
const double alphaTet2 = 1.0 - 3.0*betaTet2;
const double coordsTet2[16] = {alphaTet2, betaTet2, betaTet2, betaTet2,
                               betaTet2, alphaTet2, betaTet2, betaTet2,
                               betaTet2, betaTet2, alphaTet2, betaTet2,
                               betaTet2, betaTet2, betaTet2, alphaTet2};
const double weightsTet2[4] = { 0.25, 0.25, 0.25, 0.25 };

/// Gauss quadrature in 2D quadrilaterals: exact on cubic polynomials
const double oneSixthTet3 = 1.0/6.0;
const double coordsTet3[20] = {0.25, 0.25, 0.25, 0.25,
                               0.5, oneSixthTet3, oneSixthTet3, oneSixthTet3,
                               oneSixthTet3, 0.5, oneSixthTet3, oneSixthTet3,
                               oneSixthTet3, oneSixthTet3, 0.5, oneSixthTet3,
                               oneSixthTet3, oneSixthTet3, oneSixthTet3, 0.5};
const double weightsTet3[5] = {-0.80, 0.45, 0.45, 0.45, 0.45};

//@}

/// 1-point Gauss rule, exact on degree 1 polys
//@{
double CalculateIntegralOverEdge1(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun);

double CalculateIntegralOverEdge1(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers);
//@}
                                  
/// 2-point Gauss rule, exact on degree 3 polys
//@{
double CalculateIntegralOverEdge3(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun);

double CalculateIntegralOverEdge3(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers);
//@}

/// 3-point Gauss rule, exact on degree 5 polys
//@{
double CalculateIntegralOverEdge5(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::function<double(double,double)>& fun);

double CalculateIntegralOverEdge5(const SPoint& point1, 
                                  const SPoint& point2, 
                                  const std::vector<std::function<double(double,double)>>& funs,
                                  const std::vector<double>& multipliers);                                  
//@}
                                  
/// Gauss quadrature, N-point rule
//@{
double CalculateIntegralOverEdge(const int Rule, 
                                 const SPoint& point1, 
                                 const SPoint& point2, 
                                 const std::function<double(double,double)>& fun);

double CalculateIntegralOverEdge(const int Rule, 
                                 const SPoint& point1, 
                                 const SPoint& point2, 
                                 const std::vector<std::function<double(double,double)>>& funs,
                                 const std::vector<double>& multipliers);
//@}
                                  
/// First moment of the function
//@{
double CalculateFirstMoment5(const SPoint& point1, 
                             const SPoint& point2, 
                             const std::function<double(double,double)>& fun);

double CalculateFirstMoment5(const SPoint& point1, 
                             const SPoint& point2, 
                             const std::vector<std::function<double(double,double)>>& funs,
                             const std::vector<double>& multipliers);
//@}
                             
/// Exact on linears, central point rule
//@{
double CalculateIntegralOverTriangle1(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun);

double CalculateIntegralOverTriangle1(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers);
//@}
                                      
/// Exact on quadratics 3-point rule
//@{
double CalculateIntegralOverTriangle2(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun);

double CalculateIntegralOverTriangle2(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers);
//@}
                                      
/// Exact on qintics, 7-point rule
//@{
double CalculateIntegralOverTriangle5(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3, 
                                      const std::function<double(double,double)>& fun);

double CalculateIntegralOverTriangle5(const SPoint& point1, 
                                      const SPoint& point2, 
                                      const SPoint& point3,
                                      const std::vector<std::function<double(double,double)>>& funs,
                                      const std::vector<double>& multipliers);
//@}

}
//end namespace FEQuadrature

}
// end namespace hydrofem

#endif /** __Hydrofem_FEQuadrature_HPP__ */
