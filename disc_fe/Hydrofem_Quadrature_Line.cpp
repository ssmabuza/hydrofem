// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEQuadrature.hpp"
#include "Hydrofem_Quadrature_Line.hpp"

namespace hydrofem
{

Quadrature_Line::Quadrature_Line(const int order, const std::vector<SPoint>& line)
{
  assert(int(line.size()) == 2);
  for (auto it : line)
    assert(it.size()==static_cast<std::size_t>(1));
  buildQuadrature(order,line);
}

void Quadrature_Line::buildQuadrature(const int order, const std::vector<SPoint>& line)
{
  using namespace FEQuadrature;

  const auto& rule_ = order;
  const SPoint& line1 = line.at(0);
  const SPoint& line2 = line.at(1);
  
  auto& _q_points = get_q_pointsView();
  auto& _q_weights = get_q_weightsView();
  
  _q_points.clear();
  _q_weights.clear();
  
  if (rule_ == 1)
  {
    
    SPoint AB = line2 - line1;
    double length = AB.norm();
    _q_points.push_back(0.5*(line1 + line2));
    _q_weights.push_back(length);
    
  } else if (rule_ == 2) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 2; ++k)
    {
      x = cx + coordsGauss1D3[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D3[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 3) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 3; ++k)
    {
      x = cx + coordsGauss1D5[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D5[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 4) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 4; ++k)
    {
      x = cx + coordsGauss1D7[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D7[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 5) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 5; ++k)
    {
      x = cx + coordsGauss1D9[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D9[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 6) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 6; ++k)
    {
      x = cx + coordsGauss1D11[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D11[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 7) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 7; ++k)
    {
      x = cx + coordsGauss1D13[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D13[k]*sqrt(hx*hx));
    }

  } else if (rule_ == 8) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 8; ++k)
    {
      x = cx + coordsGauss1D15[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D15[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 9) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 9; ++k)
    {
      x = cx + coordsGauss1D17[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D17[k]*sqrt(hx*hx));
    }

  } else if (rule_ == 10) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 10; ++k)
    {
      x = cx + coordsGauss1D19[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D19[k]*sqrt(hx*hx));
    }
    
  } else if (rule_ == 11) {
    
    double hx = 0.5*(line2.x()-line1.x());
    double cx = 0.5*(line2.x()+line1.x());
    double x;
    for (int k = 0; k < 11; ++k)
    {
      x = cx + coordsGauss1D21[k]*hx;
      _q_points.push_back(SPoint(x));
      _q_weights.push_back(weightsGauss1D21[k]*sqrt(hx*hx));
    }
    
  } else {
    
    std::cerr << "Invaling Gauss quadrature rule with " << rule_ << " points, chosen.\n";
    std::cerr << "Only rules up to " << 11 << " points allowed." << std::endl;
    exit(1);
    
  }
}

}
// end namespace hydrofem
