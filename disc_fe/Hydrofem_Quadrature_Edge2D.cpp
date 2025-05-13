// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEQuadrature.hpp"
#include "Hydrofem_Quadrature_Edge2D.hpp"

namespace hydrofem
{

Quadrature_Edge2D::Quadrature_Edge2D(const int order, const std::vector<SPoint>& edge)
{
  assert(int(edge.size()) == 2);
  for (auto it : edge)
    assert(it.size()==static_cast<std::size_t>(2));
  buildQuadrature(order,edge);
}

void Quadrature_Edge2D::buildQuadrature(const int order, const std::vector<SPoint>& edge)
{
  using namespace FEQuadrature;
  
  const SPoint& point1 = edge.at(0);
  const SPoint& point2 = edge.at(1);

  const auto& rule_ = order;
  auto& _q_points = get_q_pointsView();
  auto& _q_weights = get_q_weightsView();
  
  _q_points.clear();
  _q_weights.clear();
  
  if (rule_ == 1)
  {
    
    SPoint AB = point2 - point1;
    double length = AB.norm();
    _q_points.push_back(0.5*(point1 + point2));
    _q_weights.push_back(length);
    
  } else if (rule_ == 2) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 2; ++k)
    {
      x = cx + coordsGauss1D3[k]*hx;
      y = cy + coordsGauss1D3[k]*hy;      
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D3[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 3) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 3; ++k)
    {
      x = cx + coordsGauss1D5[k]*hx;
      y = cy + coordsGauss1D5[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D5[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 4) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 4; ++k)
    {
      x = cx + coordsGauss1D7[k]*hx;
      y = cy + coordsGauss1D7[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D7[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 5) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 5; ++k)
    {
      x = cx + coordsGauss1D9[k]*hx;
      y = cy + coordsGauss1D9[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D9[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 6) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 6; ++k)
    {
      x = cx + coordsGauss1D11[k]*hx;
      y = cy + coordsGauss1D11[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D11[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 7) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 7; ++k)
    {
      x = cx + coordsGauss1D13[k]*hx;
      y = cy + coordsGauss1D13[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D13[k]*sqrt(hx*hx+hy*hy));
    }

  } else if (rule_ == 8) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 8; ++k)
    {
      x = cx + coordsGauss1D15[k]*hx;
      y = cy + coordsGauss1D15[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D15[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 9) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 9; ++k)
    {
      x = cx + coordsGauss1D17[k]*hx;
      y = cy + coordsGauss1D17[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D17[k]*sqrt(hx*hx+hy*hy));
    }

  } else if (rule_ == 10) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 10; ++k)
    {
      x = cx + coordsGauss1D19[k]*hx;
      y = cy + coordsGauss1D19[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D19[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else if (rule_ == 11) {
    
    double hx = 0.5*(point2.x()-point1.x());
    double hy = 0.5*(point2.y()-point1.y());
    double cx = 0.5*(point2.x()+point1.x());
    double cy = 0.5*(point2.y()+point1.y());
    double x, y;
    for (int k = 0; k < 11; ++k)
    {
      x = cx + coordsGauss1D21[k]*hx;
      y = cy + coordsGauss1D21[k]*hy;
      _q_points.push_back(SPoint(x,y));
      _q_weights.push_back(weightsGauss1D21[k]*sqrt(hx*hx+hy*hy));
    }
    
  } else {
    
    std::cerr << "Invaling Gauss quadrature rule with " << rule_ << " points, chosen.\n";
    std::cerr << "Only rules up to " << 11 << " points allowed." << std::endl;
    exit(1);
    
  }
}

}
// end namespace hydrofem
