// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Bernstein_Basis_Triangle_HPP__
#define __Hydrofem_Bernstein_Basis_Triangle_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"

namespace hydrofem
{

/// evaluate basis function
class Bernstein_Basis_Triangle
  : 
  public FEBasis
{
  
  int p_;
  int m_i;
  int m_j;
  
public:
  
  Bernstein_Basis_Triangle(const int order,
                           const int i,
                           const int j) : p_(order), 
                                          m_i(i), m_j(j) {}
                                                      
  virtual ~Bernstein_Basis_Triangle() {}
  
  virtual inline int p() const override { return p_; }

  virtual double operator()(const SPoint& point, const std::vector<SPoint>& element) const override
  { return eval(point.x(),point.y(),element.at(0),element.at(1),element.at(2),m_i,m_j,p_ - m_i - m_j); }
  
  virtual SPoint grad(const SPoint& point, const std::vector<SPoint>& element) const override
  { return grad(point.x(),point.y(),element.at(0),element.at(1),element.at(2),m_i,m_j,p_ - m_i - m_j); }
  
  using FEBasisPtr = std::shared_ptr<FEBasis>;
  static std::shared_ptr<std::vector<FEBasisPtr>>
  buildElementBasis(const int order)
  {
    std::shared_ptr<std::vector<FEBasisPtr>> basis = std::make_shared<std::vector<FEBasisPtr>>();
    basis->clear();
    for (int i = 0; i <= order; ++i)
      for (int j = 0; j <= order - i; ++j)
        basis->push_back(std::make_shared<Bernstein_Basis_Triangle>(order,i,j));
    return basis;
  }
  
private:
  
  
  /// evaluation on reference triangle
  double eval(const double x,
              const double y,
              const int i,
              const int j,
              const int k) const;

  /// evaluation on arbitrary triangle using barycentric coordinates
  double eval(const double u,
              const double v,
              const double w,
              const int i,
              const int j,
              const int k) const;
  
  /// evaluation on arbitrary triangle using barycentric coordinates
  double eval(const double x,
              const double y,
              const SPoint& tri1,
              const SPoint& tri2,
              const SPoint& tri3,
              const int i,
              const int j,
              const int k) const;
  
  double eval_u(const double u,
                const double v,
                const double w,
                const int i,
                const int j,
                const int k) const;
  
  double eval_v(const double u,
                const double v,
                const double w,
                const int i,
                const int j,
                const int k) const;
  
  double eval_w(const double u,
                const double v,
                const double w,
                const int i,
                const int j,
                const int k) const;
  
  double eval_x(const double x,
                const double y, 
                const SPoint& tri1,
                const SPoint& tri2,
                const SPoint& tri3, 
                const int i,
                const int j,
                const int k) const;
  
  double eval_y(const double x,
                const double y, 
                const SPoint& tri1,
                const SPoint& tri2,
                const SPoint& tri3, 
                const int i,
                const int j,
                const int k) const;
  
  SPoint grad(const double x,
                const double y, 
                const SPoint& tri1,
                const SPoint& tri2,
                const SPoint& tri3, 
                const int i,
                const int j,
                const int k) const;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Bernstein_Basis_Triangle_HPP__ */
