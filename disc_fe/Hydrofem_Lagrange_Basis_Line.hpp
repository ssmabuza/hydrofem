// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Lagrange_Basis_Line_HPP__
#define __Hydrofem_Lagrange_Basis_Line_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_FEBasis.hpp"

namespace hydrofem
{

class Lagrange_Basis_Line
  :
  public FEBasis
{
  
  int m_p;
  int m_k;

  /// evaluation on arbitrary interval
  double eval(const double x,
              const SPoint& line1,
              const SPoint& line2,
              const int k) const;

  /// derivative on arbitrary interval              
  double eval_x(const double x,
                const SPoint& line1,
                const SPoint& line2,
                const int k) const;
                
public:
  
  Lagrange_Basis_Line(const int order,
                       const int k) : m_p(order), m_k(k) {}
                       
  virtual ~Lagrange_Basis_Line() {}
  
  virtual inline int p() const override { return m_p; }
  
  inline virtual double operator()(const SPoint& point, const std::vector<SPoint>& element) const override
  { return eval(point.x(),element.at(0),element.at(1), m_k); }
  
  inline virtual SPoint grad(const SPoint& point, const std::vector<SPoint>& element) const override
  { return SPoint(eval_x(point.x(),element.at(0),element.at(1), m_k)); }
  
  using FEBasisPtr = FEBasis::FEBasisPtr;
  static std::shared_ptr<std::vector<FEBasisPtr>>
  buildElementBasis(const int order)
  {
    std::shared_ptr<std::vector<FEBasisPtr>> basis = std::make_shared<std::vector<FEBasisPtr>>();
    for (int i = 0; i <= order; ++i)
      basis->push_back(std::make_shared<Lagrange_Basis_Line>(order,i));
    return basis;
  }
  
  template <typename ScalarT>
  inline ScalarT evaluateFunction(const ScalarT x, const int k) const 
  {
    assert((k <= m_p)&&(k >= 0));
    ScalarT h = 1.0/m_p;
    ScalarT num = 1.0;
    ScalarT denom = 1.0;
    const int& i = k;
    for (int j = 0; j <= m_p; ++j)
      if (j != i) {
        num *= (x - j*h);
        denom *= (i - j)*h;
      }
    return num/denom;
  }
  
  /// evaluation on Bezier interval
  double eval(const double x,
              const int k) const;

  /// derivative on Bezier interval
  double eval_x(const double x,
                const int k) const;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Lagrange_Basis_Line_HPP__ */
