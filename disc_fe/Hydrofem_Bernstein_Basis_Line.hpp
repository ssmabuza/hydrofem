// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Bernstein_Basis_Line_HPP__
#define __Hydrofem_Bernstein_Basis_Line_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_FEBasis.hpp"

namespace hydrofem
{

class Bernstein_Basis_Line
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
  
  Bernstein_Basis_Line(const int order,
                       const int k) : m_p(order), m_k(k) {}
                       
  virtual ~Bernstein_Basis_Line() {}
  
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
      basis->push_back(std::make_shared<Bernstein_Basis_Line>(order,i));
    return basis;
  }
  
  template <typename ScalarT>
  inline ScalarT evaluateFunction(const ScalarT x, const int k) const 
  {
    assert((k <= m_p)&&(k >= 0));
    return static_cast<ScalarT>(Intern::factorial(m_p)/(Intern::factorial(k)*Intern::factorial(m_p-k)))
                   * std::pow(x,ScalarT(k)) * std::pow(static_cast<ScalarT>(1.0)-x,ScalarT(m_p-k));
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

#endif /** __Hydrofem_Bernstein_Basis_Line_HPP__ */
