// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FunctionElement_HPP__
#define __Hydrofem_FunctionElement_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{
  
template <typename ScalarT>
class FunctionElement
{
  
public:

  /// Coefficient type 
  using CoeffType = typename FEArray<ScalarT>::Node;

  /// Const Coefficient type 
  using CoeffTypeConst = typename FEArray<const ScalarT>::Node;
  
  //! \brief Ctor
  FunctionElement(const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis)
  {
    m_basis = basis;
    m_p = basis->at(0)->p();
    m_ndof = m_basis->size();
    m_c_e = createKArray<CoeffType>(m_ndof);
  }
  
  //! \brief Dtor
  virtual ~FunctionElement() {}
  
  /// The basis order
  virtual inline int p() const
  { return m_p; }
  
  /// set coefficients
  virtual void set_coefficients(const CoeffType& c_e);
  
  /// return coefficients
  virtual const CoeffType& coeffs() const;
  
  virtual CoeffType& coeffs();
  
  /// return number of degrees of freedom
  virtual inline int ndof() const
  { return m_ndof; }
  
  /// eval function at cartesian coordinates
  virtual ScalarT operator()(const SPoint& P, const std::vector<SPoint>& element) const;

  /// eval function at cartesian coordinates
  virtual ScalarT operator()(const SPoint& P, const std::vector<SPoint>& element, const CoeffTypeConst& c_e) const;
  
  /// compute the gradient at cartesian coordinates
  virtual CoeffType grad(const SPoint&P, const std::vector<SPoint>& element) const;
  
  /// calculate error on element
  virtual ScalarT error_l1(const std::function<ScalarT(const SPoint&)>& func,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<SPoint>& element) const;
  
  /// calculate error on element
  virtual ScalarT error_l2(const std::function<ScalarT(const SPoint&)>& func,
                           const std::shared_ptr<Quadrature>& quadrature,
                           const std::vector<SPoint>& element) const;

protected:

  /// basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;

  /// Coefficients
  CoeffType m_c_e;
  
  /// Order
  int m_p;
  
  /// Num of local dofs
  int m_ndof;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_FunctionElement_HPP__ */
