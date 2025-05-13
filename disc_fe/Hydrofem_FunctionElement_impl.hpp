// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#pragma once

#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_FunctionElement.hpp"

namespace hydrofem
{

template <typename ScalarT>
void FunctionElement<ScalarT>::
set_coefficients(const FunctionElement::CoeffType& c_e)
{
  assert(int(c_e.size()) == m_ndof);
  for (Eigen::Index i = 0; i < c_e.size(); ++i)
    m_c_e(i) = c_e(i);
}

template <typename ScalarT>
const typename FunctionElement<ScalarT>::CoeffType& 
FunctionElement<ScalarT>::coeffs() const
{
  return m_c_e;
}

template <typename ScalarT>
typename FunctionElement<ScalarT>::CoeffType& 
FunctionElement<ScalarT>::coeffs()
{
  return m_c_e;
}

template <typename ScalarT>
ScalarT FunctionElement<ScalarT>::
operator()(const SPoint& P, const std::vector<SPoint>& element) const
{
  ScalarT res = 0.0;
  for (Eigen::Index i = 0; i < m_c_e.size(); ++i)
    res += m_c_e[i]*(*(m_basis->at(i)))(P,element);
  return res;
}

template <typename ScalarT>
ScalarT FunctionElement<ScalarT>::
operator()(const SPoint& P, const std::vector<SPoint>& element, const FunctionElement<ScalarT>::CoeffTypeConst& c_e) const
{
  ScalarT res = 0.0;
  for (Eigen::Index i = 0; i < m_c_e.size(); ++i)
    res += c_e[i]*(*(m_basis->at(i)))(P,element);
  return res;
}

template <typename ScalarT>
typename FunctionElement<ScalarT>::CoeffType 
FunctionElement<ScalarT>::grad(const SPoint& P, const std::vector<SPoint>& element) const
{
  CoeffType res = createKArray<CoeffType>(P.size());
  for (Eigen::Index dim = 0; dim < res.dimension(0); ++dim)
    res(dim) = 0.0;
  for (Eigen::Index i = 0; i < m_c_e.size(); ++i)
  {
    auto grad_loc = (*(m_basis->at(i))).grad(P,element);
    for (Eigen::Index dim = 0; dim < res.size(); ++dim)
      res(dim) += m_c_e[i]*grad_loc(dim);
  }
  return res;
}

template <typename ScalarT>
ScalarT FunctionElement<ScalarT>::
error_l1(const std::function<ScalarT(const SPoint&)>& func,
         const std::shared_ptr<Quadrature>& quadrature,
         const std::vector<SPoint>& element) const
{
  const auto f = [&](const SPoint& P)->ScalarT
    {
      return std::abs(func(P)-this->operator()(P,element));
    };

  const auto& qweights = quadrature->get_q_weights();
  const auto& qpoints = quadrature->get_q_points();
  
  ScalarT res = 0.0;
  for (std::size_t qp = 0; qp < qpoints.size(); ++qp)
    res += qweights.at(qp) * f(qpoints.at(qp));
  return res;
}

template <typename ScalarT>
ScalarT FunctionElement<ScalarT>::
error_l2(const std::function<ScalarT(const SPoint&)>& func,
         const std::shared_ptr<Quadrature>& quadrature, const std::vector<SPoint>& element) const
{
  const auto f = [&](const SPoint& P)->ScalarT
    {
      return std::pow(func(P)-this->operator()(P,element),2.0);
    };

  const auto& qweights = quadrature->get_q_weights();
  const auto& qpoints = quadrature->get_q_points();
  
  ScalarT res = 0.0;
  for (std::size_t qp = 0; qp < qpoints.size(); ++qp)
    res += qweights.at(qp) * f(qpoints.at(qp));
  return res;
}

}
// end namespace hydrofem
