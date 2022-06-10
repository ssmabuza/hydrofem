// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_Linear_impl_HPP__
#define __Hydrofem_Equation_Linear_impl_HPP__

namespace hydrofem
{

template <typename ScalarT>
typename Equation_Linear<ScalarT>::LVec 
Equation_Linear<ScalarT>::
flux(const SPoint& x, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  LVec res = createKArray<LVec>(1);
  res(0) = U(0) * ((*m_velocity)(x) * normal);
  return res;
}

template <typename ScalarT>
typename Equation_Linear<ScalarT>::LMat 
Equation_Linear<ScalarT>::
flux_u(const SPoint& x, const double /*t*/, const LVec& /*U*/, const SPoint& normal) const
{
  LMat res = createKArray<LMat>(1,1);
  res(0,0) = (*m_velocity)(x) * normal;
  return res;
}

template <typename ScalarT>
ScalarT Equation_Linear<ScalarT>::
lambda_max(const SPoint& x, const double /*t*/, const LVec& /*U*/, const SPoint& normal) const
{
  return static_cast<ScalarT>(std::fabs(normal * (*m_velocity)(x)));
}

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_Linear_impl_HPP__ */
