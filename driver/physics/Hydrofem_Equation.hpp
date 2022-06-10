// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_HPP__
#define __Hydrofem_Equation_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"

#define EQN_SAFETY_CHECK( neq, ndims ) \
  { \
    assert((neq)>0); \
    assert((ndims)>0); \
  }

namespace hydrofem
{
  
//! \brief base class for scalar hyperbolic equations
template <typename ScalarT>
class Equation
{
public:
  
  using LVec = LVEC_<ScalarT>;
  using LMat = LMAT_<ScalarT>;

  //! \brief number of dimensions
  virtual inline int dim() const
  { return ndims; }
  
  //! \brief size of the system of equations
  virtual inline int ns() const
  { return nEq; }
  
  //! \brief Ctor
  Equation() {}
  
  //! \brief Dtor
  virtual ~Equation() {}

  //! \brief flux
  virtual LVec flux(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LVec>(ns()); }

  //! \brief flux jacobian
  virtual LMat flux_u(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LMat>(ns(),ns()); }
  
  //! \brief diffusive flux
  virtual LVec flux_d(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const LMat& /*grad_U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LVec>(ns()); }
  
  //! \brief diffusive flux jacobian
  virtual LMat flux_d_u(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const LMat& /*grad_U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LMat>(ns(),ns()); }
  
  //! \brief max eigenvalue
  virtual ScalarT lambda_max(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const = 0;
  
  //! \brief eigenvalues 
  virtual LVec lambda(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LVec>(ns()); }
  
  //! \brief eigenvector matrix wrt conserved variables
  virtual LMat eigenvector(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LMat>(ns(),ns()); }
  
  //! \brief inverse eigenvector matrix wrt conserved variables
  virtual LMat inverse_eigenvector(const SPoint& /*x*/, const double /*t*/, const LVec& /*U*/, const SPoint& /*normal*/) const
  { EQN_SAFETY_CHECK(nEq,ndims); return createKArray<LMat>(ns(),ns()); }
  
  //! \brief get the field corresponding to eq 
  [[nodiscard]] std::string getFieldName(const int eq) const { return m_field_names.at(eq); }
  
protected:

  //! field names 
  std::vector<std::string> m_field_names;
  
  //! \brief size of the system of equations
  int nEq = -1;
  
  //! \brief number of dimensions
  int ndims = -1;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_HPP__ */
