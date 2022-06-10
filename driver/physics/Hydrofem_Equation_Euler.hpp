// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_Euler_HPP__
#define __Hydrofem_Equation_Euler_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Equation.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

/**
 * \brief A class for the Euler equations
 */
template <typename ScalarT>
class Equation_Euler
  :
  public Equation<ScalarT>
{
public:

  using LVec = typename Equation<ScalarT>::LVec;
  using LMat = typename Equation<ScalarT>::LMat;
  
  //! \brief Empty Ctor : default value for gamma is 1.4 (air)
  Equation_Euler() : m_gamma(1.4)
  {
    this->ndims = 1;
    this->nEq = this->ndims+2;
  }
  
  //! \brief Ctor from adiabatic constant
  Equation_Euler(const double gamma, const int dims) : m_gamma(gamma)
  {
    this->ndims = dims;
    this->nEq = this->ndims+2;
    if (this->ndims==1)
      this->m_field_names = {"rho","rhou","rhoE"};
    else if (this->ndims==2)
      this->m_field_names = {"rho","rhou","rhov","rhoE"};
    else if (this->ndims==3)
      this->m_field_names = {"rho","rhou","rhov","rhow","rhoE"};
  }
  
  //! \brief Dtor
  virtual ~Equation_Euler() {}
  
  //! \brief 
  virtual LVec flux(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief 
  virtual LMat flux_u(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief max eigenvalue
  virtual ScalarT lambda_max(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief eigenvalues 
  virtual LVec lambda(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief eigenvector matrix wrt conserved variables
  virtual LMat eigenvector(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief inverse eigenvector matrix wrt conserved variables
  virtual LMat inverse_eigenvector(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const override;
  
  //! \brief get the value for gamma
  inline double gamma() const { return m_gamma; }
  
  //! \brief get the pressure
  inline ScalarT p(const LVec& U) const { return internal::pressure(U,m_gamma,this->ndims); }
  
  inline ScalarT c(const LVec& U) const { return internal::c(U,m_gamma,this->ndims); }

  inline ScalarT S(const LVec& U) const { return internal::S(U,m_gamma,this->ndims); }
  
  using Equation<ScalarT>::ns;
  using Equation<ScalarT>::dim;
  
  
private:

  /** \brief The internal implementation of Euler quantities */
  class internal
  {
  public:
    
    static ScalarT pressure(const LVec& U, const double gamma, const int nDims);
    
    //! \brief Physical flux wrt. conserved variables in the direction of normal
    static LVec flux(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    //! \brief Eigenvalues wrt conserved variables
    static LVec lambda(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    //! \brief Eigenvalue matrix wrt conserved variables
    static LMat lambda_mat(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    //! \brief Maximum eigenvalue wrt conserved variables
    static ScalarT lambda_max(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    static LMat eigenvector(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    static LMat inverse_eigenvector(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
    static ScalarT c(const LVec& U, const double gamma, const int nDims);
    
    static ScalarT H(const LVec& U, const double gamma, const int nDims);
    
    static ScalarT S(const LVec& U, const double gamma, const int nDims);
    
    static void conservedtoroe(LVec& V, const LVec& U, const double gamma, const int nDims);
    
    static void roe_mean(LVec& Ulr, const LVec& Ul, const LVec& Ur, const double gamma, const int nDims);
    
    static LMat abs_flux_u(const LVec& U, const SPoint& normal, const double gamma, const int nDims);
    
  };
  
  //! adiabatic constant
  double m_gamma;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_Euler_HPP__ */
