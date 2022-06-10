// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_Linear_HPP__
#define __Hydrofem_Equation_Linear_HPP__

#include "Hydrofem_Equation.hpp"

namespace hydrofem
{

/**
 * \brief class describing linear advection equation
 * 
 */
template <typename ScalarT>
class Equation_Linear
  :
  public Equation<ScalarT> 
{
  
public:
  
  using LVec = typename Equation<ScalarT>::LVec;
  using LMat = typename Equation<ScalarT>::LMat;
  
  //! \brief Ctor
  Equation_Linear()
  {
    this->nEq = 1;
    this->ndims = 2;
    this->m_field_names = {"u"};
  }

  //! \brief Ctor from a given characteristic velocity field  
  Equation_Linear(const std::shared_ptr<std::function<SPoint(SPoint)>>& velocity)
    : m_velocity(velocity)
  {
    this->nEq = 1;
    this->ndims = 2;
    this->m_field_names = {"u"};
  }
  
  //! \brief Dtor
  virtual ~Equation_Linear() {}
  
  //! \brief flux
  virtual LVec flux(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const;

  //! \brief flux jacobian
  virtual LMat flux_u(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const;
  
  //! \brief max eigenvalue
  virtual ScalarT lambda_max(const SPoint& x, const double t, const LVec& U, const SPoint& normal) const;
  
protected:

  //! \brief set the velocity field 
  inline void set_velocity(const std::shared_ptr<std::function<SPoint(SPoint)>>& velocity)
  { m_velocity = velocity; }  
  
private:

  //! \brief the main definition of the velocity field
  std::shared_ptr<std::function<SPoint(SPoint)>> m_velocity;  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_Linear_HPP__ */
