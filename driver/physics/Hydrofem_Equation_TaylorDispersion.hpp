// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_TaylorDispersion_HPP__
#define __Hydrofem_Equation_TaylorDispersion_HPP__

#include "Hydrofem_Equation_Linear.hpp"

namespace hydrofem
{

/**
 * \brief the class for the linear advection equation with
 *       a Poiseuille velocity for advection as in the Taylor dispersion study
 */
template <typename ScalarT>
class Equation_TaylorDispersion
  :
  public Equation_Linear<ScalarT>
{
public:
  
  //! \brief Ctor
  Equation_TaylorDispersion(double V_max, double y_length);

  //! \brief Dtor
  virtual ~Equation_TaylorDispersion() {}
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_TaylorDispersion_HPP__ */
