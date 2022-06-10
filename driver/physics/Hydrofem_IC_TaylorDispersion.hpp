// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IC_TaylorDispersion_HPP__
#define __Hydrofem_IC_TaylorDispersion_HPP__

#include "Hydrofem_InitialCondition.hpp"

namespace hydrofem
{

/**
 * @brief The concrete implementation of IC for taylor dispersion
 */
class IC_TaylorDispersion
  :
  public InitialCondition
{
public:
  
  using LVec = InitialCondition::LVec;
  
  //! @brief Ctor
  IC_TaylorDispersion();
  
  //! @brief spatial evaluation call
  LVec evaluate(const SPoint&) const override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IC_TaylorDispersion_HPP__ */
