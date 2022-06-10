// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IC_ReactiveFlowWithBlob_HPP__
#define __Hydrofem_IC_ReactiveFlowWithBlob_HPP__

#include "Hydrofem_InitialCondition.hpp"

namespace hydrofem
{

class IC_ReactiveFlowWithBlob
  :
  public InitialCondition
{
public:
  
  using LVec = InitialCondition::LVec;
  
  //! \brief Ctor
  IC_ReactiveFlowWithBlob() {}

  //! \brief Dtor
  virtual ~IC_ReactiveFlowWithBlob() {}
  
  //! \brief spatial evaluation call
  LVec evaluate(const SPoint&) const override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IC_ReactiveFlowWithBlob_HPP__ */
