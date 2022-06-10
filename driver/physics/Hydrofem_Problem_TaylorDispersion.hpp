// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_TaylorDispersion_HPP__
#define __Hydrofem_Problem_TaylorDispersion_HPP__

#include "Hydrofem_Problem_HyperbolicSystem.hpp"

namespace hydrofem
{

/**
 * Taylor dispersion experiment problem as presented in:
 * O. Boyarkine et al./ JCP 230 (2011) 2896-2914
 * Convection dominant flow to show the performance of the
 * new monotone LPS scheme and the LED schemes.
 */
class Problem_TaylorDispersion
  :
  public Problem_HyperbolicSystem
{
public:

  // creates a standard LPS stabilized problem
  explicit Problem_TaylorDispersion(const std::string& name = "taylor-dispersion") : public Problem_HyperbolicSystem()
  {
    setName(name);
  }

  ~Problem_TaylorDispersion() override = default;
  
  void init() override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_TaylorDispersion_HPP__ */
