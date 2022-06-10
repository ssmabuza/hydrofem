// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Problem_Bioseparation_HPP__
#define __Hydrofem_Problem_Bioseparation_HPP__

#include "Hydrofem_Problem.hpp"

namespace hydrofem
{

/**
 * Taylor dispersion experiment problem as presented in:
 * O. Boyarkine et al./ JCP 230 (2011) 2896-2914
 * Convection dominant flow to show the performance of the
 * new monotone LPS scheme and the LED schemes.
 */
class Problem_Bioseparation
  :
  public Problem
{
public:

  // creates a standard LPS stabilized problem
  explicit Problem_Bioseparation(const std::string& name = "bioseparation") : public Problem()
  {
    setName(name);
    setDofNames({{"conc"}});
  }

  ~Problem_Bioseparation() override = default;
  
  void init() override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_Bioseparation_HPP__ */
