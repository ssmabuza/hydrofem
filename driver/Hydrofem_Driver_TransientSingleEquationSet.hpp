// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Driver_TransientSingleEquationSet_HPP__
#define __Hydrofem_Driver_TransientSingleEquationSet_HPP__

#include "Hydrofem_Driver.hpp"
#include "Hydrofem_OptionHandler.hpp"


namespace hydrofem
{

class Driver_TransientSingleEquationSet
  :
  public Driver
{
public:

  explicit Driver_TransientSingleEquationSet(const std::shared_ptr<OptionHandler>& option_handler) : Driver(option_handler) {}
  
  virtual ~Driver_TransientSingleEquationSet() {}

  void setup() {}
  
  void solve() {}
  
private:

  /** @brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config)
  {
    // nothing to parse
  }
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Driver_TransientSingleEquationSet_HPP__ */
