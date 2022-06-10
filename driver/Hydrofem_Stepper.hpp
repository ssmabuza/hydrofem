// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Stepper_HPP__
#define __Hydrofem_Stepper_HPP__

#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_OptionHandler.hpp"


namespace hydrofem
{

/**
 * @brief A base time stepper class
 */
class Stepper
  :
  public Optionable
{
public:
  
  //! @brief Ctor
  explicit Stepper(const std::shared_ptr<OptionHandler>& option_handler)
  {
    m_option_handler = option_handler;
  }

  //! @brief Dtor
  virtual ~Stepper() = default;
  
  //! @brief solve one time step
  virtual void solveStep() = 0;
  
  //! @brief get current time (helps when doing adaptive time stepping)
  [[nodiscard]] virtual double time() const = 0;
  
  //! @brief get the current solution (for printing to file)
  [[nodiscard]] virtual std::shared_ptr<FEVector> getCurrentSolution() const = 0;
  
protected:

  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;
  

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Stepper_HPP__ */
