// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Stepper_HPP__
#define __Hydrofem_Stepper_HPP__

#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_OptionHandler.hpp"
#include "Hydrofem_Assembler_Base.hpp"


namespace hydrofem
{

/**
 * \brief A base time stepper class
 */
class Stepper
  :
  public Optionable
{
public:
  
  //! \brief Ctor
  explicit Stepper(const std::shared_ptr<OptionHandler>& option_handler,
                   const std::shared_ptr<Assembler_Base>& assembler) : Optionable(option_handler)
  {
    option_handler->parse();
    m_option_handler = option_handler;
    m_assembler = assembler;
  }

  //! \brief Dtor
  virtual ~Stepper() = default;
  
  //! \brief solve one time step
  virtual void solveStep() = 0;
  
  //! \brief get current time (helps when doing adaptive time stepping)
  [[nodiscard]] virtual double time() const { return m_time; }

  //! \brief current time step size
  [[nodiscard]] virtual double dt() const { return m_delta_t; }

  //! \brief initial time
  [[nodiscard]] virtual double t0() const { return m_t0; }

  //! \brief final time
  [[nodiscard]] virtual double tf() const { return m_tf; }

  //! @brief get the current solution (for printing to file)
  [[nodiscard]] virtual std::shared_ptr<const FEVector> getCurrentSolution() const = 0;
  
protected:

  virtual void addOptionsCallback(po::options_description& config)
  {
    config.add_options()
      ("stepper-t0",po::value<double>(&m_t0)->default_value(0.0),"Initial time")
      ("stepper-tf",po::value<double>(&m_tf)->default_value(0.0),"Final time")
      ("stepper-dt",po::value<double>(&m_delta_t)->default_value(NAN),"Time step");
  }

  //! \brief routine for computing delta t from CFL
  [[maybe_unused]] virtual void computeDeltaT() const {}

  //! common stepper parameters 
  //@{
  //! initial time
  double                           m_t0;
  //! final time
  double                           m_tf;
  //! delta t
  mutable double                   m_delta_t;
  //! current time
  mutable double                   m_time;
  //! time step number
  mutable int                      m_stepnum = 0;
  //@}
  
  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;
  // tool for assembling residual and Jacobian
  std::shared_ptr<Assembler_Base> m_assembler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Stepper_HPP__ */
