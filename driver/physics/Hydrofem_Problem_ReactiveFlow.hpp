// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Problem_ReactiveFlow_HPP__
#define __Hydrofem_Problem_ReactiveFlow_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

/**
 * \brief The base class that describes a continuous problem
 */
class Problem_ReactiveFlow
  :
  public Problem
{
public:
  
  explicit Problem_ReactiveFlow(const std::shared_ptr<OptionHandler>& opt_handler)
  { 
    m_name = "Reactive Flow";
    m_dof_names = {"conc_f","conc_w"};
  }
  
  ~Problem_ReactiveFlow() override = default;
  
  void init() override;
  
  //! \brief get functions
  //@{
  
  //! @brief The gives an analytical fluid velocity
  [[nodiscard]]
  std::shared_ptr<std::function<SPoint(SPoint)>>
  getFluidVelocity() const
  { return m_fluid_velocity; }
  
  //! @brief get the desorption rate
  [[nodiscard]]
  double getDesorptionRate() const
  { return k_d; }

  //! @brief get the diffusion coefficient
  [[nodiscard]]
  double getDiffusionCoefficient() const
  { return m_diff; }

  //! @brief get the adsorption isotherm
  [[nodiscard]]
  std::shared_ptr<std::function<double(double)>>
  getAdsorptionIsotherm() const
  { return m_ads; }
  
  //! @brief get the inflow concentration
  [[nodiscard]]
  std::shared_ptr<std::function<double(SPoint,double)>>
  getInflowConcentration() const
  { return m_inflow; }
  //@}

private:

  // parameter list
  std::shared_ptr<OptionHandler> m_opt_handler;
  // analytical fluid velocity function v(x)
  std::shared_ptr<std::function<SPoint(SPoint)>> m_fluid_velocity;
  // adsorption isotherm R(c)
  std::shared_ptr<std::function<double(double)>> m_ads;
  // desorption constant
  double k_d;
  // inflow function c_in(x,t)
  std::shared_ptr<std::function<double(SPoint,double)>> m_inflow = Teuchos::null;
  // diffusion coefficient
  double m_diff;
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_ReactiveFlow_HPP__ */
