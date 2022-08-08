// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_Bioseparation.hpp"

#include "Hydrofem_BC_Bioseparation.hpp"
#include "Hydrofem_InitialCondition.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{
  
void Problem_Bioseparation::init()
{
  m_ic = std::make_shared<ConstantScalarInitialCondition>(0.0);
  assert(ic);
  auto u_in_ = std::make_shared<AnalyticalExpression>();
  assert(u_in_);
  this->set_Uin(u_in_);
  auto bc = std::make_shared<BC_Scalar>();
  assert(bc);
  bc->setEquation(m_equation);
  this->set_bc(bc);
}

void Problem_Bioseparation::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("prob-biosep-omega",po::value<double>(&m_omega)->default_value(0.84),"Porosity value")
    ("prob-biosep-rho_s",po::value<double>(&m_rho_s)->default_value(1.0),"Density of membrane")
    ("prob-biosep-q_max",po::value<double>(&m_q_max)->default_value(150.0),"Maximum binding capacity")
    ("prob-biosep-K_eq",po::value<double>(&m_K_eq)->default_value(2.06),"Langmuir equilibrium constant");
}

}
// end namespace hydrofem
