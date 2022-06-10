// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_Bioseparation.hpp"

#include "Hydrofem_BC_Bioseparation.hpp"
#include "Hydrofem_IC_Bioseparation.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{
  
void Problem_Bioseparation::init()
{
  auto ic = std::make_shared<IC_Bioseparation>();
  assert(ic);
  this->set_ic(ic);
  auto exact_ = std::make_shared<AnalyticalExpression>();
  assert(exact_);
  this->set_exact(exact_);
  auto u_in_ = std::make_shared<AnalyticalExpression>();
  assert(u_in_);
  this->set_Uin(u_in_);
  auto bc = std::make_shared<BC_Scalar>();
  assert(bc);
  bc->setEquation(m_equation);
  this->set_bc(bc);
}
/// end void init()

}
// end namespace hydrofem
