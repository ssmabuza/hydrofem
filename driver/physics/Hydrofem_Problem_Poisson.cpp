// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_Poisson.hpp"

#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{
  
const static std::function<double(SPoint)> exact_function = 
[](SPoint x)->double { return std::sin(2.0*M_PI*x.x())*std::sin(2.0*M_PI*x.y()); };

const static std::function<double(SPoint)> rhs_function =
[](SPoint x)->double { return - 8.0*M_PI*M_PI*std::sin(2.0*M_PI*x.x())*std::sin(2.0*M_PI*x.y()); };


class Exact_Function
  :
  public ScalarAnalyticalExpression
{
public:
  
  double evaluate(const hydrofem::SPoint &x) override
  { assert(!m_is_transient); return exact_function(x); }
  
  double operator()(const hydrofem::SPoint &x) override
  { return this->evaluate(x); }
  
};

class RHS_Function
  :
  ScalarAnalyticalExpression
{

  double evaluate(const hydrofem::SPoint &x) override
  { assert(!m_is_transient); return rhs_function(x); }

  double operator()(const hydrofem::SPoint &x) override
  { return this->evaluate(x); }

};
  
void Problem_Poisson::init()
{
  auto exact_ = std::make_shared<Exact_Function>(); assert(exact_);
  this->setExactSolution(exact_);

  m_bc = std::make_shared<BC_Scalar>();
  std::dynamic_pointer_cast<BC_Scalar>(m_bc)->initializeBoundaryPointsToDirichletEverywhere();

  m_dirichlet_bc_fnc = exact_;
}
/// end void init()

}
// end namespace hydrofem
