// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_Poisson.hpp"

#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{

// exact solution
const static std::function<double(SPoint)> exact_function =
[](SPoint x)->double { return std::sin(2.0*M_PI*x.x())*std::sin(2.0*M_PI*x.y()); };

// rhs function in Poisson equation
const static std::function<double(SPoint)> rhs_function =
[](SPoint x)->double { return  8.0*M_PI*M_PI*std::sin(2.0*M_PI*x.x())*std::sin(2.0*M_PI*x.y()); };

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
  public ScalarAnalyticalExpression
{

  double evaluate(const hydrofem::SPoint &x) override
  { assert(!m_is_transient); return rhs_function(x); }

  double operator()(const hydrofem::SPoint &x) override
  { return this->evaluate(x); }

};
  
void Problem_Poisson::init()
{
  if (m_is_initialized)
    return;

  auto exact_ = std::make_shared<Exact_Function>(); assert(exact_);
  this->setExactSolution(exact_);

  m_rhs_fnc = std::make_shared<RHS_Function>(); assert(m_rhs_fnc);

  // initialize once mesh is created or set
  if (m_bc->mesh())
    std::dynamic_pointer_cast<BC_Scalar>(m_bc)->initializeBoundaryPointsToDirichletEverywhere();

  m_dirichlet_bc_fnc = exact_;
  m_is_initialized = true;

}

}
// end namespace hydrofem
