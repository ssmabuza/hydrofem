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
  
static
double boundaryFunction(const double /*x*/, const double /*y*/) { return 1.0; }

class BoundaryExpr : public ScalarAnalyticalExpression
{
public:

  BoundaryExpr() { m_is_transient = true; }

  ~BoundaryExpr() = default;

  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& x)
  { assert(!m_is_transient); return boundaryFunction(x.x(),x.y()); }

  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& x, const double /*t*/)
  { assert(m_is_transient); return boundaryFunction(x.x(),x.y()); }

  //! \brief make this a functor
  virtual double operator()(const SPoint& x)
  { assert(!m_is_transient); return boundaryFunction(x.x(),x.y()); }

  //! \brief make this a functor
  virtual double operator()(const SPoint& x, const double /*t*/)
  { assert(m_is_transient); return boundaryFunction(x.x(),x.y()); }

};

void Problem_Bioseparation::init()
{
  if (m_is_initialized) return;
  m_ic = std::make_shared<ConstantScalarInitialCondition>(0.0);
  m_u_in = std::make_shared<BoundaryExpr>();
  auto _bc = std::make_shared<BC_Bioseparation>(); // TODO: set the mesh
  const double pi = M_PI;
  const double fr = m_flowrate;
  const double width = m_xf - m_x0;
  
  std::function<SPoint(SPoint)> vel = [&pi,&fr,&width](SPoint x)->SPoint
  { return SPoint(0.0,3*fr/(4*pi*std::pow(width/2.0,3))*(width-x.x())*x.x()); };
  
  _bc->setFluidVelocity(std::make_shared<std::function<SPoint(SPoint)>>(vel));
  m_bc = _bc;
  m_is_initialized = true;
}

void Problem_Bioseparation::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("prob-biosep-x0",po::value<double>(&m_x0)->default_value(0.0),"Left value for x.")
    ("prob-biosep-xf",po::value<double>(&m_xf)->default_value(1.6),"Right value for x.")
    ("prob-biosep-omega",po::value<double>(&m_omega)->default_value(0.84),"Porosity value")
    ("prob-biosep-rho_s",po::value<double>(&m_rho_s)->default_value(1.0),"Density of membrane")
    ("prob-biosep-q_max",po::value<double>(&m_q_max)->default_value(150.0),"Maximum binding capacity")
    ("prob-biosep-flowrate",po::value<double>(&m_flowrate)->default_value(0.1),"Maximum flowrate of the carrying fluid")
    ("prob-biosep-K_eq",po::value<double>(&m_K_eq)->default_value(2.06),"Langmuir equilibrium constant")
    ("prob-biosep-alphaL",po::value<double>(&m_alphaL)->default_value(1.37),"Langmuir equilibrium constant")
    ("prob-biosep-alphaT",po::value<double>(&m_alphaT)->default_value(0.137),"Langmuir equilibrium constant")
    ("prob-biosep-d0",po::value<double>(&m_d0)->default_value(0.0000228),"Langmuir equilibrium constant");
}

}
// end namespace hydrofem
