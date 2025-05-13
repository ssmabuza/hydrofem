// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_CDR.hpp"

#include "Hydrofem_BC_CDR.hpp"
#include "Hydrofem_InitialCondition.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{

namespace {

static
double boundaryFunction(const double /*x*/, const double /*y*/)
{
  return 0.0;
}

auto initcondFunction = [](SPoint p)->decltype(auto)
{
  return boundaryFunction(p.x(),p.y());
};


class BoundaryExpr : public ScalarAnalyticalExpression
{
public:

  BoundaryExpr() { m_is_transient = true; }

  ~BoundaryExpr() = default;

  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& x)
  { return boundaryFunction(x.x(),x.y()); }

  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& x, const double /*t*/)
  { return boundaryFunction(x.x(),x.y()); }

  //! \brief make this a functor
  virtual double operator()(const SPoint& x)
  { return boundaryFunction(x.x(),x.y()); }

  //! \brief make this a functor
  virtual double operator()(const SPoint& x, const double /*t*/)
  { return boundaryFunction(x.x(),x.y()); }

};

}

void Problem_CDR::init()
{
  if (m_is_initialized) return;
  m_ic = std::make_shared<FunctionalScalarInitialCondition>(initcondFunction);
  m_u_in = std::make_shared<BoundaryExpr>();
  auto _bc = std::make_shared<BC_CDR>(); // TODO: set the mesh
  const double H = m_H;
  const double vMax = m_v_max;
  
  std::function<SPoint(SPoint)> vel = [=](SPoint x)->SPoint
  { return SPoint(x.y()/H * (1.0 - x.y()/H) * 4.0 * vMax,0.0); };
  
  _bc->setFluidVelocity(vel);
  m_bc = _bc;
  m_is_initialized = true;
}

void Problem_CDR::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("prob-cdr-L",po::value<double>(&m_L)->default_value(0.0),"Characteristic length of 2D channel.")
    ("prob-cdr-H",po::value<double>(&m_H)->default_value(1.6),"Characteristic height of 2D channel")
    ("prob-cdr-diff",po::value<double>(&m_diff)->default_value(0.84),"Diffusion coefficient")
    ("prob-cdr-v_max",po::value<double>(&m_v_max)->default_value(0.0000228),"Maximum velocity in parabolic profile");
}

}
// end namespace hydrofem
