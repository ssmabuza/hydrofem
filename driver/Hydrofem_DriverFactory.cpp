// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_DriverFactory.hpp"

#include "Hydrofem_Driver_SteadySingleEquationSet.hpp"
#include "Hydrofem_Driver_TransientSingleEquationSet.hpp"

namespace hydrofem
{

void DriverFactory::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("driver",po::value<std::string>(&m_driver_name)->default_value("steady-single-equationset"),"The driver name");
}

std::shared_ptr<Driver> DriverFactory::
buildDriver()
{
  std::shared_ptr<Driver> driver;
  if (m_driver_name == "steady-single-equationset")
    driver = std::make_shared<Driver_SteadySingleEquationSet>(m_option_handler);
  else if (m_driver_name == "transient-single-equationset")
    driver = std::make_shared<Driver_TransientSingleEquationSet>(m_option_handler);
  assert(driver);
  return driver;
}

}
// end namespace hydrofem
