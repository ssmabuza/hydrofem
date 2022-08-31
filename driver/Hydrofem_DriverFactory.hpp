// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_DriverFactory_HPP__
#define __Hydrofem_DriverFactory_HPP__

#include "Hydrofem_Driver.hpp"

namespace hydrofem
{

/**
 * \brief A factory class for creating the proper driver
 */
class DriverFactory
  :
  public Optionable
{
public:
  
  DriverFactory(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Optionable(option_handler),
    m_option_handler(option_handler)
  {
    option_handler->parse();
  }
  
  virtual ~DriverFactory() {}
  
  std::shared_ptr<Driver> buildDriver();

private:

  /** \brief options to be parsed for solver */
  void addOptionsCallback(po::options_description &config);
  
  // name of the driver
  std::string m_driver_name;
  // option handler 
  std::shared_ptr<OptionHandler> m_option_handler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_DriverFactory_HPP__ */
