// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_AssemblerFactory_HPP__
#define __Hydrofem_AssemblerFactory_HPP__

#include "Hydrofem_OptionHandler.hpp"
#include "Hydrofem_Assembler_Base.hpp"

namespace hydrofem
{

class AssemblerFactory
  :
  public Optionable
{
public:

  //! \brief Ctor
  AssemblerFactory(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Optionable(option_handler)
  {
    m_option_handler = option_handler;
  }
  
  //! \brief Dtor
  ~AssemblerFactory() {}
  
  //! \brief build function
  std::shared_ptr<Assembler_Base> 
  build() const;
  
private:
  
  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("Assembler",po::value<std::string>(&m_name)->default_value("Poisson"),"Assembler name.");
  }
  
  // name of the assembler 
  std::string m_name;
  // option handler to pass to assembler
  std::shared_ptr<OptionHandler> m_option_handler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_AssemblerFactory_HPP__ */
