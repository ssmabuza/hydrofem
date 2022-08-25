// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_BCFactory_HPP__
#define __Hydrofem_BCFactory_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem 
{

class BC;
class Mesh;
template <typename ScalarT>
class Equation;

class BCFactory
  :
  public Optionable
{
public:
  
  BCFactory(const std::shared_ptr<OptionHandler>& option_handler,
            const std::shared_ptr<Mesh>& mesh)
    :
    Optionable(option_handler)
  {
    m_option_handler = option_handler;
    m_mesh = mesh;
    m_option_handler->parse();
  }
  
  virtual ~BCFactory() {}
  
  std::shared_ptr<BC> build();

  void addOptionsCallback(po::options_description& config)
  {
    config.add_options()
      ("bc",po::value<std::string>(&m_bc_name)->default_value("bc-scalar"),"The BC name");
  }
  
private:

  std::shared_ptr<OptionHandler> m_option_handler;
  std::shared_ptr<Mesh> m_mesh;
  std::string m_bc_name;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BCFactory_HPP__ */
