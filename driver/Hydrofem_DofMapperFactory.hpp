// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_DofMapperFactory_HPP__
#define __Hydrofem_DofMapperFactory_HPP__

#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

class DofMapperFactory
  :
  public Optionable
{
public:
  
  /** @brief Ctor */
  DofMapperFactory() = default;
  
  /**
   * @brief Ctor from parameter list
   * @param pl - input parameters
   * */
  explicit DofMapperFactory(const std::shared_ptr<OptionHandler>& option_handler) : Optionable(option_handler)
  {
    m_option_handler = option_handler;
    m_option_handler->parse();
  }

  /** @brief Dtor */
  virtual ~DofMapperFactory() = default;

  /** @brief builds dofmapper using the existing parameter list */
  [[nodiscard]] std::shared_ptr<DofMapper> buildDofMapper(const std::shared_ptr<Mesh>& mesh) const;
  
private:

  /** @brief builds dofmapper from mesh and order of basis */
  static std::shared_ptr<DofMapper> buildDofMapper(const std::shared_ptr<Mesh>& mesh,
                                                   const int basis_order);
  
  /** @brief options to be parsed for solver */
  void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("basisOrder",po::value<int>(&m_basis_order)->default_value(1),"The finite element basis order");
  }
  
  // basis order
  int m_basis_order;
  
  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_DofMapperFactory_HPP__ */
