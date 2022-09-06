// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_AssemblerFactory_HPP__
#define __Hydrofem_AssemblerFactory_HPP__

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Problem.hpp"
#include "Hydrofem_Quadrature.hpp"
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
  AssemblerFactory(const std::shared_ptr<OptionHandler>& option_handler,
                   const std::shared_ptr<Problem>& problem,
                   const std::shared_ptr<DofMapper>& dofmapper,
                   const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                   const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature)
    :
    Optionable(option_handler)
  {
    option_handler->parse();
    m_dofmapper = dofmapper;
    m_basis = basis;
    m_quadrature = quadrature;
    m_problem = problem;
  }
  
  //! \brief Dtor
  ~AssemblerFactory() {}
  
  //! \brief build function
  std::shared_ptr<Assembler_Base> 
  build() const;
  
private:
  
  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config) override
  {
    config.add_options()
      ("afc-enabled",po::value<bool>(&m_afc_enabled)->default_value(false),"Flag for AFC, goes into assembly.")
      ("assembler",po::value<std::string>(&m_name)->default_value("poisson"),"Assembler name.");
  }
  
  // name of the assembler 
  std::string m_name;
  // problem defining the Poisson equation
  std::shared_ptr<Problem> m_problem;
  // the dof manager
  std::shared_ptr<DofMapper> m_dofmapper;
  // basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // global quadrature
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  // flag for doing AFC
  bool m_afc_enabled;
};

}
// end namespace hydrofem

#endif /** __Hydrofem_AssemblerFactory_HPP__ */
