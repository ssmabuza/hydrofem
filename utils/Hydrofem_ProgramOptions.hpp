// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ProgramOptions_HPP__
#define __Hydrofem_ProgramOptions_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

/**
 * \class ProgramOptions A class that provides a simple inderface for boost programoptions
 * 
 */
class ProgramOptions
{
public:
  
  //! @brief Ctor
  ProgramOptions(int argc, char **argv)
    :
    m_desc(),
    m_vm(),
    m_argc(argc),
    m_argv(argv) 
  {}
  
  //! @brief Copy Ctor
  explicit ProgramOptions(const ProgramOptions& otherOpts)
    :
    m_desc(otherOpts.m_desc),
    m_argc(otherOpts.m_argc), 
    m_argv(otherOpts.m_argv)
  {}

  //! @brief Dtor
  ~ProgramOptions() = default;
  
  //! @brief parse
  inline void verifyAllOptions()
  {
    boost::program_options::parsed_options parsed 
      = boost::program_options::command_line_parser(m_argc, m_argv).options(m_desc).allow_unregistered().run();
    boost::program_options::store(parsed,m_vm);
    boost::program_options::notify(m_vm);    
  }
  
  //! get the options description
  inline boost::program_options::options_description& desc()
  { return m_desc; }
  
private:
  
  boost::program_options::options_description m_desc;
  boost::program_options::variables_map m_vm;
  int m_argc;
  char** m_argv;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_ProgramOptions_HPP__ */
