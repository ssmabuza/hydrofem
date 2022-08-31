// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_MeshFactory_HPP__
#define __Hydrofem_MeshFactory_HPP__

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

class MeshFactory
  :
  public Optionable
{
public:
  
  /** @brief Ctor */
  MeshFactory(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Optionable(option_handler)
  {
    m_option_handler = option_handler;
    option_handler->parse();
  }
  
  /** @brief Dtor */
  virtual ~MeshFactory() {}
  
  /** @brief The mesh builder routine */
  std::shared_ptr<Mesh> buildMesh() const;

  /** \brief Flag for writing to MATLAB */
  bool writeOutputToMATLAB() const 
  { return m_write_to_matlab; }
  
private:

  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;
  
  /** get the valid parameters */
  void addOptionsCallback(po::options_description &config);
  
  
  int m_ndims;
  std::string m_input;
  std::string m_elem_type;
  std::string m_triangulation;
  bool m_write_to_matlab;
  double m_x0;
  double m_xf;
  double m_y0;
  double m_yf;
  double m_z0;
  double m_zf;
  int m_nx;
  int m_ny;
  int m_nz;
  std::string m_filename;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_MeshFactory_HPP__ */
