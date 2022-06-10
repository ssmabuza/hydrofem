// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOBase_HPP__
#define __Hydrofem_IOBase_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_FunctionElement.hpp"

namespace hydrofem
{
  
/** 
 * \brief A base class for writing a scalar field info to file
 */
class IOBase
{
public:

  /** \brief Ctor */
  IOBase(const std::shared_ptr<DofMapper>& dofmapper,
         const std::shared_ptr<FunctionElement<double>>& fe_shape)
  {
    m_dofmapper = dofmapper;
    m_fe_shape = fe_shape;
    m_field_name = "";
    m_directory_name = "";
    m_time_final = - std::numeric_limits<double>::max();
  }
        
  /** \brief Dtor */
  virtual ~IOBase() = default;

  /** \brief Main routine for writing high order solution into file : cell basis format */
  virtual void writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const {}
  
  /** \brief Main routine for writing high order solution into file : distr vector format */
  virtual void writeSolutionTofile(const std::shared_ptr<FEVector>& sol, const double time) const {}

  void setFieldName(const std::string& field_name)
  { m_field_name = field_name; }
  
  void setDIRName(const std::string& dir_name)
  { m_directory_name = dir_name; }

  void setFinalTime(const double time)
  { m_time_final = time; }

protected:
  
  std::string getPVUFilename() const
  { return m_directory_name + "/" + m_field_name + ".pvd"; }
  
  void createPVUOutputFile() const;
  void updatePVUOutputFile(const std::string file_name, const double time) const;
  void finalizePVUOutputFile() const;

  // Dof manager with mesh and spatial information
  std::shared_ptr<DofMapper> m_dofmapper;
  
  // Basis functions for basis point evaluations (in subcell formation)
  mutable std::shared_ptr<FunctionElement<double>> m_fe_shape;

  // name of the field being saved
  std::string m_field_name;
  
  // path of directory where field is saved
  std::string m_directory_name;

  // counter for output
  mutable int m_counter = 0;

  // final time
  double m_time_final;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOBase_HPP__ */
