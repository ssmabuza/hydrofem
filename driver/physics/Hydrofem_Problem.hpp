// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Problem_HPP__
#define __Hydrofem_Problem_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_BC.hpp"

namespace hydrofem
{

/**
 * \brief The base class that describes a continuous problem
 */
class Problem
{
public:
  
  Problem() = default;
  
  virtual ~Problem() = default;
  
  virtual void init() = 0;
  
  //! \brief get functions
  //@{
  //! \brief get the name of this problem
  std::string name() const
  { return m_name; }
  
  std::vector<std::string>
  dofNames() const
  { return m_dof_names; }
  
  //@}
  
  //! \brief All the main get functions for building up the solver
  //@{
  void setName(const std::string name)
  { m_name = name; }
  
  void setDofNames(const std::vector<std::string>& dof_names)
  { m_dof_names = dof_names; }
  //@}
  
  [[nodiscard]] virtual std::shared_ptr<BC> bc() const
  { return nullptr; }
  
protected:
  
  // name
  std::string m_name;
  // field names
  std::vector<std::string> m_dof_names;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Problem_HPP__ */
