// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_BC_HPP__
#define __Hydrofem_BC_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_Equation.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"

namespace hydrofem
{

class Assembler_Base;  
  
/**
 * \brief A base class for the boundary conditions
 *        BC info is described on boundary edges
 */
class BC 
{
public:
  
  //! \brief Ctor
  BC() {}
  
  //! \brief Ctor
  explicit BC(const std::shared_ptr<Mesh>& mesh)
  {
    m_mesh = mesh;
  }
  
  //! \brief Dtor
  virtual ~BC() {}
  
  //! \brief get the name
  inline std::string name() { return m_name; }
  
  //! \brief get the underlying mesh describing the domain
  inline std::shared_ptr<Mesh> mesh() const 
  { return m_mesh; }
  
  //! \brief set the mesh to be used
  inline void setMesh(const std::shared_ptr<Mesh>& mesh)
  { m_mesh = mesh; }
  
  //! \brief main initialization routine to set BC info
  virtual void initialize() = 0;
  
  // class where boundary calculations happen
  friend class Assembler_Base;
  
protected:
  
  //! @brief the BC name, default is "BC", this must match the identifier in the input file
  std::string m_name = "BC Base";
  
  //! @brief the global mesh
  std::shared_ptr<Mesh> m_mesh;
  
};

class DirichletBCFunction
{
public:
  
  //!
  DirichletBCFunction(const std::shared_ptr<ScalarAnalyticalExpression>& bc_fnc)
  {
    m_bc_fnc = bc_fnc;
  }
  
  double operator()(SPoint x) const
  { return m_bc_fnc(x); }
  
protected:

  // boundary function, can be evaluated at x & t
  std::shared_ptr<ScalarAnalyticalExpression> m_bc_fnc;
  
};

class VectorDirichletBCFunction
{
public:
  
  //! \brief Ctor 
  VectorDirichletBCFunction(const std::shared_ptr<AnalyticalExpression>& bc_fnc)
  {
    m_bc_fnc = bc_fnc;
  }

  LVec operator() (SPoint x) const
  { return m_bc_fnc(x); }

protected:

  // boundary function, can be evaluated at x & t
  std::shared_ptr<AnalyticalExpression> m_bc_fnc;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BC_HPP__ */

