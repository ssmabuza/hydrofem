// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_BC_EulerEquations_HPP__
#define __Hydrofem_BC_EulerEquations_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_BC.hpp"
#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Equation_Euler.hpp"

namespace hydrofem
{

/**
 * \brief The main implementation of boundary conditions for Euler equations
 *
 * The user has to specify the mesh on which the problem is solved.
 * The user also has to state the solution values to determine the
 * boundary condition types at each time step in the solver.
 */
class BC_EulerEquations
  :
  public BC
{
protected:

  /**
   * \brief Boundary condition labels
   */
  enum BCType
  {
    // nodes + edges
    typeGenericBoundary,
    typeInternal,
    
    // nodes
    typeNeumann,
    typePeriodic,
    typeDirichlet,
    typeWeakDirichlet,
    
    // edges (full system)
    typeSupersonicInlet,
    typeSupersonicOutlet,
    typeSubsonicInlet,
    typeSubsonicOutlet,
    typeSolidWall
  };  
  
  struct BCInfo {
    /// Boundary condition type
    BCType m_boundary_condition_type;
    /// Index on this type of BC
    int m_index_in_system;
  };
  
  struct BCInfoPeriodic {
    // BC information
    BCInfo m_bc;
    // corresponding index
    int m_index_match;
  };
  
public:

  //! Ctor
  BC_EulerEquations(const std::shared_ptr<Mesh>& mesh) : BC(mesh)
  {
    m_name = "bc-euler-equations";
  }
  
  //! Ctor
  BC_EulerEquations(const std::shared_ptr<Mesh>& mesh,
                    const std::shared_ptr<Equation_Euler<double>>& equation) : BC(mesh)
  {
    m_name = "bc-euler-equations";
    m_equation = equation;
  }
  
  //! Dtor
  virtual ~BC_EulerEquations() {}
  
protected:
  
  void initialize() override;
  
  std::shared_ptr<Equation_Euler<double>> m_equation;

  // this is for bc calculations over boundary edges   
  std::unordered_map<int,std::unique_ptr<BCInfo>> m_bc_info_edges;
  // this is for point/nodal based boundary conditions
  std::unordered_map<int,std::unique_ptr<BCInfo>> m_bc_info_points;
  // this is for point/nodal based boundary conditions periodic
  std::unordered_map<int,std::unique_ptr<BCInfoPeriodic>> m_bc_info_points_per;
  
  //! Setup style for BC
  std::string m_bcForm;

  // flag for fixed bcs 
  bool m_fixed_bc_set;
  
};
// end class BC_EulerEquations

}
// end namespace hydrofem

#endif /** __Hydrofem_BC_EulerEquations_HPP__ */
