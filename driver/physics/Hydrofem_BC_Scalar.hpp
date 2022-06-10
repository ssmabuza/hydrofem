// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_BC_Scalar_HPP__
#define __Hydrofem_BC_Scalar_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_BC.hpp"
#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_Equation.hpp"

namespace hydrofem
{

/**
 * \brief A class for general boundary 
 *        conditions for scalar problems
 *        on some nodal scalar finite element.
 */
class BC_Scalar
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
    
    // edges
    typeGammaZero,
    typeGammaPlus,
    typeGammaMinus,
    typeSigma,
    typeInterface
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
  
  //! \brief Ctor
  explicit BC_Scalar() : BC() {}

  //! \brief Ctor from \p mesh
  explicit BC_Scalar(const std::shared_ptr<Mesh>& mesh) : BC(mesh)
  {
    m_name = "Scalar BC";
  }
  
  //! \brief Dtor 
  virtual ~BC_Scalar() {}
  
  // special case: set Dirichlet everywhere
  void initializeBoundaryPointsToDirichletEverywhere();
  
  // iterator over points/edges 
  using BCIt = std::unordered_map<int,std::unique_ptr<BCInfo>>::iterator;
  // constant iterator over points/edges
  using BCConstIt = std::unordered_map<int,std::unique_ptr<BCInfo>>::const_iterator;
  // container for BC data
  using BCIterable = std::unordered_map<int,std::unique_ptr<BCInfo>>;
  
  
  BCIterable bcPoints() const
  { return m_bc_info_points; }

  BCIterable bcEdges() const
  { return m_bc_info_edges; }
  
  BCConstIt startBCPoints() const 
  {
    return m_bc_info_points.begin();
  }

  BCConstIt endBCPoints() const 
  {
    return m_bc_info_points.end();
  }
  
  BCConstIt startBCEdges() const 
  {
    return m_bc_info_edges.begin();
  }

  BCConstIt endBCEdges() const 
  {
    return m_bc_info_edges.end();
  }
  
  bool isDirichlet(BCConstIt it) const
  { return (it->second->m_boundary_condition_type == typeDirichlet); }
  
  bool isGammaPlus(BCConstIt it) const
  { return (it->second->m_boundary_condition_type == typeGammaPlus); }

  bool isGammaMinus(BCConstIt it) const
  { return (it->second->m_boundary_condition_type == typeGammaMinus); }
  
  bool isGammaZero(BCConstIt it) const
  { return (it->second->m_boundary_condition_type == typeGammaZero); }
  
  bool isSigma(BCConstIt it) const
  { return (it->second->m_boundary_condition_type == typeSigma); }
  
protected:
  
  void initialize() override;

  // this is for weak/integral based boundary calculations in 2D
  std::unordered_map<int,std::unique_ptr<BCInfo>> m_bc_info_edges;
  // this is for point/nodal based boundary conditions
  std::unordered_map<int,std::unique_ptr<BCInfo>> m_bc_info_points;
  // this is for point/nodal based boundary conditions periodic
  std::unordered_map<int,std::unique_ptr<BCInfoPeriodic>> m_bc_info_points_per;
  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BC_Scalar_HPP__ */
