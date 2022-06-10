// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER
 
#ifndef __Hydrofem_BC_Bioseparation_HPP__
#define __Hydrofem_BC_Bioseparation_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_BC_Scalar.hpp"
#include "Hydrofem_Equation.hpp"

namespace hydrofem
{

/**
 * \brief A class for general boundary 
 *        conditions for scalar problems
 */
class BC_Bioseparation
  :
  public BC_Scalar
{
public:
  
  //! \brief Ctor from \p mesh 
  BC_Bioseparation(const std::shared_ptr<Mesh>& mesh)
    :
    BC_Scalar(mesh)
  {
    m_name = "bc-bioseparation";
  }

  //! \brief Ctor from \p mesh and \p velocity 
  BC_Bioseparation(const std::shared_ptr<Mesh>& mesh,
                   const std::shared_ptr<std::function<SPoint(SPoint)>>& velocity)
    :
    BC_Scalar(mesh)
  {
    m_name = "bc-bioseparation";
    m_velocity = velocity;
  }
  
  //! \brief sets the fluid velocity profile
  void setFluidVelocity(const std::shared_ptr<std::function<SPoint(SPoint)>>& velocity)
  { m_velocity = velocity; }
  
  //! \brief Dtor
  virtual ~BC_Bioseparation() {}
  
private:
  
  void initialize();

  std::shared_ptr<std::function<SPoint(SPoint)>> m_velocity;  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BC_Bioseparation_HPP__ */
