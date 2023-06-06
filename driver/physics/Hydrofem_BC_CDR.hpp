// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER
 
#ifndef __Hydrofem_BC_CDR_HPP__
#define __Hydrofem_BC_CDR_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_BC_Scalar.hpp"

namespace hydrofem
{

/**
 * \brief A class for general boundary 
 *        conditions for scalar problems
 */
class BC_CDR
  :
  public BC_Scalar
{
public:

  using Ptr = std::shared_ptr<BC_CDR>;
  
  BC_CDR()
  {
    m_name = "bc-cdr";
  }

  //! \brief Ctor from \p mesh 
  BC_CDR(const std::shared_ptr<Mesh>& mesh)
    :
    BC_Scalar(mesh)
  {
    m_name = "bc-cdr";
  }

  //! \brief Ctor from \p mesh and \p velocity 
  BC_CDR(const std::shared_ptr<Mesh>& mesh,
                   const std::function<SPoint(SPoint)>& velocity)
    :
    BC_Scalar(mesh)
  {
    m_name = "bc-cdr";
    m_velocity = velocity;
  }
  
  //! \brief sets the fluid velocity profile
  void setFluidVelocity(const std::function<SPoint(SPoint)>& velocity)
  { m_velocity = velocity; }

  //! \brief gets the fluid velocity profile
  const std::function<SPoint(SPoint)>& getFluidVelocity() const
  { return m_velocity; }

  //! \brief Dtor
  virtual ~BC_CDR() {}
  
  void initialize() override;

private:
  
  std::function<SPoint(SPoint)> m_velocity;  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BC_CDR_HPP__ */
