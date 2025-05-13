// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Assembler_Base_HPP__
#define __Hydrofem_Assembler_Base_HPP__

#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_LinearObjectBuilder.hpp"

namespace hydrofem
{

class Assembler_Base
{
public:

  using LOB = LinearObjectBuilder;
  
  //! \brief Ctor
  Assembler_Base() = default;

  //! \brief Dtor
  virtual ~Assembler_Base() = default;

  /**
   * \brief The transient assembly explicit and implicit stepping including explicit RK steppers
   *
   * \param U       - current solution
   * \param U_dot   - current time derivative of solution 
   * \param res_U   - residual
   * \param jac_U   - Jacobian
   * \param time    - current time
   * \param delta_t - current time step
   * \param beta    - system multiplier (which is 1 most of the time)
   */
  virtual
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& /*U*/,
                                const std::shared_ptr<const FEVector>& /*U_dot*/,
                                const std::shared_ptr<FEVector>& /*res_U*/,
                                const std::shared_ptr<FEMatrix>& /*jac_U*/,
                                const double /*time*/,
                                const double /*delta_t*/,
                                const double /*beta*/) const { }

  /**
   * \brief The steady part of the residual in
   *
   * \param U       - current solution
   * \param U_dot   - current time derivative of solution (finite difference approx)
   * \param res_U   - residual
   * \param jac_U   - Jacobian
   * \param time    - current time
   * \param delta_t - current time step
   * \param beta    - system multiplier (which is 1 most of the time)
   */
  virtual
  void buildSteadyResidualAndJacobian(const std::shared_ptr<const FEVector>& /*U*/,
                                      const std::shared_ptr<const FEVector>& /*U_dot*/,
                                      const std::shared_ptr<FEVector>& /*res_U*/,
                                      const std::shared_ptr<FEMatrix>& /*jac_U*/,
                                      const double /*time*/,
                                      const double /*delta_t*/,
                                      const double /*beta*/) const { }

  /**
   * \brief The non-transient assembly (steady solver)
   *
   * \param U     - current solution
   * \param res_U - residual
   * \param jac_U - Jacobian
   * \param beta  - system multiplier
   */
  virtual
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& /*U*/,
                                const std::shared_ptr<FEVector>& /*res_U*/,
                                const std::shared_ptr<FEMatrix>& /*jac_U*/,
                                const double /*beta*/) const { }

  /**
   * \brief The non-transient assembly (steady solver)
   *
   * \param U     - current solution
   * \param res_U - residual
   * \param jac_U - Jacobian
   * \param beta  - system multiplier
   */
  virtual
  void buildSteadyResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                      const std::shared_ptr<FEVector>& res_U,
                                      const std::shared_ptr<FEMatrix>& jac_U,
                                      const double beta) const 
  {
    buildResidualAndJacobian(U,res_U,jac_U,beta);
  }

  /**
   * \brief Applies Dirichlet BC to global system
   * 
   * 
   */
  virtual
  void applyDirichletBC(const std::shared_ptr<FEVector>& /*res_U*/,
                        const std::shared_ptr<FEMatrix>& /*jac_U*/) const { }
                                
  /**
   * \brief this will properly call fill complete in the Jacobian and do any post-processing
   *
   * \param res_U
   * \param jac_U
   */
  virtual void finalizeAssembly(const std::shared_ptr<FEVector>& res_U,
                                const std::shared_ptr<FEMatrix>& jac_U) const
  {
    applyDirichletBC(res_U,jac_U);
  }
                                
  /** \brief */                                
  virtual std::shared_ptr<LOB> lob() { return m_lob; }

  /** \brief Get the matrix from dR/dU_dot. This method is required by some classic steppers and some explicit steppers. */
  std::shared_ptr<FEMatrix> getTimeDerivativeMatrix() const 
  { return nullptr; }

protected:
  
  // user built LOB
  std::shared_ptr<LOB> m_lob;
  // apply Dirichlet only once
  bool m_apply_dirichlet_only_once = false;
  // dirichlet has been applied
  bool m_dirichlet_applied = false;
  // flag to check if we have a steady problem or not

  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Assembler_Base_HPP__ */
