// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Assembler_Bioseparation_MatrixForm_HPP__
#define __Hydrofem_Assembler_Bioseparation_MatrixForm_HPP__


#include "Hydrofem_Assembler_Bioseparation.hpp"

namespace hydrofem
{



/**
 * \brief A matrix form of the theta scheme for the bioseparation problem with and without AFC stabilization
 *
 */
class Assembler_Bioseparation_MatrixForm
  :
  public Assembler_Bioseparation
{
public:

  using LOB = Assembler_Bioseparation::LOB;

  /** \brief Ctor */
  Assembler_Bioseparation_MatrixForm(const std::shared_ptr<Problem>& problem,
                          const std::shared_ptr<DofMapper>& dofmapper,
                          const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& basis,
                          const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature, bool afc_enabled)
    :
    Assembler_Bioseparation(problem,dofmapper,basis,quadrature,afc_enabled)
  {

  }

  /** \brief Dtor */
  ~Assembler_Bioseparation_MatrixForm() override = default;

  /** \brief  */
  void buildResidualAndJacobian(const std::shared_ptr<const FEVector>& U,
                                const std::shared_ptr<const FEVector>& U_dot,
                                const std::shared_ptr<FEVector>& res_U,
                                const std::shared_ptr<FEMatrix>& jac_U,
                                const double time,
                                const double delta_t,
                                const double beta) const override;


private:

  //! update m_M_mat
  void updateMassMatrix();
  //! compute m_S_mat once
  void buildStiffnessMatrix();
  //! compute m_K_mat once
  void buildConvectionMatrix();

  //
  bool m_stiffness_matrix_computed = false;
  bool m_convection_matrix_computed = false;

  // global stiffness matrix
  std::shared_ptr<FEMatrix> m_S_mat;
  // global convection matrix
  std::shared_ptr<FEMatrix> m_K_mat;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Assembler_Bioseparation_MatrixForm_HPP__ */

