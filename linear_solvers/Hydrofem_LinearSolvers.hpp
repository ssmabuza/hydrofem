// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LinearSolvers_Tpetra_HPP__
#define __Hydrofem_LinearSolvers_Tpetra_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

class LinearSolvers
{
  
  static const int maxiter = 10000;
  constexpr static double tol = 1.0e-06;
  constexpr static double eps = 1.0e-15;

public:

  //! \brief method solves a linear system using the Jacobi preconditioned CG
  static void solveSystemCG(const std::shared_ptr<const FEMatrix>& A,
                            const std::shared_ptr<const FEVector>& rhs,
                            const std::shared_ptr<FEVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned CG
  static void solveSystemBiCGSTAB(const std::shared_ptr<const FEMatrix>& A,
                                  const std::shared_ptr<const FEVector>& rhs,
                                  const std::shared_ptr<FEVector>& sol);

  //! \brief method solves a linear system using the Jacobi preconditioned GMRES
  static void solveSystemGMRES(const std::shared_ptr<const FEMatrix>& A,
                               const std::shared_ptr<const FEVector>& rhs,
                               const std::shared_ptr<FEVector>& sol);
  
};


class LinearSolverInterface
{
public:
  
  LinearSolverInterface(std::string method_name) : m_name(method_name) {}
  
  inline
  void solve(const std::shared_ptr<const FEMatrix>& A,
             const std::shared_ptr<const FEVector>& rhs,
             const std::shared_ptr<FEVector>& sol) const
  {
    if (m_name=="cg")
      LinearSolvers::solveSystemCG(A,rhs,sol);
    else if (m_name=="bicgstab")
      LinearSolvers::solveSystemBiCGSTAB(A,rhs,sol);
    else if (m_name=="gmres")
      LinearSolvers::solveSystemGMRES(A,rhs,sol);
    else
    {
      std::stringstream ss;
      ss << "Error in LinearSolverInterface!!! method name " << m_name << " invalid, aborting!!";
      throw std::runtime_error(ss.str());
    }
  }
  
private:

  std::string m_name;  
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_LinearSolvers_Tpetra_HPP__ */
