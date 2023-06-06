// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LinearSolvers.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

namespace hydrofem
{

void
LinearSolvers::
solveSystemCG(const std::shared_ptr<const FEMatrix>& A,
              const std::shared_ptr<const FEVector>& rhs,
              const std::shared_ptr<FEVector>& sol)
{
  // create the CG solver object
  Eigen::ConjugateGradient<FEMatrix, Eigen::Lower|Eigen::Upper> solver;
  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);
  solver.compute(*A);
  *sol = solver.solve(*rhs);
  if (solver.info() == Eigen::ComputationInfo::Success)
  {
    std::cout << "PCG solver converged\n" 
              << "Number of iterations = " << solver.iterations() << '\n'
              << "Error at convergence = " << solver.error() << "." 
              << std::endl;
  } else {
    std::stringstream ss;
    ss << "PCG solver did not converge\n"
       << "Solver failed after " << solver.iterations() << " iterations." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void
LinearSolvers::
solveSystemBiCGSTAB(const std::shared_ptr<const FEMatrix>& A,
                    const std::shared_ptr<const FEVector>& rhs,
                    const std::shared_ptr<FEVector>& sol)
{
  // create the BiCGSTAB solver object
  Eigen::BiCGSTAB<FEMatrix , Eigen::IncompleteLUT<double>> solver;
  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);
  solver.compute(*A);
  *sol = solver.solve(*rhs);
  if (solver.info() == Eigen::ComputationInfo::Success)
  {
    std::cout << "BiCGSTAB solver converged\n" 
              << "Number of iterations = " << solver.iterations() << '\n'
              << "Error at convergence = " << solver.error() << "." 
              << std::endl;
  } else {
    std::stringstream ss;
    ss << "BiCGSTAB solver did not converge\n"
       << "Solver failed after " << solver.iterations() << " iterations." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

void
LinearSolvers::
solveSystemGMRES(const std::shared_ptr<const FEMatrix>& A,
                 const std::shared_ptr<const FEVector>& rhs,
                 const std::shared_ptr<FEVector>& sol)
{
  // build GMRes Solver
  Eigen::GMRES<FEMatrix, Eigen::IncompleteLUT<double>> solver;
  solver.setTolerance(tol);
  solver.setMaxIterations(maxiter);
  solver.compute(*A);
  *sol = solver.solve(*rhs);
  if (solver.info() == Eigen::ComputationInfo::Success)
  {
    std::cout << "GMRES solver converged\n" 
              << "Number of iterations = " << solver.iterations() << '\n'
              << "Error at convergence = " << solver.error() << "." 
              << std::endl;
  } else {
    std::stringstream ss;
    ss << "GMRES solver did not converge\n"
       << "Solver failed after " << solver.iterations() << " iterations." << std::endl;
    throw std::runtime_error(ss.str());
  }
}

}
// end namespace hydrofem
