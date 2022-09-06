// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Driver_TransientSingleEquationSet_HPP__
#define __Hydrofem_Driver_TransientSingleEquationSet_HPP__

#include "Hydrofem_Driver.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_GlobalGather.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_OptionHandler.hpp"

namespace hydrofem
{

class Mesh;
class IOBase;
class FEBasis; 
class Stepper;
class Problem;
class DofMapper;
class Quadrature;
class NewtonSolver;
class Assembler_Base;

/**
 * \brief The transient driver for solving time dependent problems.
 *
 */
class Driver_TransientSingleEquationSet
  :
  public Driver
{
public:

  /** \brief Ctor */
  explicit Driver_TransientSingleEquationSet(const std::shared_ptr<OptionHandler>& option_handler)
    :
    Driver(option_handler)
  {
    option_handler->parse();
    m_option_handler = option_handler;
  }

  /** \brief Dtor */
  virtual ~Driver_TransientSingleEquationSet() = default;

  /** \brief The main routine that calls all solvers */
  virtual void solve();

  /** \brief The setup for the solvers */
  virtual void setup();

protected:

  // assembler for the steady problem
  std::shared_ptr<Assembler_Base> m_discrete_problem_assembler;
  // solver for the problem
  std::shared_ptr<NewtonSolver> m_solver;
  // time stepper for the problem
  std::shared_ptr<Stepper> m_stepper;
  // items needed for the discretization & assembly
  std::shared_ptr<Mesh> m_mesh;
  // dof managers
  std::shared_ptr<DofMapper> m_dofmapper;
  // basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // quadrature rule (domain)
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  // quadrature rule (boundary)
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature_bdry;
  // IO object
  std::shared_ptr<IOBase> m_io;
  // solution
  std::shared_ptr<FEVector> m_U;
  // gathered solution for output
  std::shared_ptr<FEArray<double>::CellBasis> m_U_gather;
  // exact solution
  std::shared_ptr<FEVector> m_U_exact;
  // IO
  std::shared_ptr<GlobalGather> m_gather;
  // underlying problem
  std::shared_ptr<Problem> m_problem;
  // other options/flags
  bool m_xml;
  bool m_write_solution_matlab;
  bool m_write_solution_vtk;
  bool m_compute_convergence_errors;

  double m_final_time;

  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("xml-out",po::value<bool>(&m_xml)->default_value(false),"Write in XML format for VTK")
      ("compute-convergence-errors",po::value<bool>(&m_compute_convergence_errors)->default_value(false),"Compute convergence error")
      ("write-output-vtk",po::value<bool>(&m_write_solution_vtk)->default_value(true),"Write the solution for output in ParaView");
  }

  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Driver_TransientSingleEquationSet_HPP__ */
