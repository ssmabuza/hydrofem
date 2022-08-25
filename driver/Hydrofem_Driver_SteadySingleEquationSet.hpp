// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Driver_SteadySingleEquationSet_HPP__
#define __Hydrofem_Driver_SteadySingleEquationSet_HPP__

#include "Hydrofem_Driver.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_GlobalGather.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

class Mesh;
class IOBase;
class FEBasis; 
class Problem;
class DofMapper;
class Quadrature;
class NewtonSolver;
class Assembler_Base;

/**
 * \brief A pure base class for the driver 
 */
class Driver_SteadySingleEquationSet
  :
  public Driver
{
public:
  
  /** \brief Ctor */
  explicit Driver_SteadySingleEquationSet(const std::shared_ptr<OptionHandler>& option_handler)
    : Driver(option_handler)
  {
    m_option_handler = option_handler;
    option_handler->parse();
  }
  
  /** \brief Dtor */
  virtual ~Driver_SteadySingleEquationSet() = default;

  /** \brief The main routine that calls all solvers */
  virtual void solve();

  /** \brief The setup for the solvers */
  virtual void setup();
  
protected:
  
  // assembler for the steady problem
  std::shared_ptr<Assembler_Base> m_discrete_problem_assembler;
  // solver for the problem
  std::shared_ptr<NewtonSolver> m_solver;
  // items needed for the discretization & assembly
  std::shared_ptr<Mesh> m_mesh;
  // dof managers
  std::shared_ptr<DofMapper> m_dofmapper;
  // basis functions
  std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>> m_basis;
  // quadrature rule
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
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
  
  /** \brief options to be parsed for solver */
  virtual void addOptionsCallback(po::options_description &config)
  {
    config.add_options()
      ("xml-out",po::value<bool>(&m_xml)->default_value(false),"Write in XML format for VTK")
      ("compute-convergence-errors",po::value<bool>(&m_compute_convergence_errors)->default_value(false),"Compute convergence error")
      ("write-solution-matlab",po::value<bool>(&m_write_solution_matlab)->default_value(false),"Write the solution for output in MATLAB")
      ("write-solution-vtk",po::value<bool>(&m_write_solution_vtk)->default_value(true),"Write the solution for output in ParaView");
  }
  
  // the system input from bash file or command line
  std::shared_ptr<OptionHandler> m_option_handler;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_Driver_SteadySingleEquationSet_HPP__ */
