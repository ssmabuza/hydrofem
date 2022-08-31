// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_Driver_SteadySingleEquationSet.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_DofMapper.hpp"

#include "Hydrofem_IOFactory.hpp"
#include "Hydrofem_MeshFactory.hpp"
#include "Hydrofem_DofMapperFactory.hpp"

#include "Hydrofem_Quadrature_Tri.hpp"

#include "Hydrofem_Bernstein_Basis_Triangle.hpp"
#include "Hydrofem_Bernstein_Basis_Quadrilateral.hpp"

#include "Hydrofem_LinearSolvers.hpp"

#include "Hydrofem_GlobalGather.hpp"
#include "Hydrofem_Simple_IO.hpp"

#include "Hydrofem_ComputeError.hpp"

#include "Hydrofem_NewtonSolver.hpp"
#include "Hydrofem_AssemblerFactory.hpp"

#include "Hydrofem_ProblemFactory.hpp"

#include "Hydrofem_BC_Scalar.hpp"

namespace hydrofem 
{

void Driver_SteadySingleEquationSet::setup()
{
  // form the problem mesh 
  MeshFactory mesh_factory(m_option_handler);
  m_mesh = mesh_factory.buildMesh();
  m_write_solution_matlab = mesh_factory.writeOutputToMATLAB();
  // form the DOF manager
  DofMapperFactory dofmapper_factory(m_option_handler);
  m_dofmapper = dofmapper_factory.buildDofMapper(m_mesh);
  // form the basis (either triangle P1, or quad Q1)
  m_basis = Bernstein_Basis_Triangle::buildElementBasis(1);
//   m_basis = Bernstein_Basis_Quadrilateral::buildElementBasis(1);
  // form the numerical quadrature
  m_quadrature = std::make_shared<std::vector<std::shared_ptr<Quadrature>>>(m_mesh->numOfElements());
  for (int elemInd = 0; elemInd < m_mesh->numOfElements(); ++elemInd)
    m_quadrature->at(elemInd) = std::make_shared<Quadrature_Tri>(m_dofmapper->p(),m_mesh->getElementVertices(elemInd));
  
  m_gather = std::make_shared<GlobalGather>(m_dofmapper);
  
  // build the physical/continuous problem from the factory
  auto prob_factory = std::make_shared<ProblemFactory>(m_option_handler);
  // build problem without initalization
  m_problem = prob_factory->build();
  auto bc = std::make_shared<BC_Scalar>(m_mesh);
  bc->initializeBoundaryPointsToDirichletEverywhere();
  m_problem->setBoundaryCondition(bc);

  // set the mesh in the problem's BC
  // m_problem->getBoundaryCondition()->setMesh(m_mesh);
  // initialize the problem
  m_problem->init();
  
  // create IO
  if (m_write_solution_vtk)
  {
    auto iof = IOFactory();
    m_io = iof.buildIO(m_dofmapper,m_basis,m_xml,true);
    m_io->setDIRName("./");
    m_io->setFieldName(/*"u"*/ m_problem->dofNames().at(0));
    m_io->setFinalTime(0.0);
  }
  
  auto assembler_factory = std::make_shared<AssemblerFactory>(m_option_handler,m_problem,m_dofmapper,m_basis,m_quadrature);
  m_discrete_problem_assembler = assembler_factory->build();
  m_solver = std::make_shared<NewtonSolver>(m_option_handler,m_discrete_problem_assembler);
  m_solver->initialize();
  m_U = std::make_shared<FEVector>(m_dofmapper->global_ndof()); m_U->setZero();
  m_U_exact = std::make_shared<FEVector>(m_dofmapper->global_ndof()); m_U_exact->setZero();
  m_U_gather = std::make_shared<FEArray<double>::CellBasis>(m_dofmapper->nelements(),m_dofmapper->local_ndof());
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Solver Step //////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

void Driver_SteadySingleEquationSet::solve()
{
  // solve the problem
  m_solver->solve(m_U,m_U);
  // get the solution
  m_gather->doGather(*m_U_gather,m_U);
  // write the solution to file in paraview
  if (m_write_solution_vtk)
    m_io->writeSolutionTofile(*m_U_gather,0.0);
  // write the solution to MATLAB
  if (m_write_solution_matlab)
    Simple_IO<FEVector>::writeData(false,m_U->size(),m_problem->dofNames().at(0) + ".mat",*m_U);
  // compute error
  if (m_compute_convergence_errors && m_problem->getExactSolution())
  {
    const auto u_e_ptr = m_problem->getExactSolution();
    auto err = ComputeError::evaluateScalarFieldError(m_basis,m_quadrature,m_dofmapper,m_U,u_e_ptr,nullptr,Norm::L1);
    std::cout << "L1 Error = " << err << std::endl;
    err = ComputeError::evaluateScalarFieldError(m_basis,m_quadrature,m_dofmapper,m_U,u_e_ptr,nullptr,Norm::L2);
    std::cout << "L2 Error = " << err << std::endl;
    err = ComputeError::evaluateScalarFieldError(m_basis,m_quadrature,m_dofmapper,m_U,u_e_ptr,nullptr,Norm::Linf);
    std::cout << "Linf Error = " << err << std::endl;
  }
}

}
// end namespace hydrofem
