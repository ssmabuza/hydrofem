
#pragma once 
#include <gtest/gtest.h>

#include "Hydrofem_LinearSolvers.hpp"

TEST(linear_unit_tests,identity_matrix_test_cg)
{
  const int n = 100;
  // Test simple matrix and vector functionality
  auto A = std::make_shared<hydrofem::FEMatrix>(n,n);
  A->setIdentity();
  
  auto x = std::make_shared<hydrofem::FEVector>(n);
  x->setConstant(1.0);
  
  auto y = std::make_shared<hydrofem::FEVector>(n);
  hydrofem::LinearSolvers::solveSystemCG(A,x,y);
  
  for (int i = 0; i < n; ++i)
    ASSERT_NEAR((*x)[i],(*y)[i],1.0e-10);
}

TEST(linear_unit_tests,identity_matrix_test_bicgstab)
{
  const int n = 100;
  // Test simple matrix and vector functionality
  auto A = std::make_shared<hydrofem::FEMatrix>(n,n);
  A->setIdentity();
  
  auto x = std::make_shared<hydrofem::FEVector>(n);
  x->setConstant(1.0);
  
  auto y = std::make_shared<hydrofem::FEVector>(n);
  hydrofem::LinearSolvers::solveSystemBiCGSTAB(A,x,y);
  
  for (int i = 0; i < n; ++i)
    ASSERT_NEAR((*x)[i],(*y)[i],1.0e-10);
//  ASSERT_EQ(*x,*y);
}

TEST(linear_unit_tests,identity_matrix_test_gmres)
{
  const int n = 100;
  // Test simple matrix and vector functionality
  auto A = std::make_shared<hydrofem::FEMatrix>(n,n);
  A->setIdentity();
  
  auto x = std::make_shared<hydrofem::FEVector>(n);
  x->setConstant(1.0);
  
  auto y = std::make_shared<hydrofem::FEVector>(n);
  hydrofem::LinearSolvers::solveSystemGMRES(A,x,y);
  
  for (int i = 0; i < n; ++i)
    ASSERT_NEAR((*x)[i],(*y)[i],1.0e-10);
//  ASSERT_EQ(*x,*y);
}
