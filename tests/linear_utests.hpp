
#pragma once 
#include <gtest/gtest.h>

#include "Hydrofem_LinearSolvers.hpp"

TEST(linear_unit_tests,indentity_matrix)
{
  using hydrofem::FEMatrix;
  using hydrofem::FEVector;
  using hydrofem::LinearSolvers;

  auto x = std::make_shared<FEVector>(100);
  auto y = std::make_shared<FEVector>(100);
  x->setConstant(1.0);
  for (int i{0}; i < 100; ++i) {
    if (i%2==0)
      (*x)[i] *= -1.0;
  }
  
  auto A = std::make_shared<FEMatrix>(100,100);
  A->setIdentity();
  LinearSolvers::solveSystemCG(A,x,y);
  ASSERT_EQ(x->sum(), 0.0);
}


