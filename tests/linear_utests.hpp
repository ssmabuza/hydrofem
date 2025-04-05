
#pragma once 
#include <gtest/gtest.h>

#include "Hydrofem_LinearSolvers.hpp"

TEST(linear_unit_tests,bernstein_poly_test_order_1)
{

  hydrofem::FEVector x(100),y(100);
  x.setConstant(1.0);
  for (int i=0; i<100; ++i)
  {
    if (i%2==0)
      x[i] *= -1.0;
  }
  
  std::shared_ptr<hydrofem::FEMatrix> A(100,100);
  A.setIdentity();
  
  hydrofem::LinearSolvers::solveSystemCG(A,x,y);



}
