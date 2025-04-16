
#pragma once 

#include <gtest/gtest.h>
#include "Hydrofem_Bernstein_Basis_Line.hpp"
#include "Hydrofem_Bernstein_Basis_Triangle.hpp"

TEST(disc_fe_unit_tests,bernstein_poly_test_order_1)
{
  // Test the Bernstein basis functions for a line element
  std::vector<hydrofem::SPoint> element;
  element.push_back(hydrofem::SPoint(0.0));
  element.push_back(hydrofem::SPoint(1.0));
  hydrofem::SPoint x(0.5);
  constexpr int order = 1;
  auto phi = hydrofem::Bernstein_Basis_Line(order, 1);
  
  // Check if the result is within an acceptable range
  EXPECT_NEAR(phi(x,element), 0.5, 1e-6);
}

TEST(disc_fe_unit_tests,bernstein_poly_test_order_2)
{
  // Test the Bernstein basis functions for a triangle element
  std::vector<hydrofem::SPoint> element;
  element.push_back(hydrofem::SPoint(0.0, 0.0));
  element.push_back(hydrofem::SPoint(1.0, 0.0));
  element.push_back(hydrofem::SPoint(0.5, std::sqrt(3.0)/2.0));
  hydrofem::SPoint x(0.5, std::sqrt(3.0)/6.0);
  constexpr int order = 2;
  auto phi = hydrofem::Bernstein_Basis_Triangle(order, 0,0);
  
  // Check if the result is within an acceptable range
  EXPECT_NEAR(phi(x,element), 1.0/3.0, 1e-6);
}

#include "Hydrofem_Bernstein_Basis_Triangle.hpp"

