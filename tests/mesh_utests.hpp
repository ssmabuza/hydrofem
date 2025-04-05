
#pragma once

#include <gtest/gtest.h>

#include "Hydrofem_Mesh1D.hpp"

TEST(mesh_unit_tests, mesh_test_unit_line)
{
  hydrofem::Mesh1D mesh(10, 0.0, 1.0);
  double area = 0.0;
  for (int i = 0; i < mesh.numOfElements(); ++i) {
    area += mesh.getElement(i).getArea();
  }
  double expected_area = 1.0;
  // EXPECT_NEAR(area, expected_area, 1e-6);
  EXPECT_EQ(mesh.numOfElements(), 10);
  // EXPECT_EQ(mesh.numOfPoints(), 11);
  // EXPECT_EQ(mesh.numOfEdges(), 10);
  
}