
#pragma once 
#include <gtest/gtest.h>

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_Mesh1D.hpp"
#include "Hydrofem_Mesh2D.hpp"
#include "Hydrofem_Mesh2D_Tri.hpp"
#include "Hydrofem_Mesh2D_Quad.hpp"


TEST(mesh_unit_tests,mesh1d_test)
{
  const int n = 2;
  auto mesh = std::make_shared<hydrofem::Mesh1D>(n);
  ASSERT_EQ(mesh->getNumNodes(),n+1);
  ASSERT_EQ(mesh->getNumElements(),n);
  ASSERT_EQ(mesh->getDim(),1);
  ASSERT_EQ(mesh->getType(),hydrofem::MeshType::LINEAR);
}

TEST(mesh_unit_tests,mesh2d_test_quad)
{
  const int n = 2;
  auto mesh = std::make_shared<hydrofem::Mesh2D_Quad>(n,n,0,1,0,1);
  ASSERT_EQ(mesh->getNumNodes(),(n+1)*(n+1));
  ASSERT_EQ(mesh->getNumElements(),n*n);
  ASSERT_EQ(mesh->getDim(),2);
  ASSERT_EQ(mesh->getType(),hydrofem::MeshType::LINEAR);
}

TEST(mesh_unit_tests,mesh2d_test_tri)
{
  const int n = 2;
  auto mesh = std::make_shared<hydrofem::Mesh2D_Tri>(n,n,0,1,0,1);
  ASSERT_EQ(mesh->getNumNodes(),(n+1)*(n+1));
  ASSERT_EQ(mesh->getNumElements(),2*n*n);
  ASSERT_EQ(mesh->getDim(),2);
  ASSERT_EQ(mesh->getType(),hydrofem::MeshType::TRIANGULAR);
}

TEST(mesh_unit_tests,mesh2d_test_tri_2)
{
  const int n = 4;
  auto mesh = std::make_shared<hydrofem::Mesh2D>(n,n,0,1,0,1);
  ASSERT_EQ(mesh->getNumNodes(),(n+1)*(n+1));
  ASSERT_EQ(mesh->getNumElements(),2*n*n);
  ASSERT_EQ(mesh->getDim(),2);
  ASSERT_EQ(mesh->getType(),hydrofem::MeshType::TRIANGULAR);
} 

