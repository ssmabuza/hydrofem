
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
  std::shared_ptr<hydrofem::Mesh> mesh
    = std::make_shared<hydrofem::Mesh1D>(n,0.,1.);
  ASSERT_EQ(mesh->numOfPoints(),n+1);
  ASSERT_EQ(mesh->numOfElements(),n);
  ASSERT_EQ(mesh->numOfDims(),1);
  ASSERT_EQ(mesh->meshType(),hydrofem::MeshType::typeLineMesh);
}

TEST(mesh_unit_tests,mesh2d_test_quad)
{
  const int n = 2;
  std::shared_ptr<hydrofem::Mesh> mesh
    = std::make_shared<hydrofem::Mesh2D_Quad>(n,n,0,1,0,1);
  ASSERT_EQ(mesh->numOfPoints(),(n+1)*(n+1));
  ASSERT_EQ(mesh->numOfElements(),n*n);
  ASSERT_EQ(mesh->numOfDims(),2);
  ASSERT_EQ(mesh->meshType(),hydrofem::MeshType::typeQuadrilateralMesh);
}

TEST(mesh_unit_tests,mesh2d_test_tri)
{
  const int n = 2;
  auto mesh = std::make_shared<hydrofem::Mesh2D_Tri>(n,n,0,1,0,1,hydrofem::TriangulationType::typeUnionJack);
  ASSERT_EQ(mesh->numOfPoints(),(n+1)*(n+1));
  ASSERT_EQ(mesh->numOfElements(),2*n*n);
  ASSERT_EQ(mesh->numOfDims(),2);
  ASSERT_EQ(mesh->meshType(),hydrofem::MeshType::typeTriangularMesh);
}

TEST(mesh_unit_tests,mesh2d_test_tri_2)
{
  const int n = 4;
  auto mesh = std::make_shared<hydrofem::Mesh2D_Tri>(n,n,0,1,0,1,hydrofem::TriangulationType::typeUnionJack);
  ASSERT_EQ(mesh->numOfPoints(),(n+1)*(n+1));
  ASSERT_EQ(mesh->numOfElements(),2*n*n);
  ASSERT_EQ(mesh->numOfDims(),2);
  ASSERT_EQ(mesh->meshType(),hydrofem::MeshType::typeTriangularMesh);
} 

