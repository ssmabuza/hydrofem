// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_DofMapperFactory.hpp"

#include "Hydrofem_HGrad_DofMapper_Line.hpp"
#include "Hydrofem_HGrad_DofMapper_Triangle.hpp"
#include "Hydrofem_HGrad_DofMapper_Quadrilateral.hpp"

namespace hydrofem
{

std::shared_ptr<DofMapper> 
DofMapperFactory::buildDofMapper(const std::shared_ptr<Mesh>& mesh) const
{
  return buildDofMapper(mesh,m_basis_order);
}

std::shared_ptr<DofMapper>
DofMapperFactory::buildDofMapper(const std::shared_ptr<Mesh> &mesh, const int basis_order)
{
  std::shared_ptr<DofMapper> dofmapper;
  MeshType mesh_type = mesh->meshType();
  if (mesh_type == MeshType::typeLineMesh)
  {
    dofmapper = std::make_shared<HGrad_DofMapper_Line>(mesh,basis_order);
  }
  else if (mesh_type == MeshType::typeTriangularMesh)
  {
    dofmapper = std::make_shared<HGrad_DofMapper_Triangle>(mesh,basis_order);
  }
  else if (mesh_type == MeshType::typeQuadrilateralMesh)
  {
    dofmapper = std::make_shared<HGrad_DofMapper_Quadrilateral>(mesh,basis_order);
  }
  else {
    
    std::stringstream ss;
    ss << "Error in DofMapperFactory::buildDofMapper, invalid mesh provided.\n";
    throw std::logic_error(ss.str());
    
  }
  return dofmapper;
}

}
// end namespace hydrofem
