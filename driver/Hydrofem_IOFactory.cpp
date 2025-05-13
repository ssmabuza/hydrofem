// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IOFactory.hpp"
#include "Hydrofem_DofMapper.hpp"

#include "Hydrofem_IOTri.hpp"
#include "Hydrofem_IOTriXML.hpp"
#include "Hydrofem_IOLine.hpp"
#include "Hydrofem_IOQuad.hpp"

#include "Hydrofem_FunctionElement.hpp"

namespace hydrofem
{

std::shared_ptr<IOBase>
IOFactory::buildIO(const std::shared_ptr<DofMapper>& dofmapper,
                   const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis, bool xml, bool steady)
{
  std::shared_ptr<IOBase> io;
  auto mesh = dofmapper->mesh();
  auto mesh_type = mesh->meshType();
  std::shared_ptr<FunctionElement<double>> m_fe_shape = std::make_shared<FunctionElement<double>>(fe_basis);
  assert(m_fe_shape);
  if (mesh_type==MeshType::typeLineMesh)
  {
    // build the io
    io = std::make_shared<IOLine>(std::dynamic_pointer_cast<HGrad_DofMapper_Line>(dofmapper),m_fe_shape);
  } else if (mesh_type==MeshType::typeTriangularMesh) {
    // build the io
    if (xml)
      io = std::make_shared<IOTriXML>(std::dynamic_pointer_cast<HGrad_DofMapper_Triangle>(dofmapper),m_fe_shape,false,steady);
    else 
      io = std::make_shared<IOTri>(std::dynamic_pointer_cast<HGrad_DofMapper_Triangle>(dofmapper),m_fe_shape);
  } else if (mesh_type==MeshType::typeQuadrilateralMesh) {
    // build the io
    io = std::make_shared<IOQuad>(std::dynamic_pointer_cast<HGrad_DofMapper_Quadrilateral>(dofmapper),m_fe_shape);
  }
  return io;
}

}
// end namespace hydrofem
