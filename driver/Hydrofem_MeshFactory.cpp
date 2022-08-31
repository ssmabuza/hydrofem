// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_MeshFactory.hpp"

#include "Hydrofem_Mesh1D.hpp"
#include "Hydrofem_Mesh2D_Tri.hpp"
#include "Hydrofem_Mesh2D_Quad.hpp"

namespace hydrofem
{

void MeshFactory::addOptionsCallback(po::options_description &config)
{
  config.add_options()
    ("ndims",po::value<int>(&m_ndims)->default_value(1),"Number of dimensions")
    ("input",po::value<std::string>(&m_input)->default_value("Inline"),"Mesh source")
    ("element-type",po::value<std::string>(&m_elem_type)->default_value("Line"),"Mesh element type")
    ("triangulation",po::value<std::string>(&m_triangulation)->default_value("Generic"),"Triangular mesh style")
    ("write-output-matlab",po::value<bool>(&m_write_to_matlab)->default_value(false),"Write mesh to MATLAB trimesh format")
    ("x0",po::value<double>(&m_x0)->default_value(0.0),"Left value for x.")
    ("xf",po::value<double>(&m_xf)->default_value(1.0),"Right value for x.")
    ("y0",po::value<double>(&m_y0)->default_value(0.0),"Left value for y.")
    ("yf",po::value<double>(&m_yf)->default_value(1.0),"Right value for y.")
    ("z0",po::value<double>(&m_z0)->default_value(0.0),"Left value for z.")
    ("zf",po::value<double>(&m_zf)->default_value(1.0),"Right value for z.")
    ("nx",po::value<int>(&m_nx)->default_value(-10),"Number of cells in the x direction.")
    ("ny",po::value<int>(&m_ny)->default_value(-10),"Number of cells in the y direction.")
    ("nz",po::value<int>(&m_nz)->default_value(-10),"Number of cells in the z direction.")
    ("mesh-file-name",po::value<std::string>(&m_filename)->default_value("mesh.am_fmt"),"Mesh file name");
}

std::shared_ptr<Mesh>
MeshFactory::buildMesh() const
{
  std::shared_ptr<Mesh> mesh;
  
  // mesh information
  //@{
  bool invalid_dim = !((m_ndims==1)||(m_ndims==2));
  if (invalid_dim)
  {
    std::stringstream ss;
    ss << "Error in \"MeshFactory::buildMesh()\"\n. The number of dimensions = " 
       << m_ndims << " is invalid.\n Valid dimensions are 1 & 2.";
    throw std::logic_error(ss.str());
  }
  
  MeshType mesh_type;
  if (m_elem_type=="Line")
  {
    mesh_type = MeshType::typeLineMesh;
  }
  else if (m_elem_type=="Triangular")
  {
    mesh_type = MeshType::typeTriangularMesh;
  }
  else if (m_elem_type=="Quadrilateral")
  {
    mesh_type = MeshType::typeQuadrilateralMesh;
  }
  else
  {
    std::stringstream ss;
    ss << "Error in \"MeshFactory::buildMesh()\", the \"Element Type\" is invalid.";
    throw std::logic_error(ss.str());
  }
  
  TriangulationType tri_type = TriangulationType::typeGenericTriangulation;
  if (mesh_type==MeshType::typeTriangularMesh)
  {
    if (m_triangulation=="SWtoNE")
    {
      tri_type = TriangulationType::typeSWtoNE;
    }
    else if (m_triangulation=="SEtoNW")
    {
      tri_type = TriangulationType::typeSEtoNW;
    }
    else if (m_triangulation=="UnionJack")
    {
      tri_type = TriangulationType::typeUnionJack;
    }
    else if (m_triangulation=="BottomNWtoSETopSEtoNW")
    {
      tri_type = TriangulationType::typeBottomNWtoSETopSEtoNW;
    }
    else {
      std::stringstream ss;
      ss << "Error in \"MeshFactory::buildMesh()\", the \"Triangulation\" is invalid.";
      throw std::logic_error(ss.str());
    }
  }
  
  bool inline_mesh (m_input == "Inline");
  //@}
  
  if (inline_mesh)
  {
    using mesh_info_type = Mesh::mesh_info_type;
    
    std::vector<double> domain;
    std::vector<std::size_t> num_cells;
    
    if (m_ndims==1)
    {
      domain = {m_x0,m_xf};
      num_cells = {static_cast<std::size_t>(m_nx)};
    }
    else if (m_ndims==2)
    {
      domain = {m_x0,m_xf,m_y0,m_yf};
      num_cells = {static_cast<std::size_t>(m_nx),
                   static_cast<std::size_t>(m_ny)};
    }
    
    // specify mesh information
    mesh_info_type mesh_info
      = std::make_tuple(std::size_t(m_ndims),
                        std::vector<double>(domain),
                        std::vector<std::size_t>(num_cells),
                        MeshType(mesh_type),
                        TriangulationType(tri_type));
    
    if (mesh_type==MeshType::typeLineMesh)
      mesh = std::make_shared<Mesh1D>(mesh_info);
    else if (mesh_type==MeshType::typeTriangularMesh)
      mesh = std::make_shared<Mesh2D_Tri>(mesh_info);
    else if (mesh_type==MeshType::typeQuadrilateralMesh)
      mesh = std::make_shared<Mesh2D_Quad>(mesh_info);
  
  } else /* NOTE: Assumes a mesh in a file AM_FMT formart */ {
    
    assert(m_ndims > 1);
    if (mesh_type==MeshType::typeTriangularMesh)
      mesh = std::make_shared<Mesh2D_Tri>(m_filename);
    else if (mesh_type==MeshType::typeQuadrilateralMesh)
      mesh = std::make_shared<Mesh2D_Quad>(m_filename);
  }
  
  // check if mesh is built  
  assert(mesh);
  
  // Write original mesh to MATLAB if requested so we can view it solo
  if (m_write_to_matlab)
    mesh->writeMeshToFile_MATLAB("original_mesh");
    
  return mesh;
}

}
// end namespace hydrofem
