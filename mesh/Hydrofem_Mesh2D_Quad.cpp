// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Mesh2D_Quad.hpp"

namespace hydrofem
{

Mesh2D_Quad::Mesh2D_Quad(const std::string mesh_filename) : Mesh2D()
{
  
  num_dims() = 2;
  m_mesh_type = typeQuadrilateralMesh;
  readMesh_AM_FMT(mesh_filename);
  
}

Mesh2D_Quad::Mesh2D_Quad(const Mesh2D_Quad::mesh_info_type& mesh_data) : Mesh2D()
{
  
  const int& ndims = std::get<0>(mesh_data);
  assert(ndims == 2);
  num_dims() = ndims;
  m_mesh_type = typeQuadrilateralMesh;
  const std::vector<double>& sides = std::get<1>(mesh_data);
  const double& L_l = sides.at(0);
  const double& L_r = sides.at(1);
  const double& H_l = sides.at(2);
  const double& H_r = sides.at(3);
  const std::vector<std::size_t>& sizes = std::get<2>(mesh_data);
  const int& nx = sizes.at(0);
  const int& ny = sizes.at(1);
  generateUniformRectangularMesh(nx,ny,L_l,L_r,H_l,H_r);

}

Mesh2D_Quad::Mesh2D_Quad(const int nx, const int ny, const double x0, const double xf, const double y0, const double yf) : Mesh2D()
{

  num_dims() = 2;
  m_mesh_type = typeQuadrilateralMesh;
  generateUniformRectangularMesh(nx,ny,x0,xf,y0,yf);
  
}

void Mesh2D_Quad::generateUniformRectangularMesh(const int nx, const int ny, const double L, const double H)
{
  auto& points_ = this->points();
  auto& elems_ = this->elems();
  
  const int nx1 = nx + 1;
  const int ny1 = ny + 1;
  points_.resize(nx1*ny1);
  elems_.resize(nx*ny);
  const double hx = L/nx;
  const double hy = H/ny;
  
  /// Column by column
  int pointInd = 0;
  for (int i = 0; i <= nx; ++i)
  {
    for (int j = 0; j <= ny; ++j)
    {
      points_[pointInd] = SPoint(double(i)*hx,double(j)*hy);
      ++pointInd;
    }
  }

  int quadInd = 0;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      elems_[quadInd] = std::make_shared<Quadrilateral>(j+i*ny1, j+(i+1)*ny1, (j+1)+(i+1)*ny1, (j+1)+i*ny1);
      ++quadInd;
    }
  }
  
  setEdges();
  setAreaEachElement();
  setLengthEachEdge();
  setStencil();
  setNeighbors();
}

std::shared_ptr<Mesh>
Mesh2D_Quad::refineMesh()
{
  std::shared_ptr<Mesh2D_Quad> fine_grid = std::make_shared<Mesh2D_Quad>();
  
  const int nPoints_coarse = numOfPoints();
  const int nEdges_coarse = numOfEdges();
  const int nCells_coarse = numOfElements();
  
  const int nPoints = nPoints_coarse + nEdges_coarse + nCells_coarse;
  
  // Points
  auto& points_ = fine_grid->points();
  points_.resize(nPoints);
  
  // a) add old points
  for (int pointInd = 0; pointInd < nPoints_coarse; ++pointInd)
    fine_grid->point(pointInd) = getPoint(pointInd);
  
  // b) add midpoints of the coarse edges
  for (int edgeInd = 0; edgeInd < nEdges_coarse; ++edgeInd)
    fine_grid->point(nPoints_coarse + edgeInd) = 0.5*(getPoint(getEdge(edgeInd).m_nodes[0]) +
                                                      getPoint(getEdge(edgeInd).m_nodes[1]));

  for (int cellInd = 0; cellInd < nCells_coarse; ++cellInd)
    fine_grid->point(nPoints_coarse + nEdges_coarse + cellInd) = evalElementCenterPoint(cellInd);

  // Elements 
  auto& elems_ = fine_grid->elems();
  elems_.resize(nCells_coarse*4);
  
  for (int cellInd = 0; cellInd < nCells_coarse; ++cellInd)
  {
    // make sure it is really a quad
    assert(elem(cellInd).getElementType()==typeQuadrilateral);
    // add triangles in front
    const int n[4] = {getElement(cellInd).m_nodes[0],
                      getElement(cellInd).m_nodes[1],
                      getElement(cellInd).m_nodes[2],
                      getElement(cellInd).m_nodes[3]};
    
    const int m[5] = {nPoints_coarse + getElement(cellInd).m_edges[0], 
                      nPoints_coarse + getElement(cellInd).m_edges[1], 
                      nPoints_coarse + getElement(cellInd).m_edges[2],
                      nPoints_coarse + getElement(cellInd).m_edges[3],
                      nPoints_coarse + nEdges_coarse + cellInd};
    
    elems_[4*cellInd+0] = std::make_shared<Quadrilateral>(n[0],m[0],m[4],m[3]);
    elems_[4*cellInd+1] = std::make_shared<Quadrilateral>(m[0],n[1],m[1],m[4]);
    elems_[4*cellInd+2] = std::make_shared<Quadrilateral>(m[4],m[1],n[2],m[2]);
    elems_[4*cellInd+3] = std::make_shared<Quadrilateral>(m[3],m[4],m[2],n[3]);
    
    for (int i = 0 ; i < 4; ++i)
    {
      elems_[4*cellInd+i]->orientation = getElement(cellInd).orientation;
      elems_[4*cellInd+i]->m_domain_ind = getElement(cellInd).m_domain_ind;
    }
  }
  
  fine_grid->setEdges();
  fine_grid->setLengthEachEdge();
  fine_grid->setAreaEachElement();
  fine_grid->setStencil();
  fine_grid->setNeighbors();
  
  return fine_grid;
}


void Mesh2D_Quad::readMesh_AM_FMT(const std::string filename_am_fmt)
{
  auto& points_ = this->points();
  auto& elems_ = this->elems();
  
  std::ifstream* input_stream = new std::ifstream(filename_am_fmt);
  if (!input_stream->good())
    throw std::runtime_error("Error in void Mesh2D_Quad::readMesh_AM_FMT\n Can not open input file " + filename_am_fmt);
  int nPoints;
  int nCells;
  *input_stream >> nPoints >> nCells;

  char next;
  while(input_stream->get(next))
  {
    if (next == '\n')
      break;
  }

  points_.resize(nPoints,SPoint(int(2)));
  elems_.resize(nCells);
  int vert_indexes[4];
  for (int cellInd = 0; cellInd < nCells; ++cellInd)
  {
    *input_stream >> vert_indexes[0] >> vert_indexes[1] >> vert_indexes[2] >> vert_indexes[3];
    elems_[cellInd] = std::make_shared<Quadrilateral>(vert_indexes[0]-1, vert_indexes[1]-1, vert_indexes[2]-1, vert_indexes[3]-1);
  }

  for (int pointInd = 0; pointInd < nPoints; ++pointInd)
    *input_stream >> points_.at(pointInd)(0) >> points_.at(pointInd)(1);

  delete input_stream;
  input_stream = 0;
  
  setEdges();
  setAreaEachElement();
  setLengthEachEdge();
  
  setStencil();
  setNeighbors();
}

void Mesh2D_Quad::updateNodes_AM_FMT(const std::string filename_am_fmt)
{
  auto& points_ = this->points();
  
  std::ifstream* input_stream = new std::ifstream(filename_am_fmt);
  if (!input_stream->good())
  {
    std::string err_msg = "Error in void Mesh2D_Quad::updateNodes_AM_FMT\n Can not open input file " + filename_am_fmt;
    throw std::runtime_error(err_msg);
  }

  int nPoints;
  int nCells;
  *input_stream >> nPoints >> nCells;
  char next;
  while (input_stream->get(next))
  {
    if (next == '\n')
      break;
  }

  /// doing nothing but going through the file to get to the up to date points, enumeration is the same
  int vert_indexes[4];
  for (int cellInd = 0; cellInd < nCells; ++cellInd)
    *input_stream >> vert_indexes[0] >> vert_indexes[1] >> vert_indexes[2] >> vert_indexes[3];

  /// now get the points
  for (int pointInd = 0; pointInd < nPoints; ++pointInd)
    *input_stream >> points_.at(pointInd)(0) >> points_.at(pointInd)(1);

  delete input_stream;
  input_stream = 0;
  
  setAreaEachElement();
  setLengthEachEdge();
  
}

}
// end namespace hydrofem
