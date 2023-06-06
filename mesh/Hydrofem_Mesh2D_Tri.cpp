// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Mesh2D_Tri.hpp"

namespace hydrofem
{

// void Mesh2D_Tri::setEdges()
// {
//   auto& edges_ = this->edges();
//   const auto& elems_ = this->elems();

//   std::unordered_map<std::pair<int, int>, std::shared_ptr<Edge>, boost::hash<std::pair<int, int>>> edge_map;

//   // loop over all triangles
//   for (const auto& tri : elems_)
//   {
//     // loop over all edges of the triangle
//     for (int i = 0; i < 3; ++i)
//     {
//       int j = (i + 1) % 3;
//       int v0 = tri->vertex(i);
//       int v1 = tri->vertex(j);

//       // check if edge (v0,v1) already exists
//       std::pair<int, int> edge_indices(std::min(v0, v1), std::max(v0, v1));
//       auto it = edge_map.find(edge_indices);

//       if (it != edge_map.end()) // edge already exists
//       {
//         // add triangle to existing edge
//         it->second->addTriangle(tri);
//       }
//       else // edge does not exist yet
//       {
//         // create new edge and add triangle to it
//         auto edge = std::make_shared<Edge>(v0, v1);
//         edge->addTriangle(tri);
//         edges_.push_back(edge);

//         // add edge to edge_map
//         edge_map[edge_indices] = edge;
//       }
//     }
//   }
// }


Mesh2D_Tri::Mesh2D_Tri(const std::string mesh_filename): Mesh2D()
{
  num_dims() = 2;
  m_mesh_type = typeTriangularMesh;
  readMesh_AM_FMT(mesh_filename);
}

void Mesh2D_Tri::generateTriangulatedUniformRectangularMesh(const int nx, const int ny, const double L, const double H, const TriangulationType orientation)
{
  auto& points_ = this->points();
  auto& elems_ = this->elems();
  
  const int ny1 = ny + 1;
  elems_.resize(2*nx*ny);
  const double hx = L/nx;
  const double hy = H/ny;

  points_.clear();
  for (int i = 0; i <= nx; ++i)
    for (int j = 0; j <= ny; ++j)
      points_.push_back(SPoint(i*hx,j*hy));
  
  /// Column by column numeration
  int triInd = 0;
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      switch (orientation)
      {
        case typeSWtoNE:
        {
          elems_[triInd] = std::make_shared<Triangle>(j+i*ny1, j + (i+1)*ny1, (j+1) + (i+1)*ny1);
          elems_[triInd]->orientation = +1;
          ++triInd;
          elems_[triInd] = std::make_shared<Triangle>((j+1)+(i+1)*ny1, (j+1) + i*ny1, j + i*ny1);
          elems_[triInd]->orientation = +1;
          ++triInd;
        }
        break;

        case typeSEtoNW:
        {
          elems_[triInd] = std::make_shared<Triangle>((j+1)+i*ny1, j + i*ny1, j + (i+1)*ny1);
          elems_[triInd]->orientation = +1;
          ++triInd;
          elems_[triInd] = std::make_shared<Triangle>(j+(i+1)*ny1, (j+1) + (i+1)*ny1, (j+1) + i*ny1);
          elems_[triInd]->orientation = +1;
          ++triInd;

        }
        break;

        case typeUnionJack:
        {
          if ((i+j)%2 == 0)
          {
            elems_[triInd] = std::make_shared<Triangle>(j+i*ny1, j + (i+1)*ny1, (j+1) + (i+1)*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
            elems_[triInd] = std::make_shared<Triangle>((j+1)+(i+1)*ny1, (j+1) + i*ny1, j + i*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
          }
          else
          {
            elems_[triInd] = std::make_shared<Triangle>((j+1)+i*ny1, j + i*ny1, j + (i+1)*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
            elems_[triInd] = std::make_shared<Triangle>(j+(i+1)*ny1, (j+1) + (i+1)*ny1, (j+1) + i*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
          }
        }
        break;

        case typeBottomNWtoSETopSEtoNW:
        {
          if (j < ny/2)
          {
            elems_[triInd] = std::make_shared<Triangle>(j+i*ny1, j + (i+1)*ny1, (j+1) + (i+1)*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
            elems_[triInd] = std::make_shared<Triangle>((j+1)+(i+1)*ny1, (j+1) + i*ny1, j + i*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
          }
          else
          {
            elems_[triInd] = std::make_shared<Triangle>((j+1)+i*ny1, j + i*ny1, j + (i+1)*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
            elems_[triInd] = std::make_shared<Triangle>(j+(i+1)*ny1, (j+1) + (i+1)*ny1, (j+1) + i*ny1);
            elems_[triInd]->orientation = +1;
            ++triInd;
          }
        }
        break;
        
        case typeGenericTriangulation:
        {
          throw std::logic_error("The mesh generation failed, typeGenericTriangulation given!");
        }
        break;
        
      }
    }
  }
  
  setEdges();
  for (int tri_ind = 0; tri_ind < numOfElements(); ++tri_ind)
    elems_[tri_ind]->m_domain_ind = 0;
  setAreaEachElement();
  setLengthEachEdge();
  setStencil();
  setNeighbors();
}

std::shared_ptr<Mesh> Mesh2D_Tri::refineMesh()
{
  std::shared_ptr<Mesh2D_Tri> fine_grid = std::make_shared<Mesh2D_Tri>();
  
  const int nPoints_coarse = numOfPoints();
  const int nEdges_coarse = numOfEdges();
  const int nCells_coarse = numOfElements();
  const int nPoints = nPoints_coarse + nEdges_coarse;

  // Points
  auto& points_ = fine_grid->points();
  points_.resize(nPoints,SPoint(2));
  
  // a) add old points
  for (int point_ind = 0; point_ind < nPoints_coarse; ++point_ind)
    fine_grid->point(point_ind) = getPoint(point_ind);
  
  // b) add midpoints of the coarse edges
  for (int edge_ind = 0; edge_ind < nEdges_coarse; ++edge_ind)
    fine_grid->point(nPoints_coarse + edge_ind) = 0.5*(getPoint(getEdge(edge_ind).m_nodes[0]) +
                                                       getPoint(getEdge(edge_ind).m_nodes[1]));

  // Elements 
  auto& elems_ = fine_grid->elems();
  elems_.resize(nCells_coarse*4);
  
  for (int cellInd = 0; cellInd < nCells_coarse; ++cellInd)
  {
    // make sure it is really a triangle
    assert(elem(cellInd).getElementType()==typeTriangle);
    // add triangles in front
    const int n[3] = {getElement(cellInd).m_nodes[0],
                      getElement(cellInd).m_nodes[1],
                      getElement(cellInd).m_nodes[2]};
    const int m[3] = {nPoints_coarse + getElement(cellInd).m_edges[0], 
                      nPoints_coarse + getElement(cellInd).m_edges[1], 
                      nPoints_coarse + getElement(cellInd).m_edges[2]};
    for (int i0 = 0; i0 < 3; ++i0)
    {
      const int i1 = (i0+1)%3;
      const int i2 = (i1+1)%3;
      elems_.at(4*cellInd+i0) = std::make_shared<Triangle>(m[i2], m[i1], n[i0]);
    }
    elems_.at(4*cellInd+3) = std::make_shared<Triangle>(m[0], m[1], m[2]);
    
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

void Mesh2D_Tri::readMesh_AM_FMT(const std::string filename_am_fmt)
{
  
  auto& points_ = this->points();
  auto& elems_ = this->elems();
  
  auto input_stream = std::make_shared<std::ifstream>(filename_am_fmt);
  if (!input_stream->good())
    throw std::runtime_error("Error in void Mesh2D_Tri::readMesh_AM_FMT\n Can not open input file " + filename_am_fmt);
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
  int vert_indexes[3];
  for (int cellInd = 0; cellInd < nCells; ++cellInd)
  {
    *input_stream >> vert_indexes[0] >> vert_indexes[1] >> vert_indexes[2];
    elems_[cellInd] = std::make_shared<Triangle>(vert_indexes[0]-1, vert_indexes[1]-1, vert_indexes[2]-1);
  }

  for (int pointInd = 0; pointInd < nPoints; ++pointInd)
    *input_stream >> points_.at(pointInd)(0) >> points_.at(pointInd)(1);

  setEdges();
  setAreaEachElement();
  setLengthEachEdge();
  
  setStencil();
  setNeighbors();
}

void Mesh2D_Tri::updateNodes_AM_FMT(const std::string filename_am_fmt)
{
  auto& points_ = this->points();
  
  auto input_stream = std::make_shared<std::ifstream>(filename_am_fmt);
  if (!input_stream->good())
    throw std::runtime_error("Error in void Mesh2D_Tri::readMesh_AM_FMT\n Can not open input file " + filename_am_fmt);

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
  int vert_indexes[3];
  for (int cellInd = 0; cellInd < nCells; ++cellInd)
    *input_stream >> vert_indexes[0] >> vert_indexes[1] >> vert_indexes[2];

  /// now get the points
  for (int pointInd = 0; pointInd < nPoints; ++pointInd)
    *input_stream >> points_[pointInd].x() >> points_[pointInd].y();

  setAreaEachElement();
  setLengthEachEdge();
  
}

}
// end namespace hydrofem
