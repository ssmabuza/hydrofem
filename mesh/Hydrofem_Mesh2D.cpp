// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Mesh2D.hpp"
#include "Hydrofem_ElementShapeTools.hpp"

namespace hydrofem
{

void Mesh2D::setNeighbors()
{
  auto& _neighbors = neighbors();
  _neighbors.resize(numOfElements());
  
  // go over all the elements 
  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    const auto& _elem = getElement(elem_ind);
    // go over all the edges 
    _neighbors[elem_ind].resize(_elem.m_edges.size());
    for (std::size_t ledge_ind(0); ledge_ind < _elem.m_edges.size(); ++ledge_ind)
    {
      const auto& _edge = getEdge(_elem.m_edges[ledge_ind]);
      if (elem_ind == _edge.m_elems[0])
        _neighbors[elem_ind][ledge_ind] = _edge.m_elems[1];
      else 
        _neighbors[elem_ind][ledge_ind] = _edge.m_elems[0];
    }
  }
}

double Mesh2D::evalElementArea(const int elem_ind)
{
  const auto& pointInd = getElementGlobalIndexes(elem_ind);
  const ObjType thisCellType = getElement(elem_ind).getElementType();
  double area = -1.0;
  switch (thisCellType)
  {
    case typeTriangle:
    {
      const auto point1 = SPoint(getPoint(pointInd.at(0)).x(),getPoint(pointInd.at(0)).y());
      const auto point2 = SPoint(getPoint(pointInd.at(1)).x(),getPoint(pointInd.at(1)).y());
      const auto point3 = SPoint(getPoint(pointInd.at(2)).x(),getPoint(pointInd.at(2)).y());
      
      area = computeTriangleArea(point1,
                                 point2,
                                 point3);
    }
    break;

    case typeQuadrilateral:
    {
      const auto point1 = SPoint(getPoint(pointInd.at(0)).x(),getPoint(pointInd.at(0)).y());
      const auto point2 = SPoint(getPoint(pointInd.at(1)).x(),getPoint(pointInd.at(1)).y());
      const auto point3 = SPoint(getPoint(pointInd.at(2)).x(),getPoint(pointInd.at(2)).y());
      const auto point4 = SPoint(getPoint(pointInd.at(3)).x(),getPoint(pointInd.at(3)).y());
      
      area = computeTriangleArea(point1,
                                 point2,
                                 point3) +
             computeTriangleArea(point3,
                                 point4,
                                 point1);
    }
    break;
    default:
    { }
  }
  return area;
}

SPoint Mesh2D::evalEdgeNormal(const int elem_ind, const int local_edge_ind) const
{
  const SPoint elem_center = evalElementCenterPoint(elem_ind);
  const int& edge_ind = getElement(elem_ind).m_edges.at(local_edge_ind);
  const SPoint edge_center = evalEdgeCenterPoint(edge_ind);
  const SPoint test = edge_center - elem_center;
  const SPoint tangent = getPoint(getEdge(edge_ind).m_nodes.at(0)) - getPoint(getEdge(edge_ind).m_nodes.at(1));
  SPoint normal = SPoint(tangent.y(),-tangent.x());
  normal.normalize();
  if (normal * test < 0.0)
    normal.reverse();
  return normal;
}

int Mesh2D::getIndexOfTheEdge(const int node1, const int node2)
{
  if (m_search_edge.empty())
  {
    m_search_edge.resize(numOfPoints());
    return -1;
  }
  const int node = (node1 < node2) ? node1 : node2;
  if (static_cast<int>(m_search_edge.size()) < node+1)
    m_search_edge.resize(node+1);
  for (std::size_t i = 0; i < m_search_edge[node].size(); ++i)
    if (m_search_edge[node][i].isTheSameEdge(node1,node2))
      return m_search_edge[node][i].m_edge_ind;
  return -1;
}

int Mesh2D::addEdge(const int node1, const int node2)
{
  const int edge_ind = getIndexOfTheEdge(node1,node2);
  if (edge_ind != -1)
    return edge_ind;
  auto& _edges = edges();
  _edges.push_back(std::make_shared<Edge2D>(node1,node2));
  int node = node1;
  if (node > node2)
    node = node2;
  if (static_cast<int>(m_search_edge.size()) < node+1)
    m_search_edge.resize(node+1);
  m_search_edge[node].push_back(SearchEdgeStruct(node1,node2,_edges.size()-1));
  return static_cast<int>(_edges.size())-1;
}

void Mesh2D::setStencil()
{
  auto& _stencil = stencils();
  _stencil.clear();
  _stencil.resize(numOfPoints());
  
  for (int node_ind = 0; node_ind < numOfPoints(); ++node_ind)
  {
    bool is_in_stencil;
    for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
    {
      
      // check if the element is in the stencil
      int node_in_elem = -1;
      is_in_stencil = false;
      
      if (getElement(elem_ind).getElementType() == typeTriangle)
      {
        
        if (node_ind == getElement(elem_ind).m_nodes.at(0))
        {
          is_in_stencil = true;
          node_in_elem = 0;
        } else if (node_ind == getElement(elem_ind).m_nodes.at(1)) {
          is_in_stencil = true;
          node_in_elem = 1;
        } else if (node_ind == getElement(elem_ind).m_nodes.at(2)) {
          is_in_stencil = true;
          node_in_elem = 2;
        } else {
          is_in_stencil = false;
          node_in_elem = -1;
        }
        
      } else if (getElement(elem_ind).getElementType() == typeQuadrilateral) {
        
        if (node_ind == getElement(elem_ind).m_nodes.at(0))
        {
          is_in_stencil = true;
          node_in_elem = 0;
        } else if (node_ind == getElement(elem_ind).m_nodes.at(1)) {
          is_in_stencil = true;
          node_in_elem = 1;
        } else if (node_ind == getElement(elem_ind).m_nodes.at(2)) {
          is_in_stencil = true;
          node_in_elem = 2;
        } else if (node_ind == getElement(elem_ind).m_nodes.at(3)) {
          is_in_stencil = true;
          node_in_elem = 3;
        } else {
          is_in_stencil = false;
          node_in_elem = -1;
        }
        
      } else {
        
        std::cerr << "Error in Mesh2D::setStencil(), element type not valid!" << std::endl;
        exit(1);
        
      }
      
      if (is_in_stencil || node_in_elem >= 0)
        _stencil[node_ind][elem_ind] = node_in_elem;
    }
  }
}

void Mesh2D::setEdges()
{
  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    auto& _elem = elem(elem_ind);
    const auto& pointInd = getElementGlobalIndexes(elem_ind);
    if (getElement(elem_ind).getElementType() == typeTriangle)
    {
      _elem.m_edges[0] = addEdge(pointInd[1], pointInd[2]);
      _elem.m_edges[1] = addEdge(pointInd[2], pointInd[0]);
      _elem.m_edges[2] = addEdge(pointInd[0], pointInd[1]);
    }
    else if (getElement(elem_ind).getElementType() == typeQuadrilateral)
    {
      _elem.m_edges[0] = addEdge(pointInd[0], pointInd[1]);
      _elem.m_edges[1] = addEdge(pointInd[1], pointInd[2]);
      _elem.m_edges[2] = addEdge(pointInd[2], pointInd[3]);
      _elem.m_edges[3] = addEdge(pointInd[3], pointInd[0]);
    }
    else
    {
      std::cout << "In Mesh2D::setEdges() : unsupported element type!!!" << std::endl;
      exit(1);
    }
  }

  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    const auto& _elem = getElement(elem_ind);
    for (std::size_t loc_edge_ind = 0; loc_edge_ind < _elem.m_edges.size(); ++loc_edge_ind)
    {
      const int& edge_ind = _elem.m_edges[loc_edge_ind];
      auto& _edge = edge(edge_ind);
      if (_edge.m_elems[0] < 0)
        _edge.m_elems[0] = elem_ind;
      else
        _edge.m_elems[1] = elem_ind;
    }
  }

  for (int edge_ind = 0; edge_ind < numOfEdges(); ++edge_ind)
    if (getEdge(edge_ind).m_elems[1] < 0)
      edge(edge_ind).m_is_boundary = true;
}

double Mesh2D::
evalEdgeLength(const int edge_ind)
{
  const auto& pointInd = getEdgeGlobalIndexes(edge_ind);
  const double length = evalDistance(getPoint(pointInd[0]), getPoint(pointInd[1]));
  return length;
}

void Mesh2D::
setLengthEachEdge()
{
  auto& _edges = edges();
  for (int edge_ind = 0; edge_ind < numOfEdges(); ++edge_ind)
    _edges[edge_ind]->length = evalEdgeLength(edge_ind);
}

void Mesh2D::
setAreaEachElement()
{
  for (int cell_ind = 0; cell_ind < numOfElements(); ++cell_ind)
    elem(cell_ind).area = evalElementArea(cell_ind);
}

void Mesh2D::
writeMeshToFile_GMV(const std::string GMVFilename) const
{
  std::shared_ptr<std::ofstream> outr = std::make_shared<std::ofstream>(GMVFilename);
  if (!outr->good())
  {
    std::cout << "Error in void Mesh2D::writeMeshToFile_GMV" << std::endl;
    std::cout << "Can not open output file " << GMVFilename << std::endl;
    exit(1);
  }

  *outr << "gmvinput ascii" << std::endl;
  *outr << "nodev " << numOfPoints() << std::endl;
  for (int point_ind = 0 ; point_ind < numOfPoints(); ++point_ind)
    *outr << getPoint(point_ind).x() << " " << getPoint(point_ind).y() << " " <<  0.0 << " " << std::endl;

  *outr << std::endl;
  *outr << "cells  " << numOfElements() << std :: endl;
  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    if (getElement(elem_ind).getElementType() == typeTriangle)
    {
      *outr << "tri 3" << std::endl;
      for (int i = 0; i < 3; ++i)
        *outr << getElement(elem_ind).m_nodes[i]+1 << " ";
    }
    else if (getElement(elem_ind).getElementType() == typeQuadrilateral)
    {
      *outr << "quad 4" << std::endl;
      for (int i = 0; i < 4; ++i)
        *outr << getElement(elem_ind).m_nodes[i]+1 << " ";
    }
    else
    {
      std::cout << "In Mesh2D::writeMeshToFile_GMV : Unsupported element type!!!\n";
      exit(1);
    }
    *outr << std::endl;
  }

  int n_materials = 0;
  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    if (n_materials < getElement(elem_ind).m_domain_ind + 1)
      n_materials = getElement(elem_ind).m_domain_ind + 1;
  }

  if (n_materials == 0)
  {
    *outr << "materials " << numOfElements() << "  0 \n";
    for (int i = 0; i < numOfElements(); ++i)
      *outr << "mat" << i+1 << "\n";
    for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
      *outr << elem_ind+1 << " ";
  }
  else
  {
    *outr << "materials " << n_materials << "  0 \n";
    for (int i = 0; i < n_materials; ++i)
      *outr << "mat" << i+1 << std::endl;
    for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
      *outr << getElement(elem_ind).m_domain_ind + 1 << " ";
  }
  *outr << "\nendgmv\n";
}

void Mesh2D::writeMeshToFile_MATLAB(const std::string FilenamePrefix) const
{
  std::string FilenameX(FilenamePrefix);
  std::string FilenameY(FilenamePrefix);
  std::string FilenameElem(FilenamePrefix);

  FilenameX = FilenameX + "_x.mat";
  FilenameY = FilenameY + "_y.mat";
  if (m_mesh_type==typeTriangularMesh)
    FilenameElem = FilenameElem + "_tri.mat";
  else if (m_mesh_type==typeQuadrilateralMesh)
    FilenameElem = FilenameElem + "_quad.mat";
  else
    throw std::logic_error("Error in Mesh2D::writeMeshToFile_MATLAB, the triangulation of the mesh is not valid.");

  auto outX = std::make_shared<std::ofstream>(FilenameX);
  if (!outX->good())
  {
    std::cout << "Error in void Mesh2D::writeMeshToFile_MATLAB" << std::endl;
    std::cout << "Can not open output file " << FilenameX << std::endl;
    exit(1);
  }

  auto outY = std::make_shared<std::ofstream>(FilenameY);
  if (!outY->good())
  {
    std::cout << "Error in void Mesh2D::writeMeshToFile_MATLAB" << std::endl;
    std::cout << "Can not open output file " << FilenameY << std::endl;
    exit(1);
  }
  
  auto outElem = std::make_shared<std::ofstream>(FilenameElem);
  if (!outElem->good())
  {
    std::cout << "Error in void Mesh2D::writeMeshToFile_MATLAB" << std::endl;
    std::cout << "Can not open output file " << FilenameElem << std::endl;
    exit(1);
  }

  const int nPoints = numOfPoints();
  for (int point_ind = 0; point_ind < nPoints; ++point_ind)
  {
    *outX << getPoint(point_ind).x() << " ";
    *outY << getPoint(point_ind).y() << " ";
  }

  *outX << std::endl;
  *outY << std::endl;

  const int nElems = numOfElements();
  for (int elem_ind = 0; elem_ind < nElems; ++elem_ind)
  {
    const auto& _elem = getElement(elem_ind);
    const auto& pointInd = getElementGlobalIndexes(elem_ind);
    if (_elem.getElementType() == typeTriangle)
      *outElem << pointInd[0]+1 << " " << pointInd[1]+1 << " " <<  pointInd[2]+1 << std::endl;
    else if (_elem.getElementType() == typeQuadrilateral)
      *outElem << pointInd[0]+1 << " " << pointInd[1]+1 << " " << pointInd[2]+1 << " " << pointInd[3]+1 << std::endl;
  }
}

}
// end namespace hydrofem
