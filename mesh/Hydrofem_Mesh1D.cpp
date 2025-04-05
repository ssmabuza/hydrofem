// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_Mesh1D.hpp"
#include "Hydrofem_ElementShapeTools.hpp"

namespace hydrofem
{

Mesh1D::Mesh1D(const Mesh1D::mesh_info_type& mesh_info)
{
  num_dims() = 1;
  const std::size_t& ndims = std::get<0>(mesh_info);
  assert(num_dims()==static_cast<int>(ndims));
  const std::vector<double>& sides = std::get<1>(mesh_info);
  const double& L_l = sides.at(0);
  const double& L_r = sides.at(1);
  const std::vector<std::size_t>& sizes = std::get<2>(mesh_info);
  const std::size_t& nx = sizes.at(0);
  generateUniformMesh(nx,L_l,L_r);
  m_mesh_type = MeshType::typeLineMesh;
}


Mesh1D::Mesh1D(const int nx, const double L_l, const double L_r)
{
  num_dims() = 1;
  generateUniformMesh(nx,L_l,L_r);
  m_mesh_type = MeshType::typeLineMesh;
}

Mesh1D::~Mesh1D() {}

void Mesh1D::generateUniformMesh(const int nx, const double L_l, const double L_r)
{
  x0 = L_l;
  xf = L_r;
  
  auto& _points = points();
  auto& _elems = elems();
  
  _points.reserve(nx+1);
  const double delta_x = (L_r-L_l)/nx;
  for (int node_ind = 0; node_ind < nx+1; ++node_ind)
    _points.emplace_back(L_l+delta_x*node_ind);
  
  _elems.clear();
  _elems.resize(nx);
  for (auto elem_ind = 0; elem_ind < nx; ++elem_ind)
  {
    int node1 = elem_ind;
    int node2 = elem_ind+1;
    _elems[elem_ind] = std::make_shared<Element1D>(node1,node2);
    _elems[elem_ind]->orientation = +1;
    _elems[elem_ind]->m_domain_ind = 0;
    _elems[elem_ind]->m_edges.at(0) = elem_ind;
  }
  
  setAreaEachElement();
  setEdges();
  setStencil();
  setNeighbors();
}

std::shared_ptr<Mesh> Mesh1D::refineMesh()
{
  return std::make_shared<Mesh1D>(2*numOfElements(),x0,xf);
}

void Mesh1D::setEdges()
{
  auto& _elems = elems();
  auto& _edges = edges();
  _edges.resize(_elems.size());
  for (std::size_t edge_ind = 0; edge_ind < _edges.size(); ++edge_ind)
  {
    _edges.at(edge_ind) = std::make_shared<Edge1D>();
    _edges.at(edge_ind)->m_nodes[0] = edge_ind;
    _edges.at(edge_ind)->m_nodes[1] = edge_ind+1;
    _edges.at(edge_ind)->m_elems[0] = edge_ind-1;
    _edges.at(edge_ind)->m_elems[1] = edge_ind+1;
    if (int(edge_ind)==0)
      _edges.at(edge_ind)->m_is_boundary = true;
    else if (int(edge_ind)==int(_edges.size()-1))
    {
      _edges.at(edge_ind)->m_elems[1] = -1;
      _edges.at(edge_ind)->m_is_boundary = true;
    }
  }
  setLengthEachEdge();
}

void Mesh1D::writeMeshToFile_MATLAB(const std::string filename_prefix) const
{
  std::string x_filename (filename_prefix);
  x_filename += "x.mat";

  std::shared_ptr<std::ofstream> outX = std::make_shared<std::ofstream>(x_filename);  
  if (!outX->good())
  {
    std::cout << "Error in void Mesh1D::WriteMeshToFile_MATLAB" << std::endl;
    std::cout << "Can not open output file " << x_filename << std::endl;
    exit(1);
  }
  
  for (int point_ind = 0; point_ind < numOfPoints(); ++point_ind)
    *outX << getPoint(point_ind).x() << " ";
  *outX << std::endl;
}

double Mesh1D::evalEdgeLength(const int edge_ind)
{
  const auto& pointInd = getEdgeGlobalIndexes(edge_ind);
  return evalDistance(getPoint(pointInd.at(0)),getPoint(pointInd.at(1)));
}

SPoint Mesh1D::evalEdgeNormal(const int elem_ind, const int /*local_edge_ind*/) const
{
  if (elem_ind == numOfElements()-1)
    return SPoint(1.0);
  else if (elem_ind == 0)
    return SPoint(-1.0);
  else 
    return SPoint(0.0);
}

double Mesh1D::evalElementArea(const int elem_ind)
{
  const auto& pointInd = getEdgeGlobalIndexes(elem_ind);
  return evalDistance(getPoint(pointInd.at(0)),getPoint(pointInd.at(1)));
}

void Mesh1D::setAreaEachElement()
{
  for (int cell_ind = 0; cell_ind < numOfElements(); ++cell_ind)
    elem(cell_ind).area = evalElementArea(cell_ind);
}

void Mesh1D::setLengthEachEdge()
{
  auto& _edges = edges();
  for (int edge_ind = 0; edge_ind < numOfEdges(); ++edge_ind)
    _edges[edge_ind]->length = evalEdgeLength(edge_ind);
}

void Mesh1D::setNeighbors()
{
  auto& _neighbors = neighbors();
  _neighbors.resize(numOfElements());
  
  // go over all the elements 
  for (int elem_ind = 0; elem_ind < numOfElements(); ++elem_ind)
  {
    const auto& _elem = getElement(elem_ind);
    // go over all the edges 
    _neighbors[elem_ind].resize(2,-1);
    for (std::size_t ledge_ind(0); ledge_ind < _neighbors[elem_ind].size(); ++ledge_ind)
    {
      const auto& _edge = getEdge(_elem.m_edges[ledge_ind]);
      if (elem_ind == _edge.m_elems[0])
        _neighbors[elem_ind][ledge_ind] = _edge.m_elems[1];
      else 
        _neighbors[elem_ind][ledge_ind] = _edge.m_elems[0];
    }
  }
}

void Mesh1D::setStencil()
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
      if (node_ind == getElement(elem_ind).m_nodes.at(0))
      {
        is_in_stencil = true;
        node_in_elem = 0;
      } else if (node_ind == getElement(elem_ind).m_nodes.at(1)) {
        is_in_stencil = true;
        node_in_elem = 1;
      }
        
      if (is_in_stencil || node_in_elem >= 0)
        _stencil[node_ind][elem_ind] = node_in_elem;
    }
  }
}

}
// end namespace hydrofem
