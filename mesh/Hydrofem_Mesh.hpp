// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Mesh_HPP__
#define __Hydrofem_Mesh_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Edge.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Element.hpp"
#include "Hydrofem_Neighbors.hpp"
#include "Hydrofem_MeshTypes.hpp"
#include "Hydrofem_NodeStencil.hpp"

namespace hydrofem
{

/**
 * \brief A base class for meshes used throughout the code.
 * 
 */
class Mesh : std::enable_shared_from_this<Mesh>
{
public:

  using Ptr = std::shared_ptr<Mesh>;
  
  //! \brief the mesh information container for a regular polygonal domain and regular uniform meshes
  using mesh_info_type = std::tuple<std::size_t,std::vector<double>,std::vector<std::size_t>,MeshType,TriangulationType>;
  
  //! \brief Ctor
  Mesh() {}
  
  //! \brief Dtor
  virtual ~Mesh() {}
  
  //! \brief refining the mesh
  virtual std::shared_ptr<Mesh> refineMesh() = 0;
  
  //! \brief get the number of dimensions
  virtual int numOfDims() const
  { return m_num_dims; }
  
  //! \brief get number of owned points
  virtual inline int numOfPoints() const
  { return m_points.size(); }
  
  //! \brief get number of owned elements
  virtual inline int numOfElements() const
  { return m_elems.size(); }

  //! \brief get number of owned edges  
  virtual inline int numOfEdges() const
  { return m_edges.size(); }

  //! \brief get number of all points
  virtual inline int numOfPointsTotal() const
  { return numOfPoints(); };
  
  //! \brief get number of all elements
  virtual inline int numOfElementsTotal() const
  { return numOfElements(); }

  //! \brief get number of all edges  
  virtual inline int numOfEdgesTotal() const
  { return numOfEdges(); }
  
  //! \brief get number of ghosted points
  virtual inline int numOfPointsGhosted() const 
  { return 0; }
  
  //! \brief get number of ghosted elements
  virtual inline int numOfElementsGhosted() const
  { return 0; }

  //! \brief get number of ghosted edges  
  virtual inline int numOfEdgesGhosted() const
  { return 0; }
  
  virtual inline int globalElemInd(const int elem_ind) const
  { assert((elem_ind >= 0) && (elem_ind < numOfElements())); return elem_ind; }
  
  //! \brief get one mesh vertex
  virtual inline std::vector<SPoint> getElementVertices(const int elem_ind) const 
  {
    std::vector<SPoint> res;
    
    if (m_elems[elem_ind]->m_element_type == typeLine)
    {
      res.resize(2,SPoint(1));
      res[0] = m_points.at(m_elems.at(elem_ind)->m_nodes[0]);
      res[1] = m_points.at(m_elems.at(elem_ind)->m_nodes[1]);
      
    } else if (m_elems[elem_ind]->m_element_type == typeTriangle) {
      
      res.resize(3,SPoint(2));
      res[0] = m_points.at(m_elems.at(elem_ind)->m_nodes[0]);
      res[1] = m_points.at(m_elems.at(elem_ind)->m_nodes[1]);
      res[2] = m_points.at(m_elems.at(elem_ind)->m_nodes[2]);
      
    } else if (m_elems[elem_ind]->m_element_type == typeQuadrilateral) {

      res.resize(4,SPoint(2));
      res[0] = m_points.at(m_elems.at(elem_ind)->m_nodes[0]);
      res[1] = m_points.at(m_elems.at(elem_ind)->m_nodes[1]);
      res[2] = m_points.at(m_elems.at(elem_ind)->m_nodes[2]);
      res[3] = m_points.at(m_elems.at(elem_ind)->m_nodes[3]);
      
    }
    else
      throw std::logic_error("Error in Mesh::getElementVertices!! Element type not valid.");
    
    return res;
  }

  virtual inline std::vector<SPoint> getEdgeVertices(const int edge_ind) const
  {
    return std::vector<SPoint>({m_points.at(m_edges.at(edge_ind)->m_nodes[0]),m_points.at(m_edges.at(edge_ind)->m_nodes[1])});
  }
  
  //! \brief get one element
  virtual inline const Element& getElement(const int elem_ind) const
  { return *(m_elems.at(elem_ind)); }
  
  //! \brief get one edge
  virtual inline const Edge& getEdge(const int edge_ind) const
  { return *(m_edges.at(edge_ind)); }

  //! \brief get one edge
  virtual inline const SPoint& getPoint(const int point_ind) const
  { return m_points.at(point_ind); }
  
  //! \brief get element neighbors
  virtual inline const Neighbors& getElementNeighbors(const int elem_ind) const
  { return m_neighbors.at(elem_ind); }
  
  //! \brief get node to element stencil 
  virtual inline const NodeStencil& getNodeToElementStencil(const int point_ind) const
  { return m_stencil.at(point_ind); }
  
  //! \brief computes the center point of an edge
  virtual inline SPoint evalEdgeCenterPoint(const int edge_ind) const
  {
    return 0.5*(m_points.at(m_edges.at(edge_ind)->m_nodes[0]) + m_points.at(m_edges.at(edge_ind)->m_nodes[1]));
  }
  //! \brief gets the center of mass in an element
  virtual SPoint evalElementCenterPoint(const int elem_ind) const
  {
    SPoint center(0.0,0.0);
    for (std::size_t i = 0; i < m_elems[elem_ind]->m_nodes.size(); ++i)
      center += m_points[m_elems[elem_ind]->m_nodes[i]];
    center /= m_elems[elem_ind]->m_nodes.size();
    return center;
  }
    
  //! \brief computes the outward normal to an element given the edge index in 2D and 3D
  virtual SPoint evalEdgeNormal(const int /*elem_ind*/, const int /*local_edge_ind*/) const
  { return SPoint(); }
  
  //! \brief returns the area of the element
  virtual inline double getElementArea(const int elem_ind) const
  { return m_elems.at(elem_ind)->area; }
  
  //! \brief gives the length of the edges  
  virtual inline double getEdgeLength(const int edge_ind) const
  { return m_edges.at(edge_ind)->length; }
  
  //! \brief returns the vertex indexes on the element
  virtual inline const std::vector<int>& getElementGlobalIndexes(const int elem_ind) const
  { return m_elems.at(elem_ind)->m_nodes; }

  //! \brief returns the vertex indexes on the element
  virtual inline const std::vector<int>& getEdgeGlobalIndexes(const int edge_ind) const
  { return m_edges.at(edge_ind)->m_nodes; }
  
  //! \brief   
  virtual inline SPoint getMaxCoords() const
  {
    SPoint res;
    if (numOfDims()==1)
    {
      double valx = std::numeric_limits<double>::min();
      for (auto point : m_points)
        valx = std::max(valx,point(0));
      res = SPoint(valx);
    } else if (numOfDims()==2) {
      double valx = std::numeric_limits<double>::min();
      double valy = std::numeric_limits<double>::min();
      for (auto point : m_points)
      {
        valx = std::max(valx,point(0));
        valy = std::max(valy,point(1));
      }
      res = SPoint(valx,valy);
    }
    return res;
  }
  
  //! \brief 
  virtual inline SPoint getMinCoords() const
  {
    SPoint res;
    if (numOfDims()==1)
    {
      double valx = std::numeric_limits<double>::max();
      for (auto point : m_points)
        valx = std::min(valx,point(0));
      res = SPoint(valx);
    } else if (numOfDims()==2) {
      double valx = std::numeric_limits<double>::max();
      double valy = std::numeric_limits<double>::max();
      for (auto point : m_points)
      {
        valx = std::min(valx,point(0));
        valy = std::min(valy,point(1));
      }
      res = SPoint(valx,valy);
    }
    return res;
  }
  
  //! \brief compute the mesh size
  virtual double getMeshSize() const 
  {
    double mesh_size = 0.0;
    for (std::size_t edge_ind = 0; edge_ind < m_edges.size(); ++ edge_ind)
      mesh_size = std::max(this->getEdgeLength(edge_ind),mesh_size);
    return mesh_size;
  }
  
  //! \brief write the mesh to files that can be read in MATLAB trimesh
  virtual void writeMeshToFile_MATLAB(const std::string /*filename_prefix*/) const
  { }
  
  //! \brief get mesh type
  inline MeshType meshType() const { return m_mesh_type; }
  
protected:

  std::vector<SPoint>& points()
  { return m_points; }  
  
  std::vector<std::shared_ptr<Element>>& elems()
  { return m_elems; }
  
  std::vector<std::shared_ptr<Edge>>& edges()
  { return m_edges; }
  
  std::vector<NodeStencil>& stencils()
  { return m_stencil; }
  
  std::vector<Neighbors>& neighbors()
  { return m_neighbors; }
  
  int& num_dims()
  { return m_num_dims; }
  
  SPoint& point(const int point_ind)
  { return m_points.at(point_ind); }
  
  Element& elem(const int elem_ind)
  { return * m_elems.at(elem_ind); }
  
  Edge& edge(const int edge_ind)
  { return * m_edges.at(edge_ind); }
  
  NodeStencil& stencil(const int point_ind)
  { return m_stencil.at(point_ind); }
  
  Neighbors& neighbor(const int elem_ind)
  { return m_neighbors.at(elem_ind); }
  
  MeshType m_mesh_type = MeshType::typeGenericMesh;
  
private:
  
  // Attributes
  //@{
  int                                   m_num_dims;
  std::vector<SPoint>                   m_points;
  std::vector<std::shared_ptr<Element>> m_elems;
  std::vector<std::shared_ptr<Edge>>    m_edges;
  std::vector<NodeStencil>              m_stencil;
  std::vector<Neighbors>                m_neighbors;
  //@}
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Mesh_HPP__ */
