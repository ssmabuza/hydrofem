// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Mesh2D_HPP__
#define __Hydrofem_Mesh2D_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Edge2D.hpp"
#include "Hydrofem_Element2D.hpp"

namespace hydrofem
{

/** 
  * \brief A mesh class in 2D for triangular or quadrilateral elements:
  * 
  *        1) mixed element of triangles and quadrilaterals 
  *           are possible but have not been used before.
  *        2) triangular meshes of am_fmt type can be read into 
  *           the mesh
  *        3) triangular meshes can be written in trimesh compatible form
  *
  *   TODO: input and output for quadrilateral meshes using GMsh
  *   TODO: output meshes into VTK for visualization in Paraview
  */
class Mesh2D
  :
  public Mesh
{
protected:

  //! \brief Ctor
  Mesh2D() : Mesh() {}
  
public:

  using mesh_info_type = Mesh::mesh_info_type;

  //! \brief Dtor
  ~Mesh2D() override = default;
  
  std::shared_ptr<Mesh> refineMesh() override = 0;

  //! \brief function evaluate normal to the local edge \p local_edge_ind  which is outward for element \p elem_ind
  SPoint evalEdgeNormal(const int elem_ind, const int local_edge_ind) const override;

  /// Input from file
  //@{
  //! \brief reads an AM_FMT type of file to form triangular mesh
  virtual void readMesh_AM_FMT(const std::string /*filename_am_fmt*/) {}
  
  //! \brief reads an AM_FMT type of file to modify existing triangular mesh
  virtual void updateNodes_AM_FMT(const std::string /*filename_am_fmt*/) {}
  //@}
  
  /// Output
  //@{
  //! \brief write a general 2D mesh into a file read in GMV
  void writeMeshToFile_GMV(const std::string GMVFilename) const;
  
  //! \brief writes a triangular mesh into files read in MATLAB for trimesh
  virtual void writeMeshToFile_MATLAB(const std::string FilenamePrefix) const override;
  
  //@}
  
protected:
  
  //! \brief function sets the stencil for all the nodes in the mesh
  void setStencil();                                                  
  
  //! \brief function that sets the neighbors of each element
  void setNeighbors();

  //! \brief function sets edges of all elements in the mesh
  void setEdges();

  //! \brief function sets the area for all elements in the mesh
  void setAreaEachElement();

  //! \brief function sets the length of all edges in the mesh
  void setLengthEachEdge();
  
private:  
  
  // utility attribute that helps to set edge info
  std::vector<std::vector<SearchEdgeStruct>> m_search_edge;

  //! \brief function finds edge with nodes node1 and node2
  int getIndexOfTheEdge(const int node1, const int node2);

  //! \brief function adds edge with nodes node1 and node2 into the mesh (array m_edges)
  //!                  and returns the index of this edge, or returns the index of the '
  int addEdge(const int node1, const int node2);

  //! \brief computes and sets the length of each edge in the mesh
  virtual double evalEdgeLength(const int edge_ind);

  //! \brief computes and sets the area of each element in the mesh
  virtual double evalElementArea(const int elem_ind);
    
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Mesh2D_HPP__ */

