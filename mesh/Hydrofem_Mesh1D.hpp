// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Mesh1D_HPP__
#define __Hydrofem_Mesh1D_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_Edge1D.hpp"
#include "Hydrofem_Element1D.hpp"

namespace hydrofem
{

/**
 * \brief Main grid implementation for simple meshes in 1D.
 */
class Mesh1D
  :
  public Mesh
{
  
  //! Ctor
  Mesh1D() {}
  
public:
  
  using mesh_info_type = Mesh::mesh_info_type;
  
  //! Ctor
  Mesh1D(const mesh_info_type& mesh_info);
  
  //! Ctor
  Mesh1D(const int nx, const double L_l, const double L_r);
  
  //! Grid destructor
  virtual ~Mesh1D();
  
  //! Functions for generating and manipulating the grid.
  //@{
  /** Generate uniform mesh */
  void generateUniformMesh(const int nx, const double L_l, const double L_r);
  
  /** Refines and returns a uniformly refined mesh */
  std::shared_ptr<Mesh> refineMesh() override;
  //@}
  
  //! Output methods
  //@{
  /** Write mesh points into an ascii file for MATLAB plotting */
  virtual void writeMeshToFile_MATLAB(const std::string filename_prefix) const override;
  //@}

  //! \brief function evaluate normal to the local edge \p local_edge_ind  which is outward for element \p elem_ind
  virtual SPoint evalEdgeNormal(const int elem_ind, const int) const override;
  
protected:
  
  MeshType m_mesh_type = MeshType::typeLineMesh;
  
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

  //! \brief computes and sets the length of each edge in the mesh
  virtual double evalEdgeLength(const int edge_ind);

  //! \brief computes and sets the area of each element in the mesh
  virtual double evalElementArea(const int elem_ind);
  
  double x0;
  double xf;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Mesh1D_HPP__ */
