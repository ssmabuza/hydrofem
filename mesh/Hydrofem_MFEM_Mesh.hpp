// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_MFEM_Mesh_HPP__
#define __Hydrofem_MFEM_Mesh_HPP__

#ifdef HYDROFEM_USE_MFEM

#include <mfem/mesh/pmesh.hpp>

namespace hydrofem
{


class ParMesh {

public:

  ParMesh() = default;

  explicit ParMesh(const std::shared_ptr<mfem::ParMesh>& mesh) : m_mesh(mesh) {}

  ParMesh(const ParMesh&) = default;

  ParMesh& operator=(const ParMesh&) = default;

  virtual ~ParMesh() = default;

  std::shared_ptr<mfem::ParMesh> getUnderlyingMesh() const { return m_mesh; }

  //! \brief get the number of dimensions
  virtual int numOfDims() const
  { return m_mesh->SpaceDimension(); }
  
  //! \brief get number of owned points
  virtual inline int numOfPoints() const
  { return m_mesh->GetNV(); }
  
  //! \brief get number of owned elements
  virtual inline int numOfElements() const
  { return m_mesh->GetNE(); }

  //! \brief get number of owned edges  
  virtual inline int numOfEdges() const
  { return m_mesh->GetNEdges(); }
  




private:    

  std::shared_ptr<mfem::ParMesh> m_mesh;


};


}
// end namespace hydrofem

#endif /** HYDROFEM_USE_MFEM */

#endif /** __Hydrofem_MFEM_Mesh_HPP__ */