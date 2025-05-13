// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_DofMapper_HPP__
#define __Hydrofem_DofMapper_HPP__

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

class DofMapper
{
protected:

  //! the number of equations
  int m_neq;
  
  //! the order of the ansatz
  int m_p;
  
  //! non-static local number of degrees of freedom
  int m_lndof;
  
  //! total number of points in the whole mesh (including ghosts)
  int m_npoints;
  
  //! total number of edges in the proc (including ghosts)
  int m_nedges;
  
  //! total number of elements in the proc (including ghosts)
  int m_nelems;
  
  //! total number of global degrees of freedom
  int m_gndof;
    
  //! the local dof indexes
  std::vector<int> m_loc_indexes;
  
  //! the global dof indexes
  std::vector<std::vector<int>> m_glob_indexes;
  
  //! mesh
  std::shared_ptr<Mesh> m_mesh;

  //! Ctor (prevent call by non derived classes)
  DofMapper(const std::shared_ptr<Mesh>& mesh,
            const int order) :
                  m_p(order),
                  m_npoints(mesh->numOfPointsTotal()),
                  m_nedges(mesh->numOfEdgesTotal()),
                  m_nelems(mesh->numOfElementsTotal()),
                  m_mesh(mesh)
  { }
  
  
public:
  
  virtual ~DofMapper() {}
  
  inline int p() const { return m_p; }

  const std::vector<int>& getLocDofIndexes() const
  { return m_loc_indexes; }
  
  virtual inline const std::vector<int>& getGlobDofIndexes(const int i_elem) const
  { return m_glob_indexes.at(i_elem); }
  
  virtual inline const int& global(const int i_elem, const int local) const 
  { return m_glob_indexes.at(i_elem).at(local); }

  virtual int global_ndof() const
  { return m_gndof; }

  virtual int local_ndof() const
  { return m_lndof; }

  // number of owned elements
  virtual int nelements() const
  { return m_mesh->numOfElements(); }
  
  // number of ghosted elements
  virtual int nelementsGhosted() const
  { return m_mesh->numOfElementsGhosted(); }

  // number of owned edges
  int nedges() const
  { return m_mesh->numOfEdges(); }
  
  // number of ghosted edges
  int nedgesGhosted() const
  { return m_mesh->numOfEdgesGhosted(); }

  // number of owned mesh nodes
  int nnodes() const
  { return m_mesh->numOfPoints(); }
  
  // number of ghosted mesh nodes
  int nnodesGhosted() const
  { return m_mesh->numOfPointsGhosted(); }

  std::shared_ptr<Mesh> mesh() const
  { return m_mesh; }

  /** \brief set the number of equations for this (for linear algebra) */
  void set_numberOfEqns(const int n_eq) 
  { m_neq = n_eq; }
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_DofMapper_HPP__ */
