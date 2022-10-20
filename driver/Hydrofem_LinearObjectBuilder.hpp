// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_LinearObjectBuilder_HPP__
#define __Hydrofem_LinearObjectBuilder_HPP__

#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

class LinearObjectBuilder 
{
public:

  explicit LinearObjectBuilder(const std::shared_ptr<DofMapper>& dofmapper)
  { m_dofmapper = dofmapper; }

  virtual ~LinearObjectBuilder() {}

  inline std::shared_ptr<FEVector> createVector() const 
  { return std::make_shared<FEVector>(m_dofmapper->global_ndof()); }
  
  inline std::shared_ptr<FEMatrix> createSparseMatrix() const 
  { return std::make_shared<FEMatrix>(m_dofmapper->global_ndof(),m_dofmapper->global_ndof()); }

  inline void buildSparseGraph(const std::shared_ptr<FEMatrix>& mat)
  { 
    std::vector<Eigen::Triplet<double>> mat_tripletList;
    mat_tripletList.reserve(m_dofmapper->nelements() * m_dofmapper->local_ndof() * m_dofmapper->local_ndof());
    for (int elem_ind = 0; elem_ind < m_dofmapper->nelements(); ++elem_ind)
    {
      const auto& glob_ind = m_dofmapper->getGlobDofIndexes(elem_ind);
      for (std::size_t dof_i = 0; dof_i < glob_ind.size(); ++dof_i)
        for (std::size_t dof_j = 0; dof_j < glob_ind.size(); ++dof_j)
          mat_tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[dof_i], glob_ind[dof_j], 1.0));
    }
    // assemble Jacobian
    mat->setFromTriplets(mat_tripletList.begin(), mat_tripletList.end());
    // compress the matrix
    mat->makeCompressed();
    // set all entries to zero
    mat->setZero();
    // done
    return;
  }

  inline std::shared_ptr<DofMapper> dofMapper() const 
  { return m_dofmapper; }
  
private:
  
  std::shared_ptr<DofMapper> m_dofmapper;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_LinearObjectBuilder_HPP__ */
