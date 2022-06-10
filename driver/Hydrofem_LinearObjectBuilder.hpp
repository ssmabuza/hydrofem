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

  inline std::shared_ptr<DofMapper> dofMapper() const 
  { return m_dofmapper; }
  
private:
  
  std::shared_ptr<DofMapper> m_dofmapper;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_LinearObjectBuilder_HPP__ */
