// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FunctionFE_HPP__
#define __Hydrofem_FunctionFE_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_FunctionElement.hpp"

namespace hydrofem
{

/**
 * \brief An abstract base class for scalar finite element function on a dofmapper
 * 
 * The FE function is constructed for general polynomials 
 * in 1D and 2D for quadrilaterals and triangles and concrete implementations 
 * are provided for those with ample room to add other features
 */
template <typename ScalarT>
class FunctionFE
{
protected:
  
  const typename FEArray<ScalarT>::CellBasis&                   m_coeffs;
  std::shared_ptr<DofMapper>                                       m_dofmapper;
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>> m_basis;

public:

  /// Copy constructor for vector
  FunctionFE(const typename FEArray<ScalarT>::CellBasis& tcoeffs,
             const std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>>& fe_basis,
             const std::shared_ptr<DofMapper>& tdofmapper) : m_coeffs(tcoeffs)
  {
    m_dofmapper = tdofmapper;
    m_basis = fe_basis;
  }
  
  inline int p() const { return m_dofmapper->p(); }

  /// get the mesh
  std::shared_ptr<Mesh> mesh() const
  { return m_dofmapper->mesh(); }
  
  /// get the dofmapper
  std::shared_ptr<DofMapper> dofmapper() const
  { return m_dofmapper; }
  
  /// get the view of the global coefficients
  const typename FEArray<ScalarT>::CellBasis& coeffs() const
  { return m_coeffs; }

  /// get the shape function in a given element
  virtual FunctionElement<ScalarT> create_function_element(const int i_elem) const = 0;
  
  /// calculate error
  virtual ScalarT error_l1(const std::function<ScalarT(SPoint)>& func) const = 0;

  /// calculate error
  virtual ScalarT error_l2(const std::function<ScalarT(SPoint)>& func) const = 0;

  /// calculate min
  virtual ScalarT min(int tp_) const = 0;

  /// calculate max
  virtual ScalarT max(int tp_) const = 0;

  /// set local bounds
  virtual void set_local_bounds(FunctionFE& func_min, FunctionFE& func_max, const int meth_) const = 0;

  /// set local bounds
  virtual void set_local_bounds(std::shared_ptr<FunctionFE<ScalarT>>& func_min,
                                std::shared_ptr<FunctionFE<ScalarT>>& func_max,
                                const int meth_) const = 0;
  
};

}
// end namespace hydrofem

#endif /** __Valiant_FunctionFE_HPP__ */
