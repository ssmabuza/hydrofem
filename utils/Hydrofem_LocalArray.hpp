// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_LocalArray_HPP__
#define __Hydrofem_LocalArray_HPP__

#include <unsupported/Eigen/CXX11/Tensor>

namespace hydrofem
{

//! local vector type for equation classes and local assembly
template <typename ScalarT>
using LVEC_ = Eigen::Tensor<ScalarT,1>;
//! local matrix type for equation classes and local assembly
template <typename ScalarT>
using LMAT_ = Eigen::Tensor<ScalarT,2>;
//! Array for <cell,basis>
template <typename ScalarT>
using KCellBasis = Eigen::Tensor<ScalarT,2>;
//! Array for <cell,basis,basis>
template <typename ScalarT>
using KCellBasisBasis = Eigen::Tensor<ScalarT,3>;
//! Integral array for <cell,basis>
using KIntCellBasis = Eigen::Tensor<int,2>;

template <typename ArrayT>
inline ArrayT createKArray(const int le)
{
  return ArrayT(le);
}

template <typename ArrayT>
inline ArrayT createKArray(const int le1,const int le2)
{
  return ArrayT(le1,le2);
}

template <typename ArrayT>
inline ArrayT createKArray(const int le1,const int le2,const int le3)
{
  return ArrayT(le1,le2,le3);
}

// ****************************************************************** //

template <typename ArrayT>
inline void zeroOutArray(ArrayT& array)
{
  array.setZero();
}

//! \ingroup defined Fad types as in Panzer
//@{
// scalar type
using RealType = double;
// first derivative type
// using FadType = Sacado::Fad::DFad<RealType>;
// 2nd derivative type
// using HessianType = Sacado::Fad::DFad<Sacado::Fad::SFad<RealType,1>>;
//@}

/** \brief FE local array (difference for FE basis which in only double) */
template <typename ScalarT>
struct FESolution
{
  
  using Node           = Eigen::Tensor<ScalarT,1>;
  using Cell           = Eigen::Tensor<ScalarT,1>;
  using CellQP         = Eigen::Tensor<ScalarT,2>;
  using CellQPDim      = Eigen::Tensor<ScalarT,3>;
  using CellBasis      = Eigen::Tensor<ScalarT,2>;
  using CellBasisDim   = Eigen::Tensor<ScalarT,3>;
  using CellBasisQP    = Eigen::Tensor<ScalarT,3>;
  using CellEdgeQP     = Eigen::Tensor<ScalarT,3>;
  using CellEdgeBasis  = Eigen::Tensor<ScalarT,3>;
  using CellBasisQPDim = Eigen::Tensor<ScalarT,4>;
  using CellEdgeQPDim  = Eigen::Tensor<ScalarT,4>;

  using CellBasisEqn   = Eigen::Tensor<ScalarT,3>;
  using CellQPEqn      = Eigen::Tensor<ScalarT,3>;
  using CellQPDimEqn   = Eigen::Tensor<ScalarT,4>;
  
  using CellBasisBasisEqnEqn = Eigen::Tensor<ScalarT,5>;
  
};

template <typename ScalarT>
using FEArray = FESolution<ScalarT>;

}
// end namespace hydrofem

#endif /** __Hydrofem_LocalArray_HPP__ */
