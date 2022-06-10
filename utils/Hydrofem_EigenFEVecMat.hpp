// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_EigenFEVecMat_HPP__
#define __Hydrofem_EigenFEVecMat_HPP__

/**
 * Points, dofmappers and linear algebra objects are renamed to 
 * resemble the previous serial version of the code.
 * 
 */

#include "Hydrofem.hpp"

namespace hydrofem
{

// map type (regular GO type for now)
// using TMap = Tpetra::Map<int,int>;
// vector type for global assembly
using FEVector = Eigen::VectorXd;
// sparse matrix type for global assembly
using FEMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;

// base class for an operator
class Operator
{
public:
  //! The method for operation on a vector
  virtual FEVector operator*(const FEVector& vec) const = 0;
  //! The method for transposed operation on a vector
  virtual FEVector trans_mult(const FEVector& vec) const = 0;
};

}
// end namespace hydrofem

#endif /** __Hydrofem_EigenFEVecMat_HPP__ */
