// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
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

// sparse matrix type for global assembly
using FEMatrix = Eigen::SparseMatrix<double,Eigen::RowMajor>;
// vector type for global assembly
using FEVector = Eigen::VectorXd;

// base class for an operator
class Operator
{
public:
  //! The method for operation on a vector
  virtual FEVector operator*(const FEVector& vec) const = 0;
  //! The method for transposed operation on a vector
  virtual FEVector trans_mult(const FEVector& vec) const = 0;
protected:
  //! possible internal op is the sparse matrix
  std::shared_ptr<FEMatrix> m_op_internal;
};

}
// end namespace hydrofem

#endif /** __Hydrofem_EigenFEVecMat_HPP__ */
