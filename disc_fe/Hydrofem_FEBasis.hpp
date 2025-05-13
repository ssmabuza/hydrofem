// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_FEBasis_HPP__
#define __Hydrofem_FEBasis_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_GlobalConstants.hpp"

namespace hydrofem
{

/** \brief A base class for FE basis */
class FEBasis
{
public:
  
  using FEBasisPtr = std::shared_ptr<FEBasis>;
  using ptr_t = FEBasisPtr;
  
  FEBasis() {}
  
  virtual ~FEBasis() {}
  
  //!\brief evaluates the basis function at a @p point on a cell with vertices : @p element 
  virtual double operator()(const SPoint& point, const std::vector<SPoint>& element) const = 0;
  
  //!\brief evaluates the basis function gradient at a @p point on a cell with vertices : @p element 
  virtual SPoint grad(const SPoint& point, const std::vector<SPoint>& element) const = 0;

  //! \brief get the order  
  virtual int p() const = 0;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_FEBasis_HPP__ */
