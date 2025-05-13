// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_AFC_Limiter_HPP__
#define __Hydrofem_AFC_Limiter_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"

namespace hydrofem
{

/**
 * An abstract interface for the iterative limiters
 * which computes the nodal limiters
 */
class AFC_Limiter
{
public:
  
  AFC_Limiter() {}

  virtual ~AFC_Limiter() {}

  /**
   * \brief main routine to build element limiters 
   */
  virtual void
  buildLimiter(const std::shared_ptr<FEVector>& res,
               const std::shared_ptr<const FEVector>& u) const = 0;
  
  /** \brief set the nodal gradients */
  virtual void 
  set_nodalGradients(const std::vector<const std::shared_ptr<const FEVector>>& /*grads*/)
  { }
  
};

}
// end namespace valiant

#endif /** __Hydrofem_AFC_Limiter_HPP__ */
