// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_AFC_LPSLimiter_HPP__
#define __Hydrofem_AFC_LPSLimiter_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_LocalArray.hpp"

#include "Hydrofem_AFC_Limiter.hpp"

namespace hydrofem
{

/**
 * \class AFC_LPSLimiter a utility class for computing linearity preserving limiters
 */
class AFC_LPSLimiter : public AFC_Limiter
{
  
  AFC_LPSLimiter() : m_eps(1.0e-16) {}

public:
  
  AFC_LPSLimiter(const std::shared_ptr<Mesh>& mesh,
                 const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                 const std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>>& basis);

  ~AFC_LPSLimiter() {}
  
  /**
   * \brief main routine to build element limiters 
   */
  virtual void
  buildLimiter(const std::shared_ptr<FEVector>& res,
               const std::shared_ptr<const FEVector>& u) const override;

private:
  
  // constant for scaling limiter
  const double m_eps;
  
  //! the mesh for the problem
  std::shared_ptr<Mesh> m_mesh;
  
  //! quadrature for doing integration in limiter calculation
  std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>> m_quadrature;
  
  //! basis functions
  std::shared_ptr<std::vector<std::vector<std::shared_ptr<FEBasis>>>> m_basis;
  
  //! temporaries preallocated for efficiency
  std::shared_ptr<FEVector> m_P;
  std::shared_ptr<FEVector> m_Q;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_AFC_LPSLimiter_HPP__ */
