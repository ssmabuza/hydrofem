// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ComputeError_HPP__
#define __Hydrofem_ComputeError_HPP__

#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_DofMapper.hpp"
#include "Hydrofem_Quadrature.hpp"
#include "Hydrofem_FunctionElement.hpp"
#include "Hydrofem_AnalyticalExpressions.hpp"


namespace hydrofem
{

enum class Norm
{
  L1,
  L2,
  Linf,
  H1
};
  
class ComputeError 
{
public:
  
  static double
  evaluateScalarFieldError(const std::shared_ptr<const std::vector<std::shared_ptr<FEBasis>>>& basis,
                           const std::shared_ptr<const std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                           const std::shared_ptr<const DofMapper>& dofmapper, 
                           const std::shared_ptr<const FEVector>& u_numerical,
                           const std::shared_ptr<const std::function<double(SPoint)>>& u_exact,
                           const std::shared_ptr<const std::function<SPoint(SPoint)>>& grad_u_exact = nullptr,
                           Norm normt = Norm::L1);

  static double
  evaluateScalarFieldError(const std::shared_ptr<const std::vector<std::shared_ptr<FEBasis>>>& basis,
                           const std::shared_ptr<const std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                           const std::shared_ptr<const DofMapper>& dofmapper,
                           const std::shared_ptr<const FEVector>& u_numerical,
                           const std::shared_ptr<ScalarAnalyticalExpression>& u_exact,
                           const std::shared_ptr<VectorAnalyticalExpression>& grad_u_exact = nullptr,
                           Norm normt = Norm::L1);

};


}
// end namespace hydrofem

#endif /** __Hydrofem_ComputeError_HPP__ */
