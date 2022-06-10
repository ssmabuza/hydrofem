// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_ModelEvaluator.hpp"

namespace hydrofem 
{

void
ModelEvaluator::evalModel(const ModelEvaluator::InArgs& in_args,
                          const ModelEvaluator::OutArgs& out_args) const
{
  if (in_args.supports_delta_t() && in_args.supports_time() && in_args.supports_u_dot())
    m_assembler->buildResidualAndJacobian(in_args.get_u(),in_args.get_u_dot(),out_args.get_f(),out_args.get_W(),
                                          in_args.get_time(),in_args.get_delta_t(),in_args.get_beta());
  else
    m_assembler->buildResidualAndJacobian(in_args.get_u(),out_args.get_f(),out_args.get_W(),in_args.get_beta());
}

}
// end namespace hydrofem

