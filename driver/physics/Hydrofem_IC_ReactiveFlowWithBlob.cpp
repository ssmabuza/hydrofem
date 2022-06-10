// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IC_ReactiveFlowWithBlob.hpp"

namespace hydrofem
{

IC_ReactiveFlowWithBlob::LVec
IC_ReactiveFlowWithBlob::evaluate(const SPoint &xP) const
{
  LVec u = createKArray<LVec>(2);
  const double x = xP.x();
  const double y = xP.y();
  u(0) = 0.0;
  u(1) = 0.0;
  if ((x >= 0.2 && x <= 0.8) && (y >= 0.2 && y <= 0.8))
    u(0) = 1.0;
  return u;
}

}
// end namespace hydrofem
