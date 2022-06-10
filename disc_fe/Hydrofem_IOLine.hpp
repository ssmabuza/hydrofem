// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOLine_HPP__
#define __Hydrofem_IOLine_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_HGrad_DofMapper_Line.hpp"

namespace hydrofem
{

class IOLine
  :
  public IOBase
{
public:
  
  IOLine(const std::shared_ptr<HGrad_DofMapper_Line>& dofmapper,
         const std::shared_ptr<FunctionElement<double>>& fe_shape);
  
  virtual ~IOLine() override = default;

  void writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOLine_HPP__ */




