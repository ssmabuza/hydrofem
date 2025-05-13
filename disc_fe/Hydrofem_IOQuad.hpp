// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOQuad_HPP__
#define __Hydrofem_IOQuad_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_HGrad_DofMapper_Quadrilateral.hpp"

namespace hydrofem
{

class IOQuad 
  :
  public IOBase
{
public:
  
  IOQuad(const std::shared_ptr<HGrad_DofMapper_Quadrilateral>& dofmapper,
         const std::shared_ptr<FunctionElement<double>>& fe_shape);
  
  ~IOQuad() override = default;

  void writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const override;
  
private:

  // internally built output points to speed up writing to file
  std::vector<std::vector<SPoint>> m_out_points;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOQuad_HPP__ */
