// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOTri_HPP__
#define __Hydrofem_IOTri_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_HGrad_DofMapper_Triangle.hpp"

namespace hydrofem
{

/** \brief IO for a 2D simplex mesh into VTK */
class IOTri
  : 
  public IOBase
{
public:

  IOTri(const std::shared_ptr<HGrad_DofMapper_Triangle>& dofmapper,
        const std::shared_ptr<FunctionElement<double>>& fe_shape)
    :
    IOBase(dofmapper,fe_shape)
  { }
  
  ~IOTri() override = default;

  void writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const override;

private:

  void printVTK(const std::string& filename,
                const std::string& directory) const;

};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOTri_HPP__ */
