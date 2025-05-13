// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_IOTriXML_HPP__
#define __Hydrofem_IOTriXML_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_IOBase.hpp"
#include "Hydrofem_EigenFEVecMat.hpp"
#include "Hydrofem_HGrad_DofMapper_Triangle.hpp"

namespace hydrofem
{

/** \brief IO for a 2D simplex mesh into VTK */
class IOTriXML
  :
    public IOBase
{
public:
  
  [[maybe_unused]] IOTriXML(const std::shared_ptr<HGrad_DofMapper_Triangle>& dofmapper,
                            const std::shared_ptr<FunctionElement<double>>& fe_shape,
                            bool fsi = false, bool steady = false)
    :
    IOBase(dofmapper,fe_shape)
  {
    // creating output points and values
    m_out_points.resize(2,std::make_shared<FEVector>(dofmapper->global_ndof()));
    m_out_vals = std::make_shared<FEVector>(dofmapper->global_ndof());
    m_fsi = fsi;
    m_steady = steady;
  }
  
  ~IOTriXML() override = default;
  
  void writeSolutionTofile(const std::shared_ptr<FEVector>& sol, const double time) const override;

private:
  
  void printVTK(const std::string& filename,
                const std::string& directory) const;
  
  // output points are changed only once except for FSI simulations
  mutable std::vector<std::shared_ptr<FEVector>> m_out_points;
  // output values changed every output step
  mutable std::shared_ptr<FEVector> m_out_vals;
  // for an fsi simulation the mesh changes and we update out points
  bool m_fsi;
  // steady problems
  bool m_steady;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_IOTriXML_HPP__ */
