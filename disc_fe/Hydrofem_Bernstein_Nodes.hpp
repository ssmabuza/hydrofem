// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_BernsteinNodes_HPP__
#define __Hydrofem_BernsteinNodes_HPP__


#include "Hydrofem_DofMapper.hpp"

namespace hydrofem
{

class Bernstein_Nodes 
{
public:  
  
  explicit Bernstein_Nodes(const std::shared_ptr<DofMapper>& dofmapper)
  {
    m_dofmapper = dofmapper;
  }
  
  std::vector<SPoint>
  getNodes(const int elem_ind) const;
  
private:
  
  std::shared_ptr<DofMapper> m_dofmapper;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_BernsteinNodes_HPP__ */
