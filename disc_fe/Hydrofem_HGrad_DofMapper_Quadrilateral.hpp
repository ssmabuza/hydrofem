// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_HGrad_DofMapper_Quadrilateral_HPP__
#define __Hydrofem_HGrad_DofMapper_Quadrilateral_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_DofMapper.hpp"

namespace hydrofem
{

class IOQuad;
 
class HGrad_DofMapper_Quadrilateral
  :
  public DofMapper
{
  
  //! get the global index
  int eval_global(const int i_elem, const int i, const int j) const;

  //! get the local index
  int local(const int i, const int j) const;
  
public:
  
  //! Ctor
  HGrad_DofMapper_Quadrilateral(const std::shared_ptr<Mesh>& mesh,
                                const int order) : DofMapper(mesh,order)
  {
    m_lndof = (m_p+1)*(m_p+1);
    m_gndof = m_npoints+m_nedges*(m_p-1)+m_nelems*(m_p-1)*(m_p-1);

    m_loc_indexes.resize(m_lndof,-1);
    m_glob_indexes.resize(nelements(),std::vector<int>(m_lndof,-1));
    
    int ind = 0;
    for (int i = 0; i <= p(); ++i)
    {
      for (int j = 0; j <= p(); ++j)
      {
        m_loc_indexes.at(ind) = local(i,j);
        for (int elem_ind = 0; elem_ind < nelements(); ++elem_ind)
        {
          m_glob_indexes.at(elem_ind).at(ind) = eval_global(elem_ind,i,j);
        }
        ind++;
      }
    }
  }
  
  virtual ~HGrad_DofMapper_Quadrilateral() {}
  
  friend class IOQuad;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_HGrad_DofMapper_Quadrilateral_HPP__ */
