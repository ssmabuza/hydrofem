// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_HGrad_DofMapper_Line_HPP__
#define __Hydrofem_HGrad_DofMapper_Line_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_DofMapper.hpp"

namespace hydrofem
{

class IOLine;
  
class HGrad_DofMapper_Line
  :
  public DofMapper
{

  //! get the global index
  int eval_global(const int i_elem, const int k) const;

  //! get the local index
  int local(const int k) const;

public:
  
  HGrad_DofMapper_Line(const std::shared_ptr<Mesh>& mesh,
                       const int order) : DofMapper(mesh,order)
  {
    m_lndof = m_p+1;
    m_gndof = m_npoints+m_nelems*(m_p-1);
    
    m_loc_indexes.resize(m_lndof,-1);
    m_glob_indexes.resize(nelements(),std::vector<int>(m_lndof,-1));
    
    int ind = 0;
    for (int i = 0; i <= p(); ++i)
    {
      m_loc_indexes.at(ind) = local(i);
      for (int elem_ind = 0; elem_ind < nelements(); ++elem_ind)
        m_glob_indexes.at(elem_ind).at(ind) = eval_global(elem_ind,i);
      ind++;
    }
  }
  
  virtual ~HGrad_DofMapper_Line() {}
  
  friend class IOLine;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Bernstein_DofMapper_Line_HPP__ */
