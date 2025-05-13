// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_HGrad_DofMapper_Triangle_HPP__
#define __Hydrofem_HGrad_DofMapper_Triangle_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_Mesh.hpp"
#include "Hydrofem_DofMapper.hpp"

namespace hydrofem
{

class IOTri;
class IOTriXML;

/**
 * \brief A nodal DOF manager that works with triangular elements in 2D
 * 
 * We build the DOF manager from the parallel mesh
 */
class HGrad_DofMapper_Triangle
  : 
  public DofMapper
{

  /**
   * \brief the local DOF layout
   */
  int local(const int i, const int j, const int k) const;

  /**
   * \brief the local DOF layout
   */
  inline int local(const int i, const int j) const
  { return local(i, j, m_p - i - j); }
  
  /**
   * \brief the main DofMapper routine for getting the global index for a given local element index
   */
  int eval_global(const int i_elem, const int i, const int j, const int k) const;
  
  
  /** \brief the main DofMapper routine for getting the global index for a given local element index */
  inline int eval_global(const int i_elem, const int i, const int j) const
  { return eval_global(i_elem, i, j, m_p - i - j); }
  
public:
  
  HGrad_DofMapper_Triangle(const std::shared_ptr<Mesh>& mesh,
                           const int order) : DofMapper(mesh,order)
  {
    m_lndof = int((m_p + 1) * (m_p + 2)) / int(2);
    m_gndof = m_npoints + m_nedges * int(m_p - 1) + m_nelems * (m_lndof - int(3 * m_p));
    
    m_loc_indexes.resize(m_lndof,-1);
    m_glob_indexes.resize(nelements(),std::vector<int>(m_lndof,-1));
    
    int ind = 0;
    for (int i = 0; i <= p(); ++i)
    {
      for (int j = 0; j <= p() - i; ++j)
      {
        m_loc_indexes.at(ind) = local(i,j);
        for (int elem_ind = 0; elem_ind < nelements(); ++elem_ind)
          m_glob_indexes.at(elem_ind).at(ind) = eval_global(elem_ind,i,j);
        ind++;
      }
    }
    
  }
  
  virtual ~HGrad_DofMapper_Triangle() {}
  
  friend class IOTri;
  friend class IOTriXML;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_HGrad_DofMapper_Triangle_HPP__ */
