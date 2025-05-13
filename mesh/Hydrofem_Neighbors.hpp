// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_Neighbors_HPP__
#define __Hydrofem_Neighbors_HPP__

#include "Hydrofem.hpp"

namespace hydrofem
{

/** \brief element neighbor index */
struct Neighbors
  :
  public std::vector<int>
{
  
  /** \brief Base */
  using super = std::vector<int>;

  /** \brief Ctor */
  Neighbors() {}
  
  /** \brief Dtor */
  virtual ~Neighbors() {}
  
  /** \brief Copy Ctor */
  Neighbors(const super& neighbors) : super(neighbors) {}
 
  /** \brief Copy assignment */
  Neighbors& operator=(const super& neighbors)
  { (super)(*this) = neighbors; return *this; }
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Neighbors_HPP__ */
