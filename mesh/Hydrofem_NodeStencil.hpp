// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_NodeStencil_HPP__
#define __Hydrofem_NodeStencil_HPP__

#include "Hydrofem.hpp"

namespace hydrofem
{

struct NodeStencil
  :
  public std::map<int,int>
{
  
  using super = std::map<int,int>;

  NodeStencil() {}
  
  NodeStencil(const super& stencil_) : super(stencil_) {}
  
  NodeStencil& operator=(const super& stencil_)
  { (super)(*this) = stencil_; return *this; }
  
  virtual ~NodeStencil() {}
  
  virtual inline const int& loc_ind(const int elem_ind) const
  { return this->at(elem_ind); }
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_NodeStencil_HPP__ */

