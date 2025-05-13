// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_Point_HPP__
#define __Hydrofem_Quadrature_Point_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

/**
 * For 1D integration
 */
class Quadrature_Point
  :
  public Quadrature
{
  
public:
  
  //! Ctor
  Quadrature_Point(const SPoint& point);
  
  //! Dtor
  virtual ~Quadrature_Point() {}
  
private:

  void buildQuadrature(const SPoint& point);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_Point_HPP__ */
