// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_Tri_HPP__
#define __Hydrofem_Quadrature_Tri_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

/**
 * \class Quadrature_Tri a class that builds quadrature on a triangle
 */
class Quadrature_Tri 
  :
  public Quadrature
{
public:
  
  //! Ctor
  Quadrature_Tri(const int order, const std::vector<SPoint>& tri);
  
  //! Dtor
  virtual ~Quadrature_Tri() {}
  
private:
  
  //! a method that builds quadrature on a triangle
  virtual void buildQuadrature(const int order, const std::vector<SPoint>& tri);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_Tri_HPP__ */
