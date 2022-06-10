// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_Quad_HPP__
#define __Hydrofem_Quadrature_Quad_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

/**
 * \class Quadrature_Quad a class that builds quadrature on a quadrilateral
 */
class Quadrature_Quad
  :
  public Quadrature
{
  
public:
  
  //! Ctor
  Quadrature_Quad(const int order, const std::vector<SPoint>& quad);
  
  //! Dtor
  ~Quadrature_Quad() {}

private:
  
  //! a method that builds quadrature on a quadrilateral
  virtual void buildQuadrature(const int order, const std::vector<SPoint>& quad);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_Quad_HPP__ */
