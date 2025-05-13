// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_Edge2D_HPP__
#define __Hydrofem_Quadrature_Edge2D_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

class Quadrature_Edge2D
  :
  public Quadrature
{
  
public:
  
  //! Ctor
  Quadrature_Edge2D(const int order, const std::vector<SPoint>& edge);

  //! Dtor  
  virtual ~Quadrature_Edge2D() {}
  
private:
  
  virtual void buildQuadrature(const int order, const std::vector<SPoint>& edge);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_Edge2D_HPP__ */
