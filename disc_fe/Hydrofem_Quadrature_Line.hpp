// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_Line_HPP__
#define __Hydrofem_Quadrature_Line_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_Quadrature.hpp"

namespace hydrofem
{

class Quadrature_Line
  :
  public Quadrature
{
  
public:
  
  //! Ctor
  Quadrature_Line(const int order, const std::vector<SPoint>& line);
  
  //! Dtor
  virtual ~Quadrature_Line() {}
  
private:

  void buildQuadrature(const int order, const std::vector<SPoint>& line);
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_Line_HPP__ */
