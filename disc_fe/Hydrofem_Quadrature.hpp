// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Quadrature_HPP__
#define __Hydrofem_Quadrature_HPP__

#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

class Quadrature
{
  
public:
  
  //! Default Ctor
  Quadrature() {}
  
  //! Dtor
  virtual ~Quadrature() {}
  
  //! \brief get the quadrature points
  virtual inline const std::vector<SPoint>& get_q_points() const
  { return m_q_points; }
  
  //! \brief get the quadrature weights
  virtual inline const std::vector<double>& get_q_weights() const
  { return m_q_weights; }
  
  //! \brief get the quadrature size (number of weights or points)
  virtual inline int size() const
  { return m_q_points.size(); }
  
protected:  
  
  /**
   * \brief A function that builds the quadrature rule provided the \p order and the element \p shape
   */
//   virtual void buildQuadrature(const int /*order*/, const std::vector<SPoint>& /*shape*/) {}
  
  //! \brief get a view of the points for modification
  inline virtual std::vector<SPoint>& get_q_pointsView()
  { return m_q_points; }
  
  //! \brief get a view of the weights for modification
  inline virtual std::vector<double>& get_q_weightsView()
  { return m_q_weights; }
  
private:
  
  //! quad points
  std::vector<SPoint> m_q_points;
  
  //! quad weights
  std::vector<double> m_q_weights;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Quadrature_HPP__ */
