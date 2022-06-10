// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_TaylorDispersion_impl_HPP__
#define __Hydrofem_Equation_TaylorDispersion_impl_HPP__

namespace hydrofem
{

template <typename ScalarT>
Equation_TaylorDispersion<ScalarT>::
Equation_TaylorDispersion(double V_max, double y_length)
  :
  Equation_Linear<ScalarT>()
{
  this->m_field_names = {"conc"};
  const double m_V_max = V_max;
  const double m_y_length = y_length;
  auto vel = std::make_shared<std::function<SPoint(SPoint)>>(
    [&m_V_max,&m_y_length](SPoint x)->SPoint 
    { auto y(x.y()); return SPoint(4.0*m_V_max*(y/m_length_y)*(1.0-(y/m_length_y)),0.0); }
    );
  this->set_velocity(vel);
}

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_TaylorDispersion_impl_HPP__ */
