// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_Problem_ReactiveFlow.hpp"

namespace hydrofem
{

void Problem_ReactiveFlow::init()
{
  // get the diffusion coefficient
  m_diff = m_pl.get<double>("Diffusion Coefficient");
  // build the analytical velocity profile
  const double H = m_pl.get<double>("H");
  const double V_max = m_pl.get<double>("V_max");
  std::function<SPoint(SPoint)> vel_fnc = [&H,&V_max](SPoint pt)->SPoint { return SPoint(V_max*pt.y()*(H - pt.y())/(H*H/4.0),0.0); };
  m_fluid_velocity = Teuchos::rcp(&vel_fnc);
  
  // build the adsorption isotherm
  std::string isotherm = m_pl.get<std::string>("Isotherm");
  // get the adsorption rate
  const double k_a = m_pl.get<double>("Adsorption Rate");
  // get the desorption rate
  k_d = m_pl.get<double>("Desorption Rate");
  const double k_d_2 = k_d;
  if (isotherm == "Linear")
  {
    std::function<double(double)> ads_fnc = [&k_a](double c)->double { return k_a*c; };
    m_ads = Teuchos::rcp(&ads_fnc);
  } 
  else if (isotherm == "Langmuir")
  {
    std::function<double(double)> ads_fnc = [&k_d_2,&k_a](double c)->double { return k_a*c/(1.0 + k_d_2*c); };
    m_ads = Teuchos::rcp(&ads_fnc);
  }
  
  std::function<double(SPoint,double)> inflow_fnc = [](SPoint,double)->double { return 1.0; };
  m_inflow = Teuchos::rcp(&inflow_fnc);
}

}
// end namespace hydrofem
