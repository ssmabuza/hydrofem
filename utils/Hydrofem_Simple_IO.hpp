// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Simple_IO_HPP__
#define __Hydrofem_Simple_IO_HPP__

#include "Hydrofem.hpp"

namespace hydrofem
{

template <typename VecT>
struct Simple_IO
{
  
  static inline
  void writeData(const bool append_flag,
                 const int len_u,
                 const std::string u_filename,
                 const VecT& u);

};

template <typename VecT>
void Simple_IO<VecT>::
writeData(const bool append_flag,
          const int len_u,
          const std::string u_filename,
          const VecT& u)
{
  std::shared_ptr<std::ofstream> outw;
  if (append_flag)
    outw = std::make_shared<std::ofstream>(u_filename, std::ios::app);
  else
    outw = std::make_shared<std::ofstream>(u_filename);
  
  if (!outw->good())
  {
    std::stringstream ss;
    ss << "Error in void Simple_IO::writeData" << std::endl
       << "Can not open output file " << u_filename << std::endl;
    throw std::runtime_error(ss.str());
  }
  
  outw->precision(18);
  for (int i = 0; i < len_u; ++i)
  {
    if (std::fabs(u[i]) < 1.0e-18)
      *outw << 0.0 << " ";
    else
      *outw << u[i] << " ";
  }
  *outw << std::endl;
}

}
// end namespace hydrofem

#endif /** __Hydrofem_Simple_IO_HPP__ */
