// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_GlobalConstants_HPP__
#define __Hydrofem_GlobalConstants_HPP__

#include <climits>

namespace hydrofem 
{

const double point_eps = 1.0e-10; /// we assume that two points are equal if the distance between them is less than point_eps
const double volume_eps = 1.0e-16;
const double area_eps = 1.0e-15;
const double length_eps = 1.0e-11;
const unsigned long int wrongInteger = ULONG_MAX;

}
// end namespace hydrofem

#endif /** __Hydrofem_GlobalConstants_HPP__ */
