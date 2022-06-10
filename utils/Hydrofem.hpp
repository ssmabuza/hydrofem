// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_HPP__
#define __Hydrofem_HPP__

// define this macro for windows
#define _USE_MATH_DEFINES

// include all the essential headers

#include <map>
#include <cmath>
#include <array>
#include <ctime>
#include <cstdio>
#include <string>
#include <cctype>
#include <vector>
#include <memory>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <utility>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unordered_map>
#include <boost/program_options.hpp>

// project namespace 
namespace hydrofem {}
// namespace alias to avoid typo error
namespace Hydrofem = hydrofem;

#endif /** __Hydrofem_HPP__ */
