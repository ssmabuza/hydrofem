// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_String_Utilities_HPP__
#define __Hydrofem_String_Utilities_HPP__

#include <vector>
#include <string>
#include <iomanip>

namespace hydrofem
{

//! Removes whitespace at beginning and end of string
void trim(std::string& str);

//! Tokenize a string, put tokens in a vector
void StringTokenizer(std::vector<std::string>& tokens,
                     const std::string& str,
                     const std::string delimiter = ",",
                     bool trim = false);

//! Turn a vector of tokens into a vector of doubles
void TokensToDoubles(std::vector<double>& values, const std::vector<std::string>& tokens);

//! Turn a vector of tokens into a vector of ints
void TokensToInts(std::vector<int>& values, const std::vector<std::string>& tokens);


inline
std::string intToStrSixDigits(int num)
{
  std::stringstream ss;
  ss << std::setw(6) << std::setfill('0') << num;
  return ss.str();
}

}
// end namespace hydrofem

#endif /** __Hydrofem_String_Utilities_HPP__ */
