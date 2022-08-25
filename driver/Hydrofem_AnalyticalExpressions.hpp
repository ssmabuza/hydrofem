// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_AnalyticalExpressions_HPP__
#define __Hydrofem_AnalyticalExpressions_HPP__

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

/**
 * \brief Base class for scalar analytical expression for
 *        boundary condition values and exact solutions.
 */
class ScalarAnalyticalExpression
{
public:
  
  using ScalarT = double;
  
  //! \brief Ctor
  ScalarAnalyticalExpression () = default;
  
  //! \brief Dtor
  virtual ~ScalarAnalyticalExpression () = default;
  
  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& /*x*/)
  { assert(!m_is_transient); return std::numeric_limits<double>::quiet_NaN(); }

  //! \brief can be used for exact solution/ boundary values
  virtual double evaluate(const SPoint& /*x*/, const double /*t*/)
  { assert(m_is_transient); return std::numeric_limits<double>::quiet_NaN(); }
  
  //! \brief make this a functor
  virtual double operator()(const SPoint& /*x*/)
  { assert(!m_is_transient); return std::numeric_limits<double>::quiet_NaN(); }
  
  //! \brief make this a functor
  virtual double operator()(const SPoint& /*x*/, const double /*t*/)
  { assert(m_is_transient); return std::numeric_limits<double>::quiet_NaN(); }
  
protected:

  // true if time dependent
  bool m_is_transient = false;
  
};

/**
 * \brief Base class for vector analytical expression for
 *        boundary condition values and exact solutions.
 */
class VectorAnalyticalExpression
{
public:

  //! \brief Ctor
  VectorAnalyticalExpression () = default;
  
  //! \brief Dtor
  virtual ~VectorAnalyticalExpression () = default;
  
  //! \brief can be used for exact solution/ boundary values
  virtual LVEC_<double> evaluate(const SPoint& /*x*/)
  {
    assert(!m_is_transient);
    auto res = createKArray<LVEC_<double>>(1);
    res(0) = std::numeric_limits<double>::quiet_NaN();
    return res;
  }
  
  //! \brief can be used for exact solution/ boundary values
  virtual LVEC_<double> evaluate(const SPoint& /*x*/, const double /*t*/)
  {
    assert(m_is_transient);
    auto res = createKArray<LVEC_<double>>(1);
    res(0) = std::numeric_limits<double>::quiet_NaN();
    return res;
  }
  
  //! \brief make this a functor
  virtual LVEC_<double> operator()(const SPoint& /*x*/)
  {
    assert(!m_is_transient);
    auto res = createKArray<LVEC_<double>>(1);
    res(0) = std::numeric_limits<double>::quiet_NaN();
    return res;
  }
  
  //! \brief make this a functor
  virtual LVEC_<double> operator()(const SPoint& /*x*/, const double /*t*/)
  {
    assert(m_is_transient);
    auto res = createKArray<LVEC_<double>>(1);
    res(0) = std::numeric_limits<double>::quiet_NaN();
    return res;
  }
  
protected:

  // true if time dependent
  bool m_is_transient = false;
  
};

//! \brief place holder for use with hyperbolic impl
using AnalyticalExpression = VectorAnalyticalExpression;

}
// end namespace hydrofem

#endif /** __Hydrofem_AnalyticalExpressions_HPP__ */
