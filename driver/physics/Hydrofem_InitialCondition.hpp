// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_InitialCondition_HPP__
#define __Hydrofem_InitialCondition_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{

/**
 * \class InitialCondition base class for scalar initial conditions
 */
class ScalarInitialCondition
{
public:

  using LVec = LVEC_<double>;

  //! \brief Ctor
  ScalarInitialCondition() = default;

  //! \brief Dtor
  virtual ~ScalarInitialCondition() = default;

  //! \brief evaluate value
  [[nodiscard]] virtual double evaluate(const SPoint&) const = 0;

};

/**
 * \class InitialCondition base class for vector initial conditions
 */
class VectorInitialCondition
{
public:
  
  using LVec = LVEC_<double>;
  
  //! \brief Ctor
  VectorInitialCondition() = default;
  
  //! \brief Dtor
  virtual ~VectorInitialCondition() = default;
  
  //! \brief evaluate value
  [[nodiscard]] virtual LVec evaluate(const SPoint&) const = 0;
  
};

/**
 * 
 * Some simple initial conditions
 * 
 */
class ZeroScalarInitialCondition
  :
  public ScalarInitialCondition
{
public:
  
  ZeroScalarInitialCondition() {}

  virtual ~ZeroScalarInitialCondition() {}

  double evaluate(const SPoint&) const override
  {
    return 0.;
  }
  
};

class ZeroVectorInitialCondition
  :
  public VectorInitialCondition
{
public:
  
  using LVec = VectorInitialCondition::LVec;
  
  explicit ZeroVectorInitialCondition(int n) : m_neq(n) {}

  virtual ~ZeroVectorInitialCondition() {}

  LVec evaluate(const SPoint&) const override
  {
    auto res = createKArray<LVec>(m_neq); zeroOutArray(res);
    return res;
  }
  
private:

  int m_neq;
  
};

class ConstantScalarInitialCondition
  :
  public ScalarInitialCondition
{
public:
  
  explicit ConstantScalarInitialCondition(double val) : m_val(val) {}

  virtual ~ConstantScalarInitialCondition() {}

  double evaluate(const SPoint&) const override
  {
    return m_val;
  }

private:

  double m_val;
  
};

class ConstantVectorInitialCondition
  :
  public VectorInitialCondition
{
public:
  
  using LVec = VectorInitialCondition::LVec;
  
  ConstantVectorInitialCondition(const std::vector<double>& val) : m_val(val) {}

  ~ConstantVectorInitialCondition() {}

  LVec evaluate(const SPoint&) const override
  {
    auto ret = createKArray<LVec>(m_val.size());
    for (std::size_t i = 0; i < m_val.size(); ++i) ret(i) = m_val.at(i);
    return ret;
  }
  
private:

  std::vector<double> m_val;
  
};

class FunctionalScalarInitialCondition
  :
  public ScalarInitialCondition
{
public:
  
  explicit
  FunctionalScalarInitialCondition(const std::function<double(SPoint)>& val) : m_val(val) {}

  ~FunctionalScalarInitialCondition() {}

  double evaluate(const SPoint& x) const override
  {
    return m_val(x);
  }

private:

  const std::function<double(SPoint)>& m_val;
  
};

class FunctionalVectorInitialCondition
  :
  public VectorInitialCondition
{
public:
  
  using LVec = VectorInitialCondition::LVec;
  
  explicit
  FunctionalVectorInitialCondition(const std::function<LVec(SPoint)>& val) : m_val(val) {}

  ~FunctionalVectorInitialCondition() {}

  LVec evaluate(const SPoint& x) const override
  {
    return m_val(x);
  }

private:

  const std::function<LVec(SPoint)>& m_val;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_InitialCondition_HPP__ */
