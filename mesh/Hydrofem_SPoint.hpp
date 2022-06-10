// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_SPoint_HPP__
#define __Hydrofem_SPoint_HPP__

#include <Eigen/Dense>

#include "Hydrofem_GlobalConstants.hpp"

namespace hydrofem
{

/** \brief A spatial vector */
class SPoint
{
public:

  // raw data type
  using super = Eigen::Array<double,Eigen::Dynamic,1>;  
  // raw constant data type
  using super_const = Eigen::Array<const double,Eigen::Dynamic,1>;
  
private:
  
  super m_data;

public:
  
  SPoint() = default;
  
  explicit SPoint(const int num_dims)
  {
    assert((num_dims>=1) && (num_dims<=3)); 
    m_data = super(num_dims);
    for (super::Index i = 0; i < this->size(); ++i)
      m_data(i) = 0.0;
  }
  
  explicit SPoint(const double x)
  {
    m_data = super(1);
    m_data(0) = x;
  }
  
  SPoint(const double x, const double y)
  {
    m_data = super(2);
    m_data(0) = x;
    m_data(1) = y;
  }
  
  SPoint(const double x, const double y, const double z)
  {
    m_data = super(3);
    m_data(0) = x;
    m_data(1) = y;
    m_data(2) = z;
  }
  
  SPoint(const SPoint& point)
  {
    m_data = super(point.m_data.size());
    for (super::Index i = 0; i < point.size(); ++i)
      m_data(i) = point.m_data(i);
  }

  explicit SPoint(const super& point)
  {
    assert((point.size() >= 1) && (point.size() <= 3));
    m_data = point;
  }
  
  SPoint& operator=(const SPoint& point)
  {
    // do deep copy 
    m_data = super(point.m_data.size());
    for (super::Index i = 0; i < m_data.size(); ++i)
      m_data[i] = point.m_data[i];
    
    return *this;
  }
  
  virtual ~SPoint() = default;
  
  [[nodiscard]] inline super::Index size() const
  { return m_data.size(); }
  
  [[nodiscard]] inline super data() const
  { return m_data; }

  /** \brief l2 norm */
  [[nodiscard]] inline double norm() const
  {
    double res(0.0);
    for (super::Index i = 0; i < m_data.size(); ++i)
      res += m_data(i) * m_data(i);
    return std::sqrt(res);
  }
  
  inline void normalize()
  {
    const double norm_ = norm();
    if (norm_ < length_eps)
    {
      throw std::runtime_error("In SPoint::normalize : Vector norm is 0 or close to 0!");
    }
    for (super::Index i = 0; i < m_data.size(); ++i)
      m_data(i) /= norm_;
  }
  
  inline void reverse()
  {
    for (super::Index i = 0; i < this->size(); ++i)
      m_data(i) *= -1.0;
  }

  double& operator()(const super::Index i)
  {
    return m_data(i);
  }
  
  double operator()(const super::Index i) const
  {
    return m_data(i);
  }  
  
  SPoint operator-() const
  {
    SPoint res(*this);
    res.reverse();
    return res;
  }

  bool operator==(const SPoint& point) const
  {
    assert(size()==point.size());
    double diff(0.0);
    for (super::Index i(0); i < size(); ++i)
      diff += std::fabs(point(i)-(*this)(i));
    return bool(diff < point_eps);
  }

  SPoint& operator+=(const SPoint& point)
  {
    assert(size()==point.size());
    for (super::Index i(0); i < size(); ++i)
      (*this)(i) += point(i);
    return *this;
  }

  SPoint& operator-=(const SPoint& point)
  {
    assert(size()==point.size());
    for (super::Index i(0); i < size(); ++i)
      (*this)(i) -= point(i);
    return *this;
  }

  SPoint& operator*=(const double val)
  {
    for (super::Index i(0); i < size(); ++i)
      (*this)(i) *= val;
    return *this;
  }

  SPoint& operator/=(const double val)
  {
    assert(val!=0.0);
    for (super::Index i(0); i < size(); ++i)
      (*this)(i) /= val;
    return *this;
  }

  friend double operator*(const SPoint& vec0, const SPoint& vec1)
  {
    assert(vec0.size()==vec1.size());
    double res(0.0);
    for (super::Index i(0); i < vec0.size(); ++i)
      res += vec0(i)*vec1(i);
    return res;
  }

  friend SPoint operator*(const SPoint& vec, const double alpha)
  {
    SPoint res(vec);
    for (super::Index i(0); i < vec.size(); ++i)
      res(i) *= alpha;
    return res;
  }
  
  friend SPoint operator*(const double alpha, const SPoint& vec)
  {
    SPoint res(vec);
    for (super::Index i(0); i < vec.size(); ++i)
      res(i) *= alpha;
    return res;
  }

  friend SPoint operator/(const SPoint& vec, const double alpha)
  {
    assert(alpha!=0.0);
    SPoint res(vec);
    for (super::Index i(0); i < vec.size(); ++i)
      res(i) /= alpha;
    return res;
  }

  friend SPoint operator+(const SPoint& vec0, const SPoint& vec1)
  {
    SPoint res(vec0);
    for (super::Index i(0); i < vec0.size(); ++i)
      res(i) += vec1(i);
    return res;
  }

  friend SPoint operator-(const SPoint& vec0, const SPoint& vec1)
  {
    SPoint res(int(vec0.size()));
    for (super::Index i(0); i < vec0.size(); ++i)
      res(i) = vec0(i) - vec1(i);
    return res;
  }
  
  [[nodiscard]] inline double x() const { assert(size()>=1); return this->operator()(0); }
  
  [[nodiscard]] inline double y() const { assert(size()>=2); return this->operator()(1); }
  
  [[nodiscard]] inline double z() const { assert(size()==3); return this->operator()(2); }
  
  inline double& x() { assert(size()>=1); return this->operator()(0); }
  
  inline double& y() { assert(size()>=2); return this->operator()(1); }
  
  inline double& z() { assert(size()==3); return this->operator()(2); }
  
  inline void print(std::ostream& output) const
  {
    if (size()==1)
      output << "( " << operator()(0) << " )" << std::endl;
    else if (size()==2)
      output << "( " << operator()(0) << ", " << operator()(1) << " )" << std::endl;
    else if (size()==3)
      output << "( " << operator()(0) << ", " << operator()(1) << ", " << operator()(2) << " )" << std::endl;
  }
  
};

inline std::ostream& operator<<(std::ostream& output, const SPoint& point)
{
  point.print(output);
  return output;
}

}
// end namespace hydrofem

#endif /** __Hydrofem_SPoint_HPP__ */
