// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ElementShapeTools_HPP__
#define __Hydrofem_ElementShapeTools_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"

namespace hydrofem
{

/** \brief routine to compute the area of a triangle in 2D */
inline double computeTriangleArea(const std::vector<SPoint>& triangle)
{
  assert(int(triangle.size()) == 3);
  const SPoint v1 = triangle[1] - triangle[0];
  const SPoint v2 = triangle[2] - triangle[0];
  return 0.5*std::fabs(v1.x()*v2.y() - v1.y()*v2.x());
}

/** \brief routine to compute the area of a triangle in 2D */
inline double computeTriangleArea(const SPoint& point1, const SPoint& point2, const SPoint& point3)
{
  const SPoint v1 = point2 - point1;
  const SPoint v2 = point3 - point1;
  return 0.5*std::fabs( v1.x()*v2.y() - v1.y()*v2.x() );
}

/** \brief convert Barycentric coordinates to Cartesian in 2D simplexes */
inline void convertBarycentricToCartesian(const SPoint& point0,
                                          const SPoint& point1,
                                          const SPoint& point2,
                                          double l0,
                                          double l1,
                                          double l2,
                                          double& x,
                                          double& y)
{
  x = l0*point0.x() + l1*point1.x() + l2*point2.x();
  y = l0*point0.y() + l1*point1.y() + l2*point2.y();
}

/** \brief convert Cartesian coordinates to Barycentric in 2D simplexes */
inline void convertCartesianToBarycentric(double& l0,
                                          double& l1,
                                          double& l2,
                                          const SPoint& point,
                                          const SPoint& tri0,
                                          const SPoint& tri1,
                                          const SPoint& tri2)
{
  const SPoint v0 = tri1 - tri0;
  const SPoint v1 = tri2 - tri0;
  const SPoint v2 = point - tri0;

  const double d00 = v0 * v0;
  const double d01 = v0 * v1;
  const double d11 = v1 * v1;
  const double d20 = v2 * v0;
  const double d21 = v2 * v1;
  const double denom = d00 * d11 - d01 * d01;
  l0 = (d11 * d20 - d01 * d21) / denom;
  l1 = (d00 * d21 - d01 * d20) / denom;
  l2 = double(1.0) - l0 - l1;  
}

/** \brief scalar product between 2 points */
inline double scalarProduct(const SPoint& point1, const SPoint& point2)
{
  return point1*point2;
}

/** \brief distance between 2 points */
inline double evalDistance(const SPoint& point1, const SPoint& point2)
{
  return (point1-point2).norm();
}

inline double tetrahedronVolume(const SPoint& p0, const SPoint& p1, const SPoint& p2, const SPoint& p3)
{
  using vectype = Eigen::Vector3d;    
  // calculate the volume of the tetrahedron using the formula:
  // V = (1/6) * |(p1 - p0) . ((p2 - p0) x (p3 - p0))|
  const vectype vec1 = (p1 - p0).data();
  const vectype vec2 = (p2 - p0).data();
  const vectype vec3 = (p3 - p0).data();
  const auto cross_prod = vec2.cross(vec3);
  const double dot_prod = vec1.dot(cross_prod);
  const double volume = std::fabs(dot_prod) / 6.0;
  return volume;
}

}
// end namespace hydrofem

#endif /** __Hydrofem_ElementShapeTools_HPP__ */
