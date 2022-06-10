// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_ReferenceQuadrilateral_HPP__
#define __Hydrofem_ReferenceQuadrilateral_HPP__

#include "Hydrofem.hpp"
#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_LocalArray.hpp"

namespace hydrofem
{
    
/**
 * \brief This code works for mappings between the Bezier cell [0,1]x[0,1] and a general element.
 *
 *  NOTE: This is different from the implementation in FEQuadrature_Quadrilateral.hpp which works 
 *        with [-1,1]^2. Both codes make use of ElementShapeTools for area calculation. All quads
 *        are oriented in an anticlockwise direction.
 *
 */    
    

//! \brief computes the area of a quad shape
double quadArea(const SPoint& point1,
                const SPoint& point2,
                const SPoint& point3,
                const SPoint& point4);
                          
//! \brief computes the Jacobian determinant in a given quad
double getJacobianDeterminant(const SPoint& ref_point,
                              const SPoint& point1,
                              const SPoint& point2,
                              const SPoint& point3,
                              const SPoint& point4);

//! \brief maps physical coordinates to reference cell
SPoint physicalQuadToReference(const SPoint& quad1,
                               const SPoint& quad2,
                               const SPoint& quad3,
                               const SPoint& quad4,
                               const SPoint& physicalPoint);

//! \brief maps reference coordinates to physical cell
SPoint referenceQuadToPhysical(const SPoint& quad1,
                               const SPoint& quad2,
                               const SPoint& quad3,
                               const SPoint& quad4,
                               const SPoint& referencePoint);

//! \brief computes d(x,y)/d(xi,eta) = [[dx_dxi dy_dxi];[dx_deta dy_deta]]
LMAT_<double> getJacobian(const SPoint& ref_point,
                          const SPoint& point1,
                          const SPoint& point2,
                          const SPoint& point3,
                          const SPoint& point4);

//! \brief computes d(xi,eta)/d(x,y) = [[dxi_x deta_dx];[dxi_dy deta_dy]]
LMAT_<double> getJacobianInverse(const SPoint& quad1,
                                 const SPoint& quad2,
                                 const SPoint& quad3,
                                 const SPoint& quad4,
                                 const double x, 
                                 const double y);

//! \brief computes eta_x at a physical coordinate
double deta_dx(const SPoint& quad1,
               const SPoint& quad2,
               const SPoint& quad3,
               const SPoint& quad4,
               const double x, 
               const double y);

//! \brief computes eta_y at a physical coordinate
double deta_dy(const SPoint& quad1,
               const SPoint& quad2,
               const SPoint& quad3,
               const SPoint& quad4,
               const double x, 
               const double y);

//! \brief computes xi_x at a physical coordinate
double dxi_dx(const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const double x, 
              const double y);

//! \brief computes xi_y at a physical coordinate
double dxi_dy(const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const double x, 
              const double y);

}
// end namespace hydrofem

#endif /** __Hydrofem_ReferenceQuadrilateral_HPP__ */
