// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Equation_Euler_impl_HPP__
#define __Hydrofem_Equation_Euler_impl_HPP__

#include "Hydrofem_Equation_Euler.hpp"

namespace hydrofem
{

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LVec Equation_Euler<ScalarT>::
flux(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  return internal::flux(U,normal,m_gamma,dim());
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat Equation_Euler<ScalarT>::
flux_u(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  LMat y = createKArray<LMat>(ns(),ns());

  if (dim()==1)
  {

    const double& nx = normal.x();
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT u = rhou/rho;
    const ScalarT H_val = internal::H(U,m_gamma,dim());

    y(0,0) = 0 * nx;
    y(1,0) = 0.5*(m_gamma-3.0)*u*u * nx;
    y(2,0) = u*(0.5*(m_gamma-1.0)*u*u-H_val) * nx;

    y(0,1) = 1.0 * nx;
    y(1,1) = ((3.0-m_gamma)*u) * nx;
    y(2,1) = (H_val-(m_gamma-1.0)*u*u) * nx;

    y(0,2) = 0.0 * nx;
    y(1,2) = (m_gamma-1.0) * nx;
    y(2,2) = (m_gamma*u) * nx;

  } else if (dim()==2) {

    const double& nx = normal.x();
    const double& ny = normal.y();

    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);

    const ScalarT u = rhou/rho;
    const ScalarT v = rhov/rho;
    const ScalarT vn = u*nx+v*ny;
    const ScalarT a = internal::c(U,m_gamma,this->ndims);
    const ScalarT h = a*a/(m_gamma-1.0);
    const ScalarT e_k = 0.5*(u*u+v*v);
    const ScalarT h0 = h + e_k;

    y(0,0) = 0.0;
    y(0,1) = nx;
    y(0,2) = ny;
    y(0,3) = 0.0;

    y(1,0) = (m_gamma-1.0)*e_k*nx-u*vn;
    y(1,1) = vn - (m_gamma-2.0)*u*nx;
    y(1,2) = u*ny - (m_gamma-1.)*v*nx;
    y(1,3) = (m_gamma-1.0)*nx;

    y(2,0) = (m_gamma-1.0)*e_k*ny-v*vn;
    y(2,1) = v*nx-(m_gamma-1.0)*u*ny;
    y(2,2) = vn-(m_gamma-2.0)*v*ny;
    y(2,3) = (m_gamma-1.0)*ny;

    y(3,0) = ( (m_gamma-1.0)*e_k - h0 )*vn;
    y(3,1) = h0*nx - (m_gamma-1.0)*u*vn;
    y(3,2) = h0*ny - (m_gamma-1.0)*v*vn;
    y(3,3) = m_gamma*vn;

  } else if (dim()==3) {

    const double& nx = normal.x();
    const double& ny = normal.y();
    const double& nz = normal.z();

    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhow = U(3);

    const ScalarT u = rhou/rho;
    const ScalarT v = rhov/rho;
    const ScalarT w = rhow/rho;
    const ScalarT vn = u*nx+v*ny+w*nz;
    const ScalarT a = internal::c(U,m_gamma,dim());
    const ScalarT h = a*a/(m_gamma-1.0);
    const ScalarT e_k = 0.5*(u*u+v*v+w*w);
    const ScalarT h0 = h + e_k;

    y(0,0) = 0.0;
    y(0,1) = nx;
    y(0,2) = ny;
    y(0,3) = 0.0;

    y(1,0) = (m_gamma-1.0)*e_k*nx-u*vn;
    y(1,1) = vn - (m_gamma-2.0)*u*nx;
    y(1,2) = u*ny - (m_gamma-1.)*v*nx;
    y(1,3) = (m_gamma-1.0)*nx;

    y(2,0) = (m_gamma-1.0)*e_k*ny-v*vn;
    y(2,1) = v*nx-(m_gamma-1.0)*u*ny;
    y(2,2) = vn-(m_gamma-2.0)*v*ny;
    y(2,3) = (m_gamma-1.0)*ny;

    y(3,0) = ( (m_gamma-1.0)*e_k - h0 )*vn;
    y(3,1) = h0*nx - (m_gamma-1.0)*u*vn;
    y(3,2) = h0*ny - (m_gamma-1.0)*v*vn;
    y(3,3) = m_gamma*vn;

  }

  return y;

}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LVec Equation_Euler<ScalarT>::
lambda(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  return internal::lambda(U,normal,m_gamma,dim());
}

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::
lambda_max(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  return internal::lambda_max(U,normal,m_gamma,dim());
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat Equation_Euler<ScalarT>::
eigenvector(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  return internal::eigenvector(U,normal,m_gamma,dim());
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat Equation_Euler<ScalarT>::
inverse_eigenvector(const SPoint& /*x*/, const double /*t*/, const LVec& U, const SPoint& normal) const
{
  return internal::inverse_eigenvector(U,normal,m_gamma,dim());
}


/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::
internal::pressure(const LVec& U, const double gamma, const int nDims)
{
  ScalarT p(0.0);

  if (nDims==1)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhoE = U(3);
    const ScalarT halfoverrho = 0.5/rho;
    p = (gamma-1.0)*(rhoE-halfoverrho*(rhou*rhou));
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhoE = U(3);
    const ScalarT halfoverrho = 0.5/rho;
    p = (gamma-1.0)*(rhoE-halfoverrho*(rhou*rhou+rhov*rhov));
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhow = U(3);
    const ScalarT& rhoE = U(4);
    const ScalarT halfoverrho = 0.5/rho;
    p = (gamma-1.0)*(rhoE-halfoverrho*(rhou*rhou+rhov*rhov+rhow*rhow));
  }

  return p;
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LVec
Equation_Euler<ScalarT>::
internal::flux(const LVec& U, const SPoint& normal, const double gamma, const int nDims)
{
  assert(normal.size()==Teuchos::as<std::size_t>(nDims));
  LVec result = createKArray<LVec>(U.extent(0));

  if (nDims==1)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& nx = normal.x();
    const ScalarT u = rhou/rho;
    const ScalarT vn = u*nx;
    const ScalarT p = pressure(U,gamma,nDims);
    const ScalarT a = c(U,gamma,nDims);
    const ScalarT h = a*a/(gamma-1.0);
    const ScalarT e_k = 0.5*(u*u);
    const ScalarT h0 = h + e_k;

    result(0) = rho*vn;
    result(1) = rho*u*vn + p*nx;
    result(2) = rho*h0*vn;
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT p = pressure(U,gamma,nDims);
    ScalarT a = c(U,gamma,nDims);
    ScalarT h = a*a/(gamma-1.0);
    ScalarT e_k = 0.5*(u*u+v*v);
    ScalarT h0 = h + e_k;

    result(0) = rho*vn;
    result(1) = rho*u*vn + p*nx;
    result(2) = rho*v*vn + p*ny;
    result(3) = rho*h0*vn;
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhow = U(3);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    const ScalarT& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT p = pressure(U,gamma,nDims);
    ScalarT a = c(U,gamma,nDims);
    ScalarT h = a*a/(gamma-1.0);
    ScalarT e_k = 0.5*(u*u+v*v+w*w);
    ScalarT h0 = h + e_k;

    result(0) = rho*vn;
    result(1) = rho*u*vn + p*nx;
    result(2) = rho*v*vn + p*ny;
    result(3) = rho*w*vn + p*nz;
    result(4) = rho*h0*vn;
  }

  return result;
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LVec
Equation_Euler<ScalarT>::
internal::lambda(const LVec& U, const SPoint& normal,
                 const double gamma, const int nDims)
{
  assert(normal.size()==static_cast<std::size_t>(nDims));
  LVec result = createKArray<LVec>(U.size());

  if (nDims==1)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& nx = normal.x();
    ScalarT u = rhou/rho;
    ScalarT vn = u*nx;
    ScalarT a = c(U,gamma,nDims);
    result(0) = vn-a;
    result(1) = vn;
    result(2) = vn+a;
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT a = c(U,gamma,nDims);
    result(0) = vn-a;
    result(1) = vn;
    result(2) = vn+a;
    result(3) = vn;
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhow = U(3);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    const ScalarT& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT a = c(U,gamma,nDims);
    result(0) = vn-a;
    result(1) = vn;
    result(2) = vn+a;
    result(3) = vn;
    result(4) = vn;
  }

  return result;
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat 
Equation_Euler<ScalarT>::
internal::lambda_mat(const LVec& U, const SPoint& normal,
                     const double gamma, const int nDims)
{
  assert(normal.size()==static_cast<std::size_t>(nDims));
  LMat result = createKArray<LMat>(U.size(),U.size());

  if (nDims==1)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& nx = normal.x();
    ScalarT u  = rhou/rho;
    ScalarT vn = u*nx;
    ScalarT a = c(U,gamma,nDims);

    result(0,0) = std::fabs(vn-a);
    result(1,1) = std::fabs(vn);
    result(2,2) = std::fabs(vn+a);
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    ScalarT u  = rhou/rho;
    ScalarT v  = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT a = c(U,gamma,nDims);

    result(0,0) = std::fabs(vn-a);
    result(1,1) = std::fabs(vn);
    result(2,2) = std::fabs(vn+a);
    result(3,3) = std::fabs(vn);
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& rhow = U(3);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    const ScalarT& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT a = c(U,gamma,nDims);
    result(0,0) = std::fabs(vn-a);
    result(1,1) = std::fabs(vn);
    result(2,2) = std::fabs(vn+a);
    result(3,3) = std::fabs(vn);
    result(4,4) = std::fabs(vn);
  }

  return result;
}

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::
internal::lambda_max(const LVec& U, const SPoint& normal, const double gamma, const int nDims)
{
  assert(normal.size()==Teuchos::as<std::size_t>(nDims));
  ScalarT lambda_max_ = 0.0;

  if (nDims==1)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& nx = normal.x();
    ScalarT u  = rhou/rho;
    ScalarT vn = u*nx;
    ScalarT a = c(U,gamma,nDims)*normal.norm();
    lambda_max_ = vn + a;
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U(0);
    const ScalarT& rhou = U(1);
    const ScalarT& rhov = U(2);
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    ScalarT u  = rhou/rho;
    ScalarT v  = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT a = c(U,gamma,nDims)*normal.norm();
    lambda_max_ = vn + a;
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U[0];
    const ScalarT& rhou = U[1];
    const ScalarT& rhov = U[2];
    const ScalarT& rhow = U[3];
    const ScalarT& nx = normal.x();
    const ScalarT& ny = normal.y();
    const ScalarT& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT a = c(U,gamma,nDims)*normal.norm();
    lambda_max_ = vn + a;
  }
  return lambda_max_;
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat Equation_Euler<ScalarT>::
internal::eigenvector(const LVec& U, const SPoint& normal, const double gamma, const int nDims)
{
  assert(normal.size()==static_cast<std::size_t>(nDims));
  LMat y = createKArray<LMat>(U.size(),U.size());

  if (nDims==1)
  {
    std::cout << "The case for 1D eigenvector has not been implemented yet" << std::endl;
    exit(1);
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U[0];
    const ScalarT& rhou = U[1];
    const ScalarT& rhov = U[2];
    const double& nx = normal.x();
    const double& ny = normal.y();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT a = c(U,gamma,nDims);
    ScalarT h = a*a/(gamma-1.0);
    ScalarT e_k = 0.5*(u*u+v*v);
    ScalarT h0 = h + e_k;

    y(0,0) = 1.0;
    y(0,1) = 1.0;
    y(0,2) = 1.0;
    y(0,3) = 0.0;

    y(1,0) = u - a*nx;
    y(1,1) = u;
    y(1,2) = u + a*nx;
    y(1,3) = ny;

    y(2,0) = v - a*ny;
    y(2,1) = v;
    y(2,2) = v + a*ny;
    y(2,3) = -nx;

    y(3,0) = h0 - a*vn;
    y(3,1) = e_k;
    y(3,2) = h0 + a*vn;
    y(3,3) = u*ny - v*nx;
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U[0];
    const ScalarT& rhou = U[1];
    const ScalarT& rhov = U[2];
    const ScalarT& rhow = U[3];
    const double& nx = normal.x();
    const double& ny = normal.y();
    const double& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT a = c(U,gamma,nDims);
    ScalarT h = a*a/(gamma-1.0);
    ScalarT e_k = 0.5*(u*u+v*v+w*w);
    ScalarT h0 = h + e_k;

    if (nx!=0.0) // R1
    {

      y(0,0) = 1.0;
      y(0,1) = 1.0;
      y(0,2) = 1.0;
      y(0,3) = 0.0;
      y(0,4) = 0.0;

      y(1,0) = u - a*nx;
      y(1,1) = u;
      y(1,2) = u + a*nx;
      y(1,3) = ny;
      y(1,4) = -nz;

      y(2,0) = v - a*ny;
      y(2,1) = v;
      y(2,2) = v + a*ny;
      y(2,3) = -nx;
      y(2,4) = 0.;

      y(3,0) = w - a*nz;
      y(3,1) = w;
      y(3,2) = w + a*nz;
      y(3,3) = 0.;
      y(3,4) = nx;

      y(4,0) = h0 - a*vn;
      y(4,1) = e_k;
      y(4,2) = h0 + a*vn;
      y(4,3) = u*ny - v*nx;
      y(4,4) = w*nx - u*nz;

    } else if (ny!=0.0) { // R2

      y(0,0) = 1.0;
      y(0,1) = 1.0;
      y(0,2) = 1.0;
      y(0,3) = 0.0;
      y(0,4) = 0.0;

      y(1,0) = u - a*nx;
      y(1,1) = u;
      y(1,2) = u + a*nx;
      y(1,3) = ny;
      y(1,4) = 0.0;

      y(2,0) = v - a*ny;
      y(2,1) = v;
      y(2,2) = v + a*ny;
      y(2,3) = -nx;
      y(2,4) = nz;

      y(3,0) = w - a*nz;
      y(3,1) = w;
      y(3,2) = w + a*nz;
      y(3,3) = 0.0;
      y(3,4) = -ny;

      y(4,0) = h0 - a*vn;
      y(4,1) = e_k;
      y(4,2) = h0 + a*vn;
      y(4,3) = u*ny - v*nx;
      y(4,4) = v*nz - w*ny;

    } else if (nz!=0.0) { // R3

      y(0,0) = 1.0;
      y(0,1) = 1.0;
      y(0,2) = 1.0;
      y(0,3) = 0.0;
      y(0,4) = 0.0;

      y(1,0) = u - a*nx;
      y(1,1) = u;
      y(1,2) = u + a*nx;
      y(1,3) = -nz;
      y(1,4) = 0.0;

      y(2,0) = v - a*ny;
      y(2,1) = v;
      y(2,2) = v + a*ny;
      y(2,3) = 0.0;
      y(2,4) = nz;

      y(3,0) = w - a*nz;
      y(3,1) = w;
      y(3,2) = w + a*nz;
      y(3,3) = nx;
      y(3,4) = -ny;

      y(4,0) = h0 - a*vn;
      y(4,1) = e_k;
      y(4,2) = h0 + a*vn;
      y(4,3) = w*nx - u*nz;
      y(4,4) = v*nz - w*ny;

    }
  }

  return y;
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat Equation_Euler<ScalarT>::
internal::inverse_eigenvector(const LVec& U, const SPoint& normal, const double gamma, const int nDims)
{
  assert(normal.size()==static_cast<std::size_t>(nDims));
  LMat y = createKArray<LMat>(U.size(),U.size());

  if (nDims==1)
  {
    std::cout << "The case for 1D eigenvector has not been implemented yet" << std::endl;
    exit(1);
  }

  if (nDims==2)
  {
    const ScalarT& rho  = U[0];
    const ScalarT& rhou = U[1];
    const ScalarT& rhov = U[2];
    const double& nx = normal.x();
    const double& ny = normal.y();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT vn = u*nx+v*ny;
    ScalarT a = c(U,gamma,nDims);
    ScalarT e_k = 0.5*(u*u+v*v);

    y(0,0) = ((gamma-1.0)*e_k + a*vn)/(2.0*a*a);
    y(0,1) = ((1.0-gamma)*u-a*nx)/(2.0*a*a);
    y(0,2) = ((1.0-gamma)*v-a*ny)/(2.0*a*a);
    y(0,3) = (gamma-1.0)/(2.0*a*a);

    y(1,0) = (a*a - (gamma-1.0)*e_k)/(a*a);
    y(1,1) = ((gamma-1.0)*u)/(a*a);
    y(1,2) = ((gamma-1.0)*v)/(a*a);
    y(1,3) = (1.0-gamma)/(a*a);

    y(2,0) = ( (gamma-1.0)*e_k - a*vn )/(2.0*a*a);
    y(2,1) = ( (1.0-gamma)*u + a*nx )/(2.0*a*a);
    y(2,2) = ( (1.0-gamma)*v + a*ny )/(2.0*a*a);
    y(2,3) = (gamma-1.0)/(2.0*a*a);

    y(3,0) = v*nx - u*ny;
    y(3,1) = ny;
    y(3,2) = -nx;
    y(3,3) = 0.0;
  }

  if (nDims==3)
  {
    const ScalarT& rho  = U[0];
    const ScalarT& rhou = U[1];
    const ScalarT& rhov = U[2];
    const ScalarT& rhow = U[3];
    const double& nx = normal.x();
    const double& ny = normal.y();
    const double& nz = normal.z();
    ScalarT u = rhou/rho;
    ScalarT v = rhov/rho;
    ScalarT w = rhow/rho;
    ScalarT vn = u*nx+v*ny+w*nz;
    ScalarT a = c(U,gamma,nDims);
    ScalarT e_k = 0.5*(u*u+v*v+w*w);

    if (nx!=0.0) // L1 = inv(R1)
    {

      y(0,0) = ((gamma-1.0)*e_k+a*vn)/(2*a*a);
      y(0,1) = ((1.0-gamma)*u-a*nx)/(2*a*a);
      y(0,2) = ((1.0-gamma)*v-a*ny)/(2*a*a);
      y(0,3) = ((1.0-gamma)*w-a*nz)/(2*a*a);
      y(0,4) = (gamma-1.0)/(2*a*a);

      y(1,0) = (a*a-(gamma-1.0)*e_k)/(a*a);
      y(1,1) = ((gamma-1.0)*u)/(a*a);
      y(1,2) = ((gamma-1.0)*v)/(a*a);
      y(1,3) = ((gamma-1.0)*w)/(a*a);
      y(1,4) = (1.0-gamma)/(a*a);

      y(2,0) = ((gamma-1.0)*e_k-a*vn)/(2*a*a);
      y(2,1) = ((1.0-gamma)*u+a*nx)/(2*a*a);
      y(2,2) = ((1.0-gamma)*v+a*ny)/(2*a*a);
      y(2,3) = ((1.0-gamma)*w+a*nz)/(2*a*a);
      y(2,4) = (gamma-1.0)/(2*a*a);

      y(3,0) = (v-vn*ny)/nx;
      y(3,1) = ny;
      y(3,2) = (ny*ny-1.0)/nx;
      y(3,3) = ny*nz/nx;
      y(3,4) = 0.0;

      y(4,0) = (vn*nz-w)/nx;
      y(4,1) = -nz;
      y(4,2) = -ny*nz/nx;
      y(4,3) = (1.0-nz*nz)/nx;
      y(4,4) = 0.0;

    } else if (ny!=0.0) { // L2 = inv(R2)

      y(0,0) = ((gamma-1.0)*e_k+a*vn)/(2*a*a);
      y(0,1) = ((1.0-gamma)*u-a*nx)/(2*a*a);
      y(0,2) = ((1.0-gamma)*v-a*ny)/(2*a*a);
      y(0,3) = ((1.0-gamma)*w-a*nz)/(2*a*a);
      y(0,4) = (gamma-1.0)/(2*a*a);

      y(1,0) = (a*a-(gamma-1.0)*e_k)/(a*a);
      y(1,1) = ((gamma-1.0)*u)/(a*a);
      y(1,2) = ((gamma-1.0)*v)/(a*a);
      y(1,3) = ((gamma-1.0)*w)/(a*a);
      y(1,4) = (1.0-gamma)/(a*a);

      y(2,0) = ((gamma-1.0)*e_k-a*vn)/(2*a*a);
      y(2,1) = ((1.0-gamma)*u+a*nx)/(2*a*a);
      y(2,2) = ((1.0-gamma)*v+a*ny)/(2*a*a);
      y(2,3) = ((1.0-gamma)*w+a*nz)/(2*a*a);
      y(2,4) = (gamma-1.0)/(2*a*a);

      y(3,0) = (vn*nx-u)/ny;
      y(3,1) = (1.0-nx*nx)/ny;
      y(3,2) = -nx;
      y(3,3) = -nx*nz/ny;
      y(3,4) = 0.0;

      y(4,0) = (w-vn*nz)/ny;
      y(4,1) = nx*nz/ny;
      y(4,2) = -nz;
      y(4,3) = (nz*nz-1.0)/ny;
      y(4,4) = 0.0;

    } else if (nz!=0.0) { // L3 = inv(R3)

      y(0,0) = ((gamma-1.0)*e_k+a*vn)/(2*a*a);
      y(0,1) = ((1.0-gamma)*u-a*nx)/(2*a*a);
      y(0,2) = ((1.0-gamma)*v-a*ny)/(2*a*a);
      y(0,3) = ((1.0-gamma)*w-a*nz)/(2*a*a);
      y(0,4) = (gamma-1.0)/(2*a*a);

      y(1,0) = (a*a-(gamma-1.0)*e_k)/(a*a);
      y(1,1) = ((gamma-1.0)*u)/(a*a);
      y(1,2) = ((gamma-1.0)*v)/(a*a);
      y(1,3) = ((gamma-1.0)*w)/(a*a);
      y(1,4) = (1.0-gamma)/(a*a);

      y(2,0) = ((gamma-1.0)*e_k-a*vn)/(2*a*a);
      y(2,1) = ((1.0-gamma)*u+a*nx)/(2*a*a);
      y(2,2) = ((1.0-gamma)*v+a*ny)/(2*a*a);
      y(2,3) = ((1.0-gamma)*w+a*nz)/(2*a*a);
      y(2,4) = (gamma-1.0)/(2*a*a);

      y(3,0) = (u - vn*nx)/nz;
      y(3,1) = (nx*nx-1.0)/nz;
      y(3,2) = nx*ny/nz;
      y(3,3) = nx;
      y(3,4) = 0.0;

      y(4,0) = (vn*ny - v)/nz;
      y(4,1) = -nx*ny/nz;
      y(4,2) = (1.0-ny*ny)/nz;
      y(4,3) = -ny;
      y(4,4) = 0.0;

    }
  }

  return y;
}

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::
internal::c(const LVec& U, const double gamma, const int nDims)
{
  const ScalarT& rho = U[0];
  const ScalarT p = pressure(U,gamma,nDims);
  return std::sqrt(gamma*p/rho);
}

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::internal::H(const LVec& U, const double gamma, const int nDims)
{
  const ScalarT& rho = U[0];
  const ScalarT& rhoE = U[U.size()-1];
  return rhoE/rho + pressure(U,gamma,nDims)/rho;
}

template <typename ScalarT>
ScalarT Equation_Euler<ScalarT>::internal::S(const LVec& U, const double gamma, const int nDims)
{
  const ScalarT& rho = U[0];
  return std::log(pressure(U,gamma,nDims)/std::pow(rho,gamma));
}

template <typename ScalarT>
void Equation_Euler<ScalarT>::internal::conservedtoroe(LVec& V, const LVec& U, const double gamma, const int nDims)
{
  const ScalarT& rho  = U[0];
  const ScalarT& rhou = U[1];
  auto u = rhou/rho;
  auto Hs = H(U,gamma,nDims);
  auto cs = c(U,gamma,nDims);

  V = createKArray<LVec>(U.size());
  if (nDims==1)
  {
    V[0] = u;
    V[2] = Hs;
    V[3] = cs;
  }
  if (nDims==2)
  {
    const ScalarT& rhov = U[2];
    auto v = rhov/rho;
    V[0] = u;
    V[1] = v;
    V[2] = Hs;
    V[3] = cs;
  }
  if (nDims==3)
  {
    const ScalarT& rhov = U[2];
    auto v = rhov/rho;
    const ScalarT& rhow = U[3];
    auto w = rhow/rho;
    V[0] = u;
    V[1] = v;
    V[2] = w;
    V[3] = Hs;
    V[4] = cs;
  }
}

template <typename ScalarT>
void Equation_Euler<ScalarT>::internal::
roe_mean(LVec& Ulr, const LVec& Ul, const LVec& Ur, const double gamma, const int nDims)
{
  Ulr = createKArray<LVec>(Ul.size());
  if (nDims==1)
  {
    const ScalarT& rhol = Ul[0];
    const ScalarT& rhoul = Ul[1];
    const ScalarT ul = rhoul/rhol;
    const ScalarT Hl = H(Ul,gamma,nDims);

    const ScalarT& rhor = Ur[0];
    const ScalarT& rhour = Ur[1];
    const ScalarT ur = rhour/rhor;
    const ScalarT Hr = H(Ur,gamma,nDims);

    const ScalarT rholr = std::sqrt(rhol*rhor);
    const ScalarT ulr = (std::sqrt(rhol)*ul + std::sqrt(rhor)*ur)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT Hlr = (std::sqrt(rhol)*Hl + std::sqrt(rhor)*Hr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT clr = std::sqrt( (gamma-1.0)*(Hlr - 0.5 * (ulr*ulr)) );
    const ScalarT plr = clr*clr*rholr/gamma;
    const ScalarT rhoulr = rholr*ulr;
    const ScalarT rhoElr = rholr*Hlr-plr;

    Ulr[0] = rholr;
    Ulr[1] = rhoulr;
    Ulr[2] = rhoElr;
  }

  if (nDims==2)
  {
    const ScalarT& rhol = Ul[0];
    const ScalarT& rhoul = Ul[1];
    const ScalarT& rhovl = Ul[2];
    const ScalarT ul = rhoul/rhol;
    const ScalarT vl = rhovl/rhol;
    const ScalarT Hl = H(Ul,gamma,nDims);

    const ScalarT& rhor = Ur[0];
    const ScalarT& rhour = Ur[1];
    const ScalarT& rhovr = Ur[2];
    const ScalarT ur = rhour/rhor;
    const ScalarT vr = rhovr/rhor;
    const ScalarT Hr = H(Ur,gamma,nDims);

    const ScalarT rholr = std::sqrt(rhol*rhor);
    const ScalarT ulr = (std::sqrt(rhol)*ul + std::sqrt(rhor)*ur)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT vlr = (std::sqrt(rhol)*vl + std::sqrt(rhor)*vr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT Hlr = (std::sqrt(rhol)*Hl + std::sqrt(rhor)*Hr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT clr = std::sqrt( (gamma-1.0)*(Hlr - 0.5 * (ulr*ulr+vlr*vlr)) );
    const ScalarT plr = clr*clr*rholr/gamma;
    const ScalarT rhoulr = rholr*ulr;
    const ScalarT rhovlr = rholr*vlr;
    const ScalarT rhoElr = rholr*Hlr-plr;

    Ulr[0] = rholr;
    Ulr[1] = rhoulr;
    Ulr[2] = rhovlr;
    Ulr[3] = rhoElr;
  }

  if (nDims==3)
  {
    const ScalarT& rhol = Ul[0];
    const ScalarT& rhoul = Ul[1];
    const ScalarT& rhovl = Ul[2];
    const ScalarT& rhowl = Ul[3];
    const ScalarT ul = rhoul/rhol;
    const ScalarT vl = rhovl/rhol;
    const ScalarT wl = rhowl/rhol;
    const ScalarT Hl = H(Ul,gamma,nDims);

    const ScalarT& rhor = Ur[0];
    const ScalarT& rhour = Ur[1];
    const ScalarT& rhovr = Ur[2];
    const ScalarT& rhowr = Ur[3];
    const ScalarT ur = rhour/rhor;
    const ScalarT vr = rhovr/rhor;
    const ScalarT wr = rhowr/rhor;
    const ScalarT Hr = H(Ur,gamma,nDims);

    const ScalarT rholr = std::sqrt(rhol*rhor);
    const ScalarT ulr = (std::sqrt(rhol)*ul + std::sqrt(rhor)*ur)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT vlr = (std::sqrt(rhol)*vl + std::sqrt(rhor)*vr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT wlr = (std::sqrt(rhol)*wl + std::sqrt(rhor)*wr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT Hlr = (std::sqrt(rhol)*Hl + std::sqrt(rhor)*Hr)/(std::sqrt(rhol) + std::sqrt(rhor));
    const ScalarT clr = std::sqrt( (gamma-1.0)*(Hlr - 0.5 * (ulr*ulr+vlr*vlr+wlr*wlr)) );
    const ScalarT plr = clr*clr*rholr/gamma;
    const ScalarT rhoulr = rholr*ulr;
    const ScalarT rhovlr = rholr*vlr;
    const ScalarT rhowlr = rholr*wlr;
    const ScalarT rhoElr = rholr*Hlr-plr;

    Ulr[0] = rholr;
    Ulr[1] = rhoulr;
    Ulr[2] = rhovlr;
    Ulr[3] = rhowlr;
    Ulr[4] = rhoElr;
  }
}

template <typename ScalarT>
typename Equation_Euler<ScalarT>::LMat 
Equation_Euler<ScalarT>::internal::
abs_flux_u(const LVec& U, const SPoint& normal,
           const double gamma, const int nDims)
{
  // make sure normal has same length as num of dimensions
  assert(normal.size()==static_cast<std::size_t>(nDims));

  auto R = eigenvector(U,normal,gamma,nDims);
  auto Lambda = lambda_mat(U,normal,gamma,nDims);
  auto inv_R = inverse_eigenvector(U,normal,gamma,nDims);

  const int nEq = U.size();
  LMat tmp = createKArray<LMat>(nEq,nEq), result = createKArray<LMat>(nEq,nEq);
  for (int i = 0; i < nEq; ++i)
  {
    for (int j = 0; j < nEq; ++j)
    {
      ScalarT s (0.0);
      for (int k = 0; k < nEq; ++k)
      {
        s += Lambda(i,k)*inv_R(k,j);
      }
      tmp(i,j) = s;
    }
  }

  for (int i = 0; i < nEq; ++i)
  {
    for (int j = 0; j < nEq; ++j)
    {
      ScalarT s (0.0);
      for (int k = 0; k < nEq; ++k)
      {
        s += R(i,k)*tmp(k,j);
      }
      result(i,j) = s;
    }
  }

  return result;
}

}
// end namespace hydrofem

#endif /** __Hydrofem_Equation_Euler_impl_HPP__ */
