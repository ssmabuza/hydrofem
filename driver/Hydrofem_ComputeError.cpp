// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#include "Hydrofem_ComputeError.hpp"


namespace hydrofem
{

double
ComputeError::evaluateScalarFieldError(const std::shared_ptr<const std::vector<std::shared_ptr<FEBasis>>>& basis,
                                       const std::shared_ptr<const std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                                       const std::shared_ptr<const DofMapper>& dofmapper, 
                                       const std::shared_ptr<const FEVector>& u_numerical,
                                       const std::shared_ptr<const std::function<double(SPoint)>>& u_exact,
                                       const std::shared_ptr<const std::function<SPoint(SPoint)>>& grad_u_exact,
                                       Norm normt)
{
  double err = 0.0;
  if (normt == Norm::Linf) err = std::numeric_limits<double>::min();
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    double err_e = 0.0;
    if (normt == Norm::Linf) err_e = std::numeric_limits<double>::min();
    const auto& quadrature_e = * quadrature->at(elem_ind);
    const auto& quad_p = quadrature_e.get_q_points();
    const auto& quad_w = quadrature_e.get_q_weights();
    const auto element = dofmapper->mesh()->getElementVertices(elem_ind);
    const auto& glinds = dofmapper->getGlobDofIndexes(elem_ind);
    const auto& linds = dofmapper->getLocDofIndexes();
    
    // integrate over each cell 
    for (int qp = 0; qp < quadrature_e.size(); ++qp)
    {
      // compute uh(x_q)
      double uh = 0.0;
      for (std::size_t i = 0; i < linds.size(); ++i)
        uh += (*u_numerical)[glinds[i]] * (* (basis->at(linds[i])) )(quad_p.at(qp),element);
      // compute the error at quad point
      if (normt == Norm::L1)
        err_e += quad_w[qp] * std::fabs(uh - (*u_exact)(quad_p.at(qp)));
      else if (normt == Norm::L2)
        err_e += quad_w[qp] * std::pow(uh - (*u_exact)(quad_p.at(qp)),2.0);
      else if (normt == Norm::Linf)
        err_e = std::max(err_e, std::fabs(uh - (*u_exact)(quad_p.at(qp))));
      else if (normt == Norm::H1) {
        SPoint grad_uh(0.0,0.0);
        for (std::size_t i = 0; i < linds.size(); ++i)
          grad_uh += (*u_numerical)[glinds[i]] * (* (basis->at(linds[i])) ).grad(quad_p.at(qp),element);
        err_e += quad_w[qp] * std::pow(uh - (*u_exact)(quad_p.at(qp)),2.0);
        if (grad_u_exact)
          err_e += quad_w[qp] * (grad_uh - (*grad_u_exact)(quad_p.at(qp))) * (grad_uh - (*grad_u_exact)(quad_p.at(qp)));
        else {
          throw std::runtime_error("Error when computing H1 norm, gradient of exact solution not provided.");
        }
      }
    }
    
    if ((normt == Norm::L1) || (normt == Norm::L2) || (normt == Norm::H1))
      err += err_e;
    else if (normt == Norm::Linf)
      err = std::max(err,err_e);
  }
  if ((normt == Norm::L2) || (normt == Norm::H1))
    err = std::sqrt(err);
  return err;
}

inline static SPoint getSPointFromLVec(const LVEC_<double>& vec)
{
  assert((vec.size() <= 3) && (vec.size() > 0));
  SPoint pt(static_cast<int>(vec.size()));
  pt.x() = vec[0];
  if (vec.size() > 1)
  {
    pt.y() = vec[1];
    if (vec.size() > 2)
      pt.z() = vec[2];
  }
  return pt;
}

double
ComputeError::evaluateScalarFieldError(const std::shared_ptr<const std::vector<std::shared_ptr<FEBasis>>>& basis,
                                       const std::shared_ptr<const std::vector<std::shared_ptr<Quadrature>>>& quadrature,
                                       const std::shared_ptr<const DofMapper>& dofmapper,
                                       const std::shared_ptr<const FEVector>& u_numerical,
                                       const std::shared_ptr<ScalarAnalyticalExpression>& u_exact,
                                       const std::shared_ptr<VectorAnalyticalExpression>& grad_u_exact,
                                       Norm normt)
{
  double err = 0.0;
  if (normt == Norm::Linf) err = std::numeric_limits<double>::min();
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    double err_e = 0.0;
    if (normt == Norm::Linf) err_e = std::numeric_limits<double>::min();
    const auto& quadrature_e = * quadrature->at(elem_ind);
    const auto& quad_p = quadrature_e.get_q_points();
    const auto& quad_w = quadrature_e.get_q_weights();
    const auto element = dofmapper->mesh()->getElementVertices(elem_ind);
    const auto& glinds = dofmapper->getGlobDofIndexes(elem_ind);
    const auto& linds = dofmapper->getLocDofIndexes();

    // integrate over each cell
    for (int qp = 0; qp < quadrature_e.size(); ++qp)
    {
      // compute uh(x_q)
      double uh = 0.0;
      for (std::size_t i = 0; i < linds.size(); ++i)
        uh += (*u_numerical)[glinds[i]] * (* (basis->at(linds[i])) )(quad_p.at(qp),element);
      // compute the error at quad point
      if (normt == Norm::L1)
        err_e += quad_w[qp] * std::fabs(uh - (*u_exact)(quad_p.at(qp)));
      else if (normt == Norm::L2)
        err_e += quad_w[qp] * std::pow(uh - (*u_exact)(quad_p.at(qp)),2.0);
      else if (normt == Norm::Linf)
        err_e = std::max(err_e, std::fabs(uh - (*u_exact)(quad_p.at(qp))));
      else if (normt == Norm::H1) {
        SPoint grad_uh(0.0,0.0);
        for (std::size_t i = 0; i < linds.size(); ++i)
          grad_uh += (*u_numerical)[glinds[i]] * (* (basis->at(linds[i])) ).grad(quad_p.at(qp),element);
        err_e += quad_w[qp] * std::pow(uh - (*u_exact)(quad_p.at(qp)),2.0);
        if (grad_u_exact)
          err_e += quad_w[qp] * (grad_uh - getSPointFromLVec((*grad_u_exact)(quad_p.at(qp))) ) * (grad_uh - getSPointFromLVec((*grad_u_exact)(quad_p.at(qp))) );
        else {
          throw std::runtime_error("Error when computing H1 norm, gradient of exact solution not provided.");
        }
      }
    }

    if ((normt == Norm::L1) || (normt == Norm::L2) || (normt == Norm::H1))
      err += err_e;
    else if (normt == Norm::Linf)
      err = std::max(err,err_e);
  }
  if ((normt == Norm::L2) || (normt == Norm::H1))
    err = std::sqrt(err);
  return err;
}

}
// end namespace hydrofem
