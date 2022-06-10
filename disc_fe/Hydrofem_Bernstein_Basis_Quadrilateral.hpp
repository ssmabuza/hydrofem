// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Bernstein_Basis_Quadrilateral_HPP__
#define __Hydrofem_Bernstein_Basis_Quadrilateral_HPP__

#include "Hydrofem.hpp"

#include "Hydrofem_SPoint.hpp"
#include "Hydrofem_FEBasis.hpp"
#include "Hydrofem_Bernstein_Basis_Line.hpp"
#include "Hydrofem_ReferenceQuadrilateral.hpp"

namespace hydrofem
{

class Bernstein_Basis_Quadrilateral
  :
  public FEBasis
{
  
  std::shared_ptr<Bernstein_Basis_Line> line_basis1;
  std::shared_ptr<Bernstein_Basis_Line> line_basis2;
  int m_i;
  int m_j;
  int m_p;
  
public:
  
  Bernstein_Basis_Quadrilateral(const int order,
                                const int i,
                                const int j)
  {
    line_basis1 = std::make_shared<Bernstein_Basis_Line>(order,i);
    line_basis2 = std::make_shared<Bernstein_Basis_Line>(order,j);
    m_i = i;
    m_j = j;
    m_p = order;
  }
                                              
  virtual ~Bernstein_Basis_Quadrilateral() {}
  
  virtual inline int p() const override { return m_p; }
  
  virtual double operator()(const SPoint& point, const std::vector<SPoint>& element) const override
  {
    const SPoint ref_pt = physicalQuadToReference(element.at(0),element.at(1),element.at(2),element.at(3),point);
    return eval(ref_pt.x(),ref_pt.y(),m_i,m_j);
  }
  
  virtual SPoint grad(const SPoint& point, const std::vector<SPoint>& element) const override
  {
    const SPoint ref_pt = physicalQuadToReference(element.at(0),element.at(1),element.at(2),element.at(3),point);
    return SPoint(eval_x(ref_pt.x(),ref_pt.y(),m_i,m_j),eval_y(ref_pt.x(),ref_pt.y(),m_i,m_j));
  }
  
  using FEBasisPtr = std::shared_ptr<FEBasis>;
  static std::shared_ptr<std::vector<FEBasisPtr>>
  buildElementBasis(const int order)
  {
    std::shared_ptr<std::vector<FEBasisPtr>> basis = std::make_shared<std::vector<FEBasisPtr>>((order+1)*(order+1));
    int index = 0;
    for (int i = 0; i <= order; ++i)
      for (int j = 0; j <= order; ++j)
      {
        basis->at(index) = std::make_shared<Bernstein_Basis_Quadrilateral>(order,i,j);
        assert(basis->at(index));
        ++index;
      }
    for (auto it : *basis)
      assert(it);
    return basis;
  }
  
private:
    
  /// evaluation on a reference unit element
  double eval(const double x,
              const double y,
              const int i,
              const int j) const;

  /// evaluation on arbitrary quadrilateral
  double eval(const double x,
              const double y,
              const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const int i,
              const int j) const;

  /// x derivative on reference element
  double eval_x(const double x,
                const double y,
                const int i,
                const int j) const;

  /// x derivative on arbitrary quadrilateral              
  double eval_x(const double x,
                const double y,
                const SPoint& quad1,
                const SPoint& quad2,
                const SPoint& quad3,
                const SPoint& quad4,
                const int i,
                const int j) const;

  /// y derivative on reference element
  double eval_y(const double x,
                const double y,
                const int i,
                const int j) const;

  /// y derivative on arbitrary quadrilateral              
  double eval_y(const double x,
                const double y,
                const SPoint& quad1,
                const SPoint& quad2,
                const SPoint& quad3,
                const SPoint& quad4,
                const int i,
                const int j) const;

  SPoint grad_eval(const double x,
              const double y,
              const SPoint& quad1,
              const SPoint& quad2,
              const SPoint& quad3,
              const SPoint& quad4,
              const int i,
              const int j) const;

  template <typename ScalarT>
  inline ScalarT evaluateFunction(const ScalarT x,
                                  const ScalarT y,
                                  const SPoint& a,
                                  const SPoint& b,
                                  const SPoint& c,
                                  const SPoint& d,
                                  const int i,
                                  const int j) const
  {
    ScalarT u(0.), v;
    
    ScalarT C = (ScalarT)(a.y() - y) * (d.x() - x) - (ScalarT)(a.x() - x) * (d.y() - y);
    ScalarT B = (ScalarT)(a.y() - y) * (c.x() - d.x()) + (ScalarT)(b.y() - a.y()) * (d.x() - x) - 
                (ScalarT)(a.x() - x) * (c.y() - d.y()) - (ScalarT)(b.x() - a.x()) * (d.y() - y);
    ScalarT A = (ScalarT)(b.y() - a.y()) * (c.x() - d.x()) - (ScalarT)(b.x() - a.x()) * (c.y() - d.y());
    
    ScalarT D = B * B - 4 * A * C;

    u = (-B - std::sqrt(D)) / (2 * A);

    ScalarT p1x = a.x() + (b.x() - a.x()) * u;
    ScalarT p2x = d.x() + (c.x() - d.x()) * u;
    ScalarT px = x;

    v = (px - p1x) / (p2x - p1x);  

    return line_basis1->template evaluateFunction<ScalarT>(u,i) * line_basis2->template evaluateFunction<ScalarT>(v,j);
  }
                
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Bernstein_Basis_Quadrilateral_HPP__ */
