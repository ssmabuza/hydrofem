// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "InitialSolution_ReactiveFlow.hpp"
#include "TransformWallToChannel.hpp"

namespace valiant
{

void InitialSolution_ReactiveFlow::evaluate()
{
  m_result = Teuchos::rcp(new TVector(m_lob->getOwnedMap()));
  m_result->putScalar(0.0);
  const auto& m_result_view = m_result->getDataNonConst();
  for (int i = 0; i < m_mesh_f->numOfPoints(); ++i)
    m_result_view[i] = m_ic->evaluate(m_mesh_f->getPoint(i))(0);
  for (int i = 0; i < m_mesh_w->numOfPoints(); ++i)
    m_result_view[i + m_mesh_f->numOfPoints()] = m_ic->evaluate(m_mesh_w->getPoint(i))(1);
  m_gather->doGather(m_result_gathered->at(0),m_result_gathered->at(1),m_result);
  m_is_computed = true;
}

}
// end namespace valiant
