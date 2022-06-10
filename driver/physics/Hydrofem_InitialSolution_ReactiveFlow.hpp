// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2022) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Valiant_InitialSolution_ReactiveFlow_HPP__
#define __Valiant_InitialSolution_ReactiveFlow_HPP__

#include "Mesh.hpp"
#include "InitialSolution.hpp"
#include "InitialCondition.hpp"
#include "LOB_ReactiveFlow.hpp"
#include "BoundaryFaceInfoSetter.hpp"
#include "GlobalGather_ReactiveFlow.hpp"

namespace valiant
{


class InitialSolution_ReactiveFlow
  :
  public InitialSolution
{
public:
  
  InitialSolution_ReactiveFlow(const Teuchos::RCP<Mesh>& mesh_f,
                               const Teuchos::RCP<Mesh>& mesh_w,
                               const Teuchos::RCP<InitialCondition>& ic,
                               const Teuchos::RCP<LinearObjectBuilder<int>>& lob,
                               const Teuchos::RCP<GlobalGather_ReactiveFlow>& gather)
    :
    InitialSolution(lob)
  {
    using Teuchos::rcp_dynamic_cast;
    
    m_mesh_f = mesh_f;
    m_mesh_w = mesh_w;
    m_ic = ic;
    
    // allocate the gathered solution
    m_result_gathered = Teuchos::rcp(new std::vector<MDArray>(2));
    m_result_gathered->at(0) = createKArray<MDArray>(m_mesh_f->numOfElements(),m_mesh_f->getElement(0).getNodes().size());
    m_result_gathered->at(1) = createKArray<MDArray>(m_mesh_w->numOfElements(),m_mesh_w->getElement(0).getNodes().size());
    for (auto& it : *m_result_gathered) zeroOutArray(it);
    
    // allocate the distributed IC
    m_result = Tpetra::createVector<double>(lob->getOwnedMap()); m_result->putScalar(0.0);
  }
  
  ~InitialSolution_ReactiveFlow() override = default;

  void evaluate() override;
  
  /** \brief gets the global scatter for the system */
  Teuchos::RCP<GlobalGather_ReactiveFlow>
  getGlobalSystemScatter() const
  { return m_gather; }
  
private:
  
  Teuchos::RCP<Mesh> m_mesh_f;
  Teuchos::RCP<Mesh> m_mesh_w;
  Teuchos::RCP<InitialCondition> m_ic;
  Teuchos::RCP<GlobalGather_ReactiveFlow> m_gather;
  
};

}
// end namespace valiant

#endif /** __Valiant_InitialSolution_ReactiveFlow_HPP__ */
