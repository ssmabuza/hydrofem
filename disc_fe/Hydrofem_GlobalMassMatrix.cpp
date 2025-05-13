// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LocalMassMatrix.hpp"
#include "Hydrofem_GlobalMassMatrix.hpp"

namespace hydrofem
{

void GlobalMassMatrix(const std::shared_ptr<FEMatrix>& mass_mat, 
                      const std::shared_ptr<const DofMapper>& dofmapper,
                      const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                      const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature)
{
  mass_mat->setZero();
  const auto& loc_ind = dofmapper->getLocDofIndexes();
  auto local_mass_matrix = createKArray<LMAT_<double>>(dofmapper->local_ndof(),dofmapper->local_ndof());
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(dofmapper->nelements() * dofmapper->local_ndof() * dofmapper->local_ndof());
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    LocalMassMatrix(local_mass_matrix,dofmapper->getLocDofIndexes(),dofmapper->mesh()->getElementVertices(elem_ind),quadrature->at(elem_ind),*fe_basis);
    const auto& glob_ind = dofmapper->getGlobDofIndexes(elem_ind);
    for (int i = 0; i < dofmapper->local_ndof(); ++i)
      for (int j = 0; j < dofmapper->local_ndof(); ++j)
      {
//        tripletList.emplace_back(Eigen::Triplet<double>(glob_ind[i], glob_ind[j], local_mass_matrix(loc_ind[i], loc_ind[j])));
        tripletList.emplace_back(glob_ind[i], glob_ind[j], local_mass_matrix(loc_ind[i], loc_ind[j]));
      }
  }
  mass_mat->setFromTriplets(tripletList.begin(), tripletList.end());
}

void GlobalLumpedMassMatrix(const std::shared_ptr<FEVector>& lumped_mass_mat,
                            const std::shared_ptr<const DofMapper>& dofmapper,
                            const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                            const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature)
{
  lumped_mass_mat->setZero();
  auto& lumped_mass_mat_view = *lumped_mass_mat;
  const auto& loc_ind = dofmapper->getLocDofIndexes();
  LMAT_<RealType> local_mass_matrix = createKArray<LMAT_<RealType>>(dofmapper->local_ndof(),dofmapper->local_ndof());
  for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    LocalLumpedMassMatrix(local_mass_matrix,
                          dofmapper->getLocDofIndexes(),
                          dofmapper->mesh()->getElementVertices(elem_ind),
                          quadrature->at(elem_ind),
                          *fe_basis);
    const auto& glob_ind = dofmapper->getGlobDofIndexes(elem_ind);
    for (int i = 0; i < dofmapper->local_ndof(); ++i)
      lumped_mass_mat_view[glob_ind[i]] += local_mass_matrix(loc_ind[i],loc_ind[i]);
  }
}

}
// end namespace hydrofem
