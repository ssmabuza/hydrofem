// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_LocalConvectionMatrix.hpp"
#include "Hydrofem_GlobalConvectionMatrix.hpp"

namespace hydrofem
{

void GlobalConvectionMatrices(const std::vector<std::shared_ptr<FEMatrix>>& conv_mats,
                              const std::shared_ptr<const DofMapper>& dofmapper,
                              const std::shared_ptr<std::vector<std::shared_ptr<FEBasis>>>& fe_basis,
                              const std::shared_ptr<std::vector<std::shared_ptr<Quadrature>>>& quadrature)
{
  // make sure the dimensions match
  eigen_assert(static_cast<int>(conv_mats.size())==dofmapper->mesh()->numOfDims());
  
  for (std::size_t dim = 0; dim < conv_mats.size(); ++dim)
    conv_mats.at(dim)->setZero();
  const auto& loc_ind = dofmapper->getLocDofIndexes();
  LMAT_<RealType> local_conv_matrix = createKArray<LMAT_<RealType>>(dofmapper->local_ndof(),dofmapper->local_ndof());
  LVEC_<int> indices(dofmapper->local_ndof());
  LVEC_<double> values(dofmapper->local_ndof());
  std::vector<Eigen::Triplet<double>> tripletList;
  for (std::size_t dim = 0; dim < conv_mats.size(); ++dim)
  {
    tripletList.clear();
    tripletList.reserve(dofmapper->nelements() * dofmapper->local_ndof() * dofmapper->local_ndof());
    for (int elem_ind = 0; elem_ind < dofmapper->nelements(); ++elem_ind)
    {
      LocalConvectionMatrix(local_conv_matrix,
                            dofmapper->getLocDofIndexes(),
                            dofmapper->mesh()->getElementVertices(elem_ind),
                            quadrature->at(elem_ind),
                            *fe_basis,
                            [&](SPoint)->SPoint { auto val = SPoint(static_cast<int>(conv_mats.size())); val.data()(dim) = 1.0; return val; });
      const auto& glob_ind = dofmapper->getGlobDofIndexes(elem_ind);
      for (int i = 0; i < dofmapper->local_ndof(); ++i)
      {
        for (int j = 0; j < dofmapper->local_ndof(); ++j)
        {
          tripletList.emplace_back(glob_ind[i], glob_ind[j], local_conv_matrix(loc_ind[i], loc_ind[j]));
        }
      }
    }
    conv_mats.at(dim)->setFromTriplets(tripletList.begin(), tripletList.end());
  }
}

}
// end namespace hydrofem
