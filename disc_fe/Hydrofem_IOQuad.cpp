// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IOQuad.hpp"
#include "Hydrofem_String_Utilities.hpp"
#include "Hydrofem_ReferenceQuadrilateral.hpp"

namespace hydrofem
{

IOQuad::IOQuad(const std::shared_ptr<HGrad_DofMapper_Quadrilateral>& dofmapper,
               const std::shared_ptr<FunctionElement<double>>& fe_shape)
  :
  IOBase(dofmapper,fe_shape)
{
  int m_p = dofmapper->p();
  std::vector<SPoint> ref_points;
  double hx = 1.0/m_p, hy = 1.0/m_p;
  for (int i(0); i <= int(m_p); ++i)
    for (int j(0); j <= int(m_p); ++j)
      ref_points.emplace_back(i*hx,j*hy);
  m_out_points.resize(dofmapper->nelements(),std::vector<SPoint>(ref_points.size()));
  for (int elem_ind(0); elem_ind < dofmapper->nelements(); ++elem_ind)
  {
    auto nodes = dofmapper->mesh()->getElementVertices(elem_ind);
    // map the ref point to this cell
    for (std::size_t pt (0); pt < ref_points.size(); ++pt)
      m_out_points[elem_ind][pt] = referenceQuadToPhysical(nodes[0],nodes[1],nodes[2],nodes[3],ref_points[pt]);
  }
}


void IOQuad::writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const
{
  // initially create the pvu file
  if (0==m_counter) createPVUOutputFile();
  // print out the solution
  std::string file_name = m_field_name + intToStrSixDigits(m_counter) + ".vtu";
  // update the pvu file
  updatePVUOutputFile(file_name,time);
  // get the quad dofmapper
  auto dofmapper_quad = std::dynamic_pointer_cast<HGrad_DofMapper_Quadrilateral>(m_dofmapper);
  assert(dofmapper_quad);
  
  std::string out_filename = m_directory_name + "/" + file_name;

  std::ofstream outw(out_filename);
  if (!outw.good())
  {
    std::stringstream ss;
    ss << "Error in \"printVTK\"." << std::endl
       << "Can not open output file '" << out_filename << "'." << std::endl;
    throw std::runtime_error(ss.str());
  }
  outw.precision(18);

  // ------------------------------
  outw << "# vtk DataFile Version 3.1" << std::endl;
  outw << "Bernstein polynomial output by Sibusiso Mabuza" << std::endl;
  outw << "ASCII" << std::endl;
  outw << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // ------------------------------
  const int m_p = dofmapper_quad->p();
  outw << "POINTS " << dofmapper_quad->nelements() * int(m_p + 1) * int(m_p + 1) << " FLOAT" << std::endl;
  
  for (int elem_ind(0); elem_ind < dofmapper_quad->nelements(); ++elem_ind)
    // map the ref point to this cell
    for (const auto & pt : m_out_points[elem_ind])
      outw << pt.x() << " " << pt.y() << " " << 0.0 << std::endl;
  outw << std::endl;

  // ------------------------------
  outw << "CELLS " << dofmapper_quad->nelements() * int(m_p * m_p)
       << " " << dofmapper_quad->nelements() * int(m_p * m_p) * int(5) << std::endl;

  const auto& loc_indexes = dofmapper_quad->getLocDofIndexes();
  
  for (int elem_ind(0); elem_ind < dofmapper_quad->nelements(); ++elem_ind)
  {
    int tmp(elem_ind * static_cast<int>(loc_indexes.size()));
    for (int i(0); i < int(m_p); ++i)
    {
      for (int j(0); j < int(m_p); ++j)
      {
        outw << "4 "
             << tmp + dofmapper_quad->local(i,j) << " "
             << tmp + dofmapper_quad->local(i+1,j) << " "
             << tmp + dofmapper_quad->local(i+1,j+1) << " "
             << tmp + dofmapper_quad->local(i,j+1) << std::endl;
      }
    }
  }
  outw << std::endl;

  outw << "CELL_TYPES " << dofmapper_quad->nelements() * int(m_p * m_p) << std::endl;
  for (int i(0); i < dofmapper_quad->nelements() * int(m_p * m_p); ++i)
  {
    outw << "9 ";
  }
  outw << std::endl;
  outw << std::endl;

  outw << "POINT_DATA " << dofmapper_quad->nelements() * loc_indexes.size() << std::endl;
  outw << "SCALARS " << m_field_name << " FLOAT" << std::endl;
  outw << "LOOKUP_TABLE default" << std::endl;

  for (int elem_ind(0); elem_ind < dofmapper_quad->nelements(); ++elem_ind)
  {
    // gather the solution into the shape function
    LVEC_<double> coeff(sol.dimension(1));
    for (Eigen::Index i(0); i < coeff.size(); ++i) coeff(i) = sol(elem_ind,i);
    m_fe_shape->set_coefficients(coeff);
    const auto nodes = dofmapper_quad->mesh()->getElementVertices(elem_ind);
    for (const auto & pt : m_out_points[elem_ind])
      outw << (*m_fe_shape)(pt,nodes) << std::endl;
  }
  outw << std::endl;
  
  if (time >= m_time_final)
    finalizePVUOutputFile();
  else
    m_counter++;
  
}

}
// end namespace hydrofem
