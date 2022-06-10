// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IOTri.hpp"
#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_String_Utilities.hpp"


namespace hydrofem
{

void IOTri::writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const
{
  // print out the solution
  std::string file_name = m_field_name + intToStrSixDigits(m_counter) + ".vtk";
  // update the pvu file
  updatePVUOutputFile(file_name,time);
  // Start by printing the mesh
  printVTK(file_name, m_directory_name);
  // ------------------------------
  std::string out_filename = m_directory_name + "/" + file_name;
  std::ofstream outw(out_filename, std::ofstream::app);
  if (!outw.good())
  {
    std::stringstream ss;
    ss << "Error in IOTtri::writeSolutionToFile" << std::endl
       << "Can not open output file '" << out_filename << "'.";
    throw std::runtime_error(ss.str());
  }

  outw.precision(18);
  const int m_p = m_dofmapper->p();
  outw << "POINT_DATA " << m_dofmapper->nelements() * int(m_p + 1) * int(m_p + 2) / int(2) << std::endl;
  outw << "SCALARS " << m_field_name << " FLOAT" << std::endl;
  outw << "LOOKUP_TABLE default" << std::endl;

  for (int elem_ind(0); elem_ind < m_dofmapper->nelements(); ++elem_ind)
  {
    // gather the solution into the shape function
    LVEC_<double> coeff(sol.dimension(1));
    for (Eigen::Index i(0); i < coeff.size(); ++i) coeff(i) = sol(elem_ind,i);
    m_fe_shape->set_coefficients(coeff);
    const auto nodes = m_dofmapper->mesh()->getElementVertices(elem_ind);
    for (int i(0); i <= int(m_p); ++i)
    {
      for (int k(m_p - i), j(0); j <= int(m_p - i); ++j, --k)
      {
        SPoint x((nodes[0].x() * double(i) + nodes[1].x() * double(j) + nodes[2].x() * double(k)) / double(m_p),
                 (nodes[0].y() * double(i) + nodes[1].y() * double(j) + nodes[2].y() * double(k)) / double(m_p));
        outw << (*m_fe_shape)(x,nodes) << std::endl;
      }
    }
  }
  outw << std::endl;
  
  // at the end increase counter for the next printing
  m_counter++;
}

void IOTri::printVTK(const std::string& filename,
                     const std::string& directory) const
{
  // get the triangle dofmapper
  auto dofmapper_tri = std::dynamic_pointer_cast<HGrad_DofMapper_Triangle>(m_dofmapper);
  assert(dofmapper_tri);
  
  std::string out_filename = directory + "/" + filename;

  std::ofstream outw(out_filename);
  if (!outw.good())
  {
    std::stringstream ss;
    ss << "Error in \"printVTK\"." << std::endl
       << "Can not open output file '" << out_filename << "'.";
    throw std::runtime_error(ss.str());
  }
  outw.precision(18);

  // ------------------------------
  outw << "# vtk DataFile Version 3.1" << std::endl;
  outw << "Bernstein polynomial output by Sibusiso Mabuza" << std::endl;
  outw << "ASCII" << std::endl;
  outw << "DATASET UNSTRUCTURED_GRID" << std::endl;

  // ------------------------------
  const int m_p = dofmapper_tri->p();
  outw << "POINTS " << dofmapper_tri->nelements() * int(m_p + 1) * int(m_p + 2) / int(2) << " FLOAT" << std::endl;
  for (int elem_ind(0); elem_ind < dofmapper_tri->nelements(); ++elem_ind)
  {
    auto nodes = dofmapper_tri->mesh()->getElementVertices(elem_ind);
    for (int i(0); i <= m_p; ++i)
    {
      for (int k(m_p - i), j(0); j <= m_p - i; ++j, --k)
      {
        outw << (nodes[0].x() * double(i) + nodes[1].x() * double(j) + nodes[2].x() * double(k)) / double(m_p) << " "
             << (nodes[0].y() * double(i) + nodes[1].y() * double(j) + nodes[2].y() * double(k)) / double(m_p) << " "
             << double(0) << std::endl;
      }
    }
  }
  outw << std::endl;

  // ------------------------------
  outw << "CELLS " << dofmapper_tri->nelements() * int(m_p * m_p)
       << " " << dofmapper_tri->nelements() * int(m_p * m_p) * int(4) << std::endl;

  const auto& loc_indexes = dofmapper_tri->getLocDofIndexes();
  
  for (int elem_ind(0); elem_ind < dofmapper_tri->nelements(); ++elem_ind)
  {
    int tmp(elem_ind * static_cast<int>(loc_indexes.size()));
    for (int i(0); i < int(m_p); ++i)
    {
      for (int j(0); j < int(m_p - i); ++j)
      {
        outw << "3 "
             << tmp + dofmapper_tri->local(i,j) << " "
             << tmp + dofmapper_tri->local(i,j+1) << " "
             << tmp + dofmapper_tri->local(i+1,j) << std::endl;
        if (j < int(m_p - i - 1))
        {
          outw << "3 "
               << tmp + dofmapper_tri->local(i+1,j) << " "
               << tmp + dofmapper_tri->local(i,j+1) << " "
               << tmp + dofmapper_tri->local(i+1,j+1) << std::endl;
        }
      }
    }
  }
  outw << std::endl;

  outw << "CELL_TYPES " << dofmapper_tri->nelements() * int(m_p * m_p) << std::endl;
  for (int i(0); i < dofmapper_tri->nelements() * int(m_p * m_p); ++i)
  {
    outw << "5 ";
  }
  outw << std::endl;
  outw << std::endl;
}

}
// end namespace hydrofem
