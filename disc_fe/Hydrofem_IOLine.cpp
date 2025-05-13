// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_IOLine.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_FunctionElement.hpp"

#include "Hydrofem_String_Utilities.hpp"

namespace hydrofem
{

IOLine::IOLine(const std::shared_ptr<HGrad_DofMapper_Line>& dofmapper,
               const std::shared_ptr<FunctionElement<double>>& fe_shape)
  :
  IOBase(dofmapper,fe_shape)
{
}

void IOLine::writeSolutionTofile(const FEArray<double>::CellBasis& sol, const double time) const
{
  // initially create the pvu file
  if (0==m_counter) createPVUOutputFile();
  // print out the solution
  std::string file_name = m_field_name + intToStrSixDigits(m_counter)+".vtu";
  // update the pvu file
  updatePVUOutputFile(file_name,time);
  // get the line dofmapper 
  auto dofmapper_line = std::dynamic_pointer_cast<HGrad_DofMapper_Line>(m_dofmapper);
  assert(dofmapper_line);
  
  std::string out_filename = m_directory_name + "/" + file_name;

  std::ofstream outw(out_filename);
  if (!outw.good())
  {
    std::stringstream ss;
    ss << "Error in \"writeSolutionTofile\"." << std::endl
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
  const int m_p = dofmapper_line->p();
  outw << "POINTS " << dofmapper_line->nelements() * int(m_p + 1) << " FLOAT" << std::endl;
  
  for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
  {
    const auto nodes = dofmapper_line->mesh()->getElementVertices(elem_ind);
    const double hx = dofmapper_line->mesh()->getElementArea(elem_ind);
    // map the ref point to this cell
    for (int i (0); i <= m_p; ++i)
      outw << nodes[0].x() + i*hx << " " << 0.0 << " " << 0.0 << std::endl;
  }
  outw << std::endl;

  // ------------------------------
  outw << "CELLS " << dofmapper_line->nelements() * int(m_p)
       << " " << dofmapper_line->nelements() * int(m_p) * int(3) << std::endl;

  const auto& loc_indexes = dofmapper_line->getLocDofIndexes();
  
  for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
  {
    int tmp(elem_ind * loc_indexes.size());
    for (int i(0); i < int(m_p); ++i)
    {
      outw << "2 "
            << tmp + dofmapper_line->local(i) << " "
            << tmp + dofmapper_line->local(i+1) << " " << std::endl;
    }
  }
  outw << std::endl;

  outw << "CELL_TYPES " << dofmapper_line->nelements() * int(m_p) << std::endl;
  for (int i(0); i < dofmapper_line->nelements() * int(m_p); ++i)
  {
    outw << "3 ";
  }
  outw << std::endl;
  outw << std::endl;

  outw << "POINT_DATA " << dofmapper_line->nelements() * loc_indexes.size() << std::endl;
  outw << "SCALARS " << m_field_name << " FLOAT" << std::endl;
  outw << "LOOKUP_TABLE default" << std::endl;

  for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
  {
    // gather the solution into the shape function
    LVEC_<double> coeff(sol.dimension(1));
    for (Eigen::Index i(0); i < coeff.size(); ++i) coeff(i) = sol(elem_ind,i);
    m_fe_shape->set_coefficients(coeff);
    
    const auto nodes = dofmapper_line->mesh()->getElementVertices(elem_ind);
    const double hx = dofmapper_line->mesh()->getElementArea(elem_ind);
    // map the ref point to this cell
    for (int i (0); i <= m_p; ++i)
      outw << (*m_fe_shape)(SPoint(double(nodes[0].x() + i*hx)),nodes) << std::endl;
  }
  outw << std::endl;

  if (time >= m_time_final)
    finalizePVUOutputFile();
  else
    m_counter++;
}

// void IOLine::writeSolutionTofileXML(const FEArray<double>::CellBasis& sol, const double time) const
// {
//   // initially create the pvu file
//   if (0==m_counter) createPVUOutputFile();
//   // print out the solution
//   std::string file_name = m_field_name + intToStrSixDigits(m_counter)+".vtu";
//   // update the pvu file
//   updatePVUOutputFile(file_name,time);
//   // get the line dofmapper
//   auto dofmapper_line = Teuchos::rcp_dynamic_cast<HGrad_DofMapper_Line>(m_dofmapper);
//   TEUCHOS_ASSERT(not dofmapper_line.is_null());
//   
//   std::string out_filename = m_directory_name + "/" + file_name;
//   
//   std::ofstream outw(out_filename);
//   TEUCHOS_TEST_FOR_EXCEPTION(!outw.good(),std::runtime_error,"Error in \"writeSolutionTofileXML\"." << std::endl
//                                    << "Can not open output file '" << out_filename << "'.")
// 
//   outw.precision(18);
//   
//   // ------------------------------
//   outw << "<?xml version=\"1.0\"?>" << std::endl;
//   outw << "<VTKFile type=\"UnstructuredGrid\"  version=\"0.1\"  >" << std::endl;
//   outw << "<UnstructuredGrid>" << std::endl;
//   outw << "<Piece NumberOfPoints=\"" << dofmapper_line->global_ndof(0) << "\" NumberOfCells=\"" << dofmapper_line->nelements() << "\">" << std::endl;
//   
//   // ------------------------------
//   outw << "<Points>" << std::endl;
//   outw << "<DataArray  type=\"Float64\"  NumberOfComponents=\"3\"  format=\"ascii\">";
//   const int m_p = dofmapper_line->p();
//   for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
//   {
//     const auto nodes = dofmapper_line->mesh()->getElementVertices(elem_ind);
//     const double hx = dofmapper_line->mesh()->getElementArea(elem_ind);
//     // map the ref point to this cell
//     for (int i (0); i <= m_p; ++i)
//       outw << nodes[0].x() + i*hx << " " << 0.0 << " " << 0.0 << "  ";
//   }
//   outw << "</DataArray>" << std::endl;
//   outw << "</Points>" << std::endl;
//   // ---------------------------------
// 
//   
//   // ------------------------------
//   outw << "<Cells>" << std::endl;
//   outw << "<DataArray  type=\"UInt32\"  Name=\"connectivity\"  format=\"ascii\">";
//   const auto& loc_indexes = dofmapper_line->getLocDofIndexes();
//   for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
//   {
//     int tmp(elem_ind * loc_indexes.size());
//     for (int i(0); i < int(m_p); ++i)
//     {
//       outw << tmp + dofmapper_line->local(i) << " "
//            << tmp + dofmapper_line->local(i+1) << "  ";
//     }
//   }
//   outw << "</DataArray>" << std::endl;
//   outw << "<DataArray  type=\"UInt32\"  Name=\"offsets\"  format=\"ascii\">";
//   for (int i(0); i < dofmapper_line->nelements() * int(m_p); ++i)
//     outw << (i+1)*2 << " ";
//   outw << "</DataArray>" << std::endl;
//   outw << "<DataArray  type=\"UInt8\"  Name=\"types\"  format=\"ascii\">";
//   for (int i(0); i < dofmapper_line->nelements() * int(m_p); ++i)
//     outw << "3 ";
//   outw << "</DataArray>" << std::endl;
//   outw << "</Cells>" << std::endl;
//   // ---------------------------------
//   
//   outw << "<PointData  Scalars=\"" << m_field_name << "\">" << std::endl;
//   outw << "<DataArray  type=\"Float64\"  Name=\"" << m_field_name << "\"  format=\"ascii\">";
//   for (int elem_ind(0); elem_ind < dofmapper_line->nelements(); ++elem_ind)
//   {
//     // gather the solution into the shape function
//     auto coeff = Kokkos::subview(sol,elem_ind,Kokkos::ALL());
//     m_fe_shape->set_coefficients(coeff);
//     
//     const auto nodes = dofmapper_line->mesh()->getElementVertices(elem_ind);
//     const double hx = dofmapper_line->mesh()->getElementArea(elem_ind);
//     // map the ref point to this cell
//     for (int i (0); i <= m_p; ++i)
//       outw << (*m_fe_shape)(SPoint(double(nodes[0].x() + i*hx)),nodes) << "  ";
//   }
//   outw << "</DataArray>" << std::endl;
//   outw << "</PointData>" << std::endl;
//   outw << "</Piece>" << std::endl;
//   outw << "</UnstructuredGrid>" << std::endl;
//   outw << "</VTKFile>" << std::endl;
//   
//   if (time >= m_time_final)
//     finalizePVUOutputFile();
//   else
//     m_counter++;
// }


}
// end namespace hydrofem
