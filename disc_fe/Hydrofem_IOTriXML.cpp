// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright (2016) S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Hydrofem_FEUtils.hpp"
#include "Hydrofem_IOTriXML.hpp"
#include "Hydrofem_LocalArray.hpp"
#include "Hydrofem_String_Utilities.hpp"


namespace hydrofem
{

void IOTriXML::writeSolutionTofile(const std::shared_ptr<FEVector>& sol, const double time) const
{

  // get the triangle dofmapper
  auto dofmapper_tri = std::dynamic_pointer_cast<HGrad_DofMapper_Triangle>(m_dofmapper);
  assert(dofmapper_tri);
  const int m_p = m_dofmapper->p();

  const auto& sol_view = * sol;
  auto& out_u_view = * m_out_vals;
  auto& out_points_x = * (m_out_points[0]);
  auto& out_points_y = * (m_out_points[1]);

  // compute the output values that are needed
  {
    auto coeffs = createKArray<LVEC_<double>>(m_dofmapper->local_ndof());
    for (int elem_ind(0); elem_ind < m_dofmapper->nelements(); ++elem_ind)
    {
      // get the output points for high order
      const auto nodes = m_dofmapper->mesh()->getElementVertices(elem_ind);
      if ((0==m_counter) || m_fsi)
      {
        for (int i(0); i <= int(m_p); ++i)
        {
          for (int k(m_p - i), j(0); j <= int(m_p - i); ++j, --k)
          {
            const int glob_ind = dofmapper_tri->global(elem_ind, dofmapper_tri->local(i, j, k));
            out_points_x[glob_ind] = (nodes[0].x() * double(i) + nodes[1].x() * double(j) + nodes[2].x() * double(k)) / double(m_p);
            out_points_y[glob_ind] = (nodes[0].y() * double(i) + nodes[1].y() * double(j) + nodes[2].y() * double(k)) / double(m_p);
          }
        }
      }

      // build the local coeffs and compute output values
      const auto& local = m_dofmapper->getLocDofIndexes();
      const auto& global = m_dofmapper->getGlobDofIndexes(elem_ind);
      for (std::size_t i (0); i < local.size(); ++i)
        coeffs[local[i]] = sol_view[global[i]];
      m_fe_shape->set_coefficients(coeffs);
      for (std::size_t i (0); i < local.size(); ++i)
      {
        SPoint x (out_points_x[global[i]],out_points_y[global[i]]);
        out_u_view[global[i]] = (*m_fe_shape)(x,nodes);
      }
    }
  }
  
  // initially create the pvu file
  if (0==m_counter) createPVUOutputFile();
  // print out the solution
  std::string file_name = m_field_name + intToStrSixDigits(m_counter) + ".vtu";
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
    ss << "Error in \"writeSolutionTofile\"." << std::endl
           << "Can not open output file '" << out_filename << "'.";
    throw std::runtime_error(ss.str());
  }
  
  outw.precision(18);
  outw << "<PointData  Scalars=\"" << m_field_name << "\">" << std::endl;
  outw << R"(<DataArray  type="Float64"  Name=")" << m_field_name << R"("  format="ascii">)" << std::endl;
  for (Eigen::Index i = 0; i < out_u_view.size(); ++i)
    outw << out_u_view[i] << " " << std::endl;
  outw << "</DataArray>" << std::endl;
  outw << "</PointData>" << std::endl;
  outw << "</Piece>" << std::endl;
  outw << "</UnstructuredGrid>" << std::endl;
  outw << "</VTKFile>" << std::endl;
  
  // at the end finalize the pvu file otherwise increase counter for the next printing
  if ((time >= m_time_final) || m_steady)
  {
    finalizePVUOutputFile();
    return;
  } else {
    m_counter++;
  }
  // catch it regardless
  finalizePVUOutputFile();
  
}

void IOTriXML::printVTK(const std::string& filename,
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
    ss << "Error in \"printVTKXML\"." << std::endl
          << "Can not open output file '" << out_filename << "'.";
    throw std::runtime_error(ss.str());
  }
  
  const auto& out_points_x = *(m_out_points[0]);
  const auto& out_points_y = *(m_out_points[1]);
  
  outw.precision(18);
  outw << "<?xml version=\"1.0\"?>" << std::endl;
  outw << R"(<VTKFile type="UnstructuredGrid"  version="0.1"  >)" << std::endl;
  outw << "<UnstructuredGrid>" << std::endl;
  outw << "<Piece NumberOfPoints=\"" << dofmapper_tri->global_ndof() << "\" NumberOfCells=\"" << dofmapper_tri->nelements() << "\">" << std::endl;
  outw << "<Points>" << std::endl;
  outw << R"(<DataArray  type="Float64"  NumberOfComponents="3"  format="ascii">)" << std::endl;
  // ------------------------------
  for (int point_ind = 0; point_ind < out_points_x.size(); ++point_ind)
  {
    outw << out_points_x[point_ind] << " " << out_points_y[point_ind] << " " << 0.0 << std::endl;
  }
  const int m_p = dofmapper_tri->p();
  outw << "</DataArray>" << std::endl;
  outw << "</Points>" << std::endl;
  // ---------------------------------
  
  // ------------------------------
  outw << "<Cells>" << std::endl;
  outw << R"(<DataArray  type="UInt32"  Name="connectivity"  format="ascii">)" << std::endl;
  for (int elem_ind(0); elem_ind < dofmapper_tri->nelements(); ++elem_ind)
  {
    for (int i(0); i < int(m_p); ++i)
    {
      for (int j(0); j < int(m_p - i); ++j)
      {
        outw << dofmapper_tri->eval_global(elem_ind,i,j)   << " "
             << dofmapper_tri->eval_global(elem_ind,i,j+1) << " "
             << dofmapper_tri->eval_global(elem_ind,i+1,j) << std::endl;
        if (j < int(m_p - i - 1))
        {
          outw << dofmapper_tri->eval_global(elem_ind,i+1,j)   << " "
               << dofmapper_tri->eval_global(elem_ind,i,j+1)   << " "
               << dofmapper_tri->eval_global(elem_ind,i+1,j+1) << std::endl;
        }
      }
    }
  }
  outw << "</DataArray>" << std::endl;
  outw << R"(<DataArray  type="UInt32"  Name="offsets"  format="ascii">)" << std::endl;
  for (int i(0); i < dofmapper_tri->nelements() * int(m_p * m_p); ++i)
    outw << (i+1)*3 << " ";
  outw << "</DataArray>" << std::endl;
  outw << R"(<DataArray  type="UInt8"  Name="types"  format="ascii">)" << std::endl;
  for (int i(0); i < dofmapper_tri->nelements() * int(m_p * m_p); ++i)
    outw << "5 ";
  outw << "</DataArray>" << std::endl;
  outw << "</Cells>" << std::endl;
  // ---------------------------------
  
}

}
// end namespace hydrofem
