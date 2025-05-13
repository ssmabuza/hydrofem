// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Mesh2D_Tri_HPP__
#define __Hydrofem_Mesh2D_Tri_HPP__


#include "Hydrofem_Mesh2D.hpp"

namespace hydrofem
{

class Mesh2D_Tri
  :
  public Mesh2D
{
public:

  using mesh_info_type = Mesh::mesh_info_type;

  Mesh2D_Tri() {}
  
  Mesh2D_Tri(const Mesh2D_Tri::mesh_info_type& mesh_data)
  : Mesh2D()
  {
    const std::size_t& ndims = std::get<0>(mesh_data);
    assert(ndims == 2);
    m_mesh_type = typeTriangularMesh;
    assert(m_mesh_type==std::get<3>(mesh_data));
    num_dims() = ndims;
    const std::vector<double>& sides = std::get<1>(mesh_data);
    const double& L_l = sides.at(0);
    const double& L_r = sides.at(1);
    const double& H_l = sides.at(2);
    const double& H_r = sides.at(3);
    const std::vector<std::size_t>& sizes = std::get<2>(mesh_data);
    const std::size_t& nx = sizes.at(0);
    const std::size_t& ny = sizes.at(1);
    const TriangulationType& orientation = std::get<4>(mesh_data);
    generateTriangulatedUniformRectangularMesh(nx,ny,L_l,L_r,H_l,H_r,orientation);
  }  
  
  
  Mesh2D_Tri(const std::string mesh_filename);
  
  Mesh2D_Tri(const int nx, const int ny, const double x0, const double xf, const double y0, const double yf, const TriangulationType orientation)
  {
    num_dims() = 2;
    m_mesh_type = typeTriangularMesh;
    generateTriangulatedUniformRectangularMesh(nx,ny,x0,xf,y0,yf,orientation);
  }
  
  virtual ~Mesh2D_Tri() {}
  
  virtual std::shared_ptr<Mesh> refineMesh() override;

private:

  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief function GenerateTriangulatedUniformRectangularMesh
  //! Input:
  //! \param nx --- number of mesh steps in x direction
  //! \param ny --- number of mesh steps in y direction
  //! \param L  --- length
  //! \param H  --- Height
  //! \param orientation --- triangulation type
  /////////////////////////////////////////////////////////////////////////////////////////
  void generateTriangulatedUniformRectangularMesh(const int nx,
                                                  const int ny,
                                                  const double L, 
                                                  const double H,
                                                  const TriangulationType orientation);

  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief function GenerateTriangulatedUniformRectangularMesh
  //! Input:
  //! \param nx  --- number of mesh steps in x direction
  //! \param ny  --- number of mesh steps in y direction
  //! \param L_l --- left x value; 
  //! \param L_r --- right x value; 
  //! \param H_l --- lower y value
  //! \param H_r --- upper y value
  //! \param orientation --- triangulation type
  /////////////////////////////////////////////////////////////////////////////////////////
  inline
  void generateTriangulatedUniformRectangularMesh(const int nx,
                                                  const int ny, 
                                                  const double L_l,
                                                  const double L_r, 
                                                  const double H_l,
                                                  const double H_r,
                                                  const TriangulationType orientation)
  {
    const double H = H_r - H_l;
    const double L = L_r - L_l;
    generateTriangulatedUniformRectangularMesh(nx,ny,L,H,orientation);
    for (int i = 0; i < numOfPoints(); ++i)
      points().at(i) += SPoint(L_l,H_l);
  }

  virtual void readMesh_AM_FMT(const std::string filename_am_fmt) override;
    
  virtual void updateNodes_AM_FMT(const std::string filename_am_fmt) override;
  
//  virtual void setEdges();
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Mesh2D_Tri_HPP__ */

