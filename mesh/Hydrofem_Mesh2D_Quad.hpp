// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef __Hydrofem_Mesh2D_Quad_HPP__
#define __Hydrofem_Mesh2D_Quad_HPP__

#include "Hydrofem_Mesh2D.hpp"

namespace hydrofem
{

class Mesh2D_Quad
  :
  public Mesh2D
{
public:
  
  using mesh_info_type = Mesh::mesh_info_type;

  Mesh2D_Quad() : Mesh2D() {}
  
  Mesh2D_Quad(const Mesh2D_Quad::mesh_info_type& mesh_data);
  
  Mesh2D_Quad(const std::string mesh_filename);
  
  Mesh2D_Quad(const int nx, const int ny, const double x0, const double xf, const double y0, const double yf);
  
  virtual ~Mesh2D_Quad() {}
  
  virtual std::shared_ptr<Mesh> refineMesh() override;
  
  virtual void updateNodes_AM_FMT(const std::string filename_am_fmt) override;

private:

  /// Mesh forming methods
  //@{
  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief function GenerateUniformRectangularMesh
  //! Input:
  //! \param nx --- number of mesh steps in x direction
  //! \param ny --- number of mesh steps in y direction
  //! \param L  --- length; H -- Height
  /////////////////////////////////////////////////////////////////////////////////////////
  void generateUniformRectangularMesh(const int nx,
                                      const int ny,
                                      const double L,
                                      const double H);

  /// Mesh forming methods
  //@{
  /////////////////////////////////////////////////////////////////////////////////////////
  //! \brief function GenerateUniformRectangularMesh
  //! Input:
  //! \param nx  --- number of mesh steps in x direction
  //! \param ny  --- number of mesh steps in y direction
  //! \param L_l --- left x value; 
  //! \param L_r --- right x value; 
  //! \param H_l --- lower y value
  //! \param H_r --- upper y value
  /////////////////////////////////////////////////////////////////////////////////////////
  inline void generateUniformRectangularMesh(const int nx,
                                             const int ny, 
                                             const double L_l,
                                             const double L_r, 
                                             const double H_l,
                                             const double H_r)
  {
    const double H = H_r - H_l;
    const double L = L_r - L_l;
    generateUniformRectangularMesh(nx,ny,L,H);
    
    for (int i = 0; i < numOfPoints(); ++i)
    {
      points().at(i)(0) += L_l;
      points().at(i)(1) += H_l;
    }
  }
  
  virtual void readMesh_AM_FMT(const std::string filename_am_fmt) override;
  
};

}
// end namespace hydrofem

#endif /** __Hydrofem_Mesh2D_Quad_HPP__ */
