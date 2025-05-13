// @HEADER
// ****************************************************************************
//                Hydrofem: Copyright S. Mabuza
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER


#ifndef __Hydrofem_MeshTypes_HPP__
#define __Hydrofem_MeshTypes_HPP__

namespace hydrofem
{

enum ObjType
{
  typeGenericElement,
  typeLine,
  typeTriangle,
  typeQuadrilateral
};

enum MeshType
{
  typeGenericMesh,
  
  typeLineMesh,
  typeTriangularMesh,
  typeQuadrilateralMesh
};

enum TriangulationType
{
  typeGenericTriangulation,
  typeSWtoNE,
  typeSEtoNW,
  typeUnionJack,
  typeBottomNWtoSETopSEtoNW
};

enum MeshInput
{
  Inline,
  File
};

}
// end namespace hydrofem

#endif /** __Hydrofem_MeshTypes_HPP__ */
