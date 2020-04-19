/*
 * HXT - Copyright (C) <2016-2018> <Université catholique de Louvain (UCL), Belgique>
 *
 * List of the contributors to the development of HXT: see AUTHORS file.
 * Description and complete License: see LICENSE file.
 * 	
 * This program (HXT) is free software: 
 * you can redistribute it and/or modify it under the terms of the GNU Lesser
 * General Public License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program (see COPYING and COPYING.LESSER files).  If not, 
 * see <http://www.gnu.org/licenses/>.
 */

#include "hxt_combine_c_api.h"

#include <vector>
#include "hxt_combine_cell.h"
#include "tet_mesh.h"

#include "hxt_combine_cpp_api.h"


/**
* \author Jeanne Pellerin
*/

using HXTCombine::HXTCombineCell;

struct HXTCombineCellVector {
  std::vector<HXTCombineCell> data;
};


HXTStatus HXTCombineCellVectorSize(HXTCombineCellVector *vector, HXTIndex *size)
{
  *size = vector->data.size();
  return HXT_STATUS_OK;
}

HXTStatus HXTCombineCellVectorCellVertices(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **vertices, HXTIndex* numVertices) 
{
  HXTCombineCell& cell = vector->data[cellIndex];
  *numVertices = cell.nbVertices();
  *vertices = &cell.vertexes[0];
  return HXT_STATUS_OK;
}

HXTStatus HXTCombineCellVectorCellInteriorTets(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **tets, HXTIndex* nbTets)
{
  HXTCombineCell& cell = vector->data[cellIndex];
  *nbTets = cell.nbInteriorTets();
  *tets = &cell.interiorTetrahedra[0];
  return HXT_STATUS_OK;
}

HXTStatus HXTCombineCellVectorCellBoundaryTets(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **tets, HXTIndex* nbTets)
{
  HXTCombineCell& cell = vector->data[cellIndex];
  *nbTets = cell.nbBoundaryTets();
  *tets = &cell.boundaryTetrahedra[0];
  return HXT_STATUS_OK;
}


HXTStatus hxtIdentifyAllHexahedra(
  const HXTMesh* mesh, double qualityThreshold, HXTCombineCellVector* hexahedra)
{
  HXTCombine::TetMeshForCombining fullMesh(mesh);
  HXTCombine::computePotentialHexes(fullMesh, qualityThreshold, hexahedra->data);
  return HXT_STATUS_OK;
}

HXTStatus hxtIdentifyAllPrisms(
  const HXTMesh* mesh, double qualityThreshold, HXTCombineCellVector* prisms)
{
  HXTCombine::TetMeshForCombining fullMesh(mesh);
  HXTCombine::computePotentialPrisms(fullMesh, qualityThreshold, prisms->data);
  return HXT_STATUS_OK;
}

HXTStatus hxtComputeCellQuality(
  const HXTMesh* mesh, const HXTCombineCellVector *vector, double** qualities)
{
  HXT_CHECK(hxtMalloc((void**)qualities, sizeof(**qualities)*vector->data.size()));
  HXTCombine::TetMeshWrapper wrapped(mesh);
  HXTCombine::computeCellQualityVector(wrapped, vector->data, *qualities);
  return HXT_STATUS_OK;
}

/**
* \brief Build the incompatibility graph for the input cells
* \details The \param cells are defined as combinations of tets of the \param mesh
*/
//HXTStatus hxtBuildIncompatibilityGraph(
//  HXTMesh* mesh, 
//  HXTCombineCellVector* hexahedra,
//  HXTGraph* graph)
//{
//  // TODO we do not need to build these -- the graph construction should work 
//  // on any not too stupid mesh structure
//  HXTCombine::TetMeshForCombining fullMesh(mesh);
//  
////  HXTCombine::incompatibilityGraph(fullMesh, hexahedra->data, *graph);
//
//  return HXT_STATUS_OK;
//}

