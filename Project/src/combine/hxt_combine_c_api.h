/*
 * HXT - Copyright (C) <2016-2018> <UniversitÃ© catholique de Louvain (UCL), Belgique>
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

#ifndef HXT_COMBINE_C_API_
#define HXT_COMBINE_C_API_


#ifdef __cplusplus
extern "C" {
#endif

#include <hxt_api.h>

/**
* \file hxt_combine_c_api.h C interface of the Combine module
* \author Jeanne Pellerin
*/

/** 
* Opaque structure storing an array of HXTCombineCells
*/
struct HXTCombineCellVector;
typedef HXTIndex *HXTGraph;


HXTStatus HXTCombineCellVectorSize(HXTCombineCellVector *vector, HXTIndex *size);

HXTStatus HXTCombineCellVectorCellVertices(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **vertices, HXTIndex* numVertices);

HXTStatus HXTCombineCellVectorCellInteriorTets(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **tets, HXTIndex* nbTets);

HXTStatus HXTCombineCellVectorCellBoundaryTets(
  HXTCombineCellVector *vector, HXTIndex cellIndex, HXTIndex **tets, HXTIndex* nbTets);

/**
* \defgroup Combine Combination of tetrahedra into hexahedra and prisms
*
* Algorithms to combine tetrahedra in hexahedra and prisms described in the article :
* Pellerin J, Johnen A, Remacle J-F. 
* Identifying combinations of tetrahedra into hexahedra: a vertex based strategy. Procedia Engineering. 2017;203:2–13.
* Proceedings of the 26th International Meshing Roundtable, Barcelona, Spain  [arXiv:1705.02451]
*/

/**
* \addtogroup Combine
* @{
*/

/**
* Identifies all possible combinations of tetrahedra into hexahedra in a given mesh
* that are valid and which quality is above the given theshold.
* 
* The memory for hexahedra is allocated by the function.
* Ownership is transfered to the client.
*/
HXTStatus hxtIdentifyAllHexahedra(
  const HXTMesh* mesh, double qualityThreshold, HXTCombineCellVector* hexahedra);

HXTStatus hxtIdentifyAllPrisms(
  const HXTMesh* mesh, double qualityThreshold, HXTCombineCellVector* prisms);

HXTStatus hxtComputeCellQuality(
  const HXTMesh* mesh, const HXTCombineCellVector *vector, double** qualities);


/**
* \brief Build the incompatibility graph for the input cells
* \details The \param cells are defined as combinations of tets of the \param mesh
*/
HXTStatus hxtBuildIncompatibilityGraph(
  HXTMesh* mesh, HXTCombineCellVector* hexahedra, HXTGraph* graph);


#ifdef __cplusplus
}
#endif

/**
* @}
*/

#endif
