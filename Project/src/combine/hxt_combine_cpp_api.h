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

#ifndef HXT_COMBINE_CPP_API_
#define HXT_COMBINE_CPP_API_

#include <vector>
#include <map>

#include <hxt_combine_cell.h>
#include <hxt_graph.h>
#include <hxt_api.h>
#include <hxt_mwis.h>
#include <algorithms.h>
#include <tet_mesh.h>

/**
* \file tet_mesh.h C++ API of the Combine module
* \author Jeanne Pellerin
*/

namespace HXTCombine {

/**
* \addtogroup Combine
* @{
*/


/**
* API structure to create, store and select potential cells 
*
* \warning Implementation of member functions is distributed in several files
*/
class HXTCombineCellStore {
public:
  enum CELLTYPE {HEX, PRISM, PYRAMID, TET};
  static const unsigned int nbCellTypes = 4;

  HXTCombineCellStore(const TetMeshForCombining& mesh):
    mesh_(mesh)
  {
    cells_.resize(nbCellTypes-1);
    selectedCells_.resize(nbCellTypes);
  }

double computeHexesVolume();

  void computeHexes(double minQuality);
  void computePrisms(double minQuality);
  void computePyramids(double minQuality);

  void selectCellsGreedy(std::array<bool,nbCellTypes> cellType);
  void selectCellsGreedyLocal(std::array<bool,nbCellTypes> cellType);
  void selectCellsGraph (std::array<bool,nbCellTypes> cellType,  const TetMeshWrapper& tets);
  
  //void computeIncompatibilityGraph(HXTGraph& graph, std::array<bool, nbCellTypes> cellType);
    void computeIncompatibilityGraph(HXTGraph& graph, std::array<bool, nbCellTypes> cellType,const TetMeshWrapper& tets);

  unsigned int nbSelectedHexes() const;
  unsigned int nbSelectedPrisms() const;
  unsigned int nbSelectedPyramids() const;
  unsigned int nbSelectedTets() const;

  void saveMSH (const std::string& filename, std::array<bool,nbCellTypes> cellType);
  void saveMESH(const std::string& filename, std::array<bool,nbCellTypes> cellType);

public:
  std::vector<HXTCombineCell>& hexes() { return cells_[HEX]; }
  std::vector<HXTCombineCell>& prisms() { return cells_[PRISM]; }
  std::vector<HXTCombineCell>& pyramids() { return cells_[PYRAMID]; }

  std::vector<bool>& selectedHexes() { return selectedCells_[HEX]; }
  std::vector<bool>& selectedPrisms() { return selectedCells_[PRISM]; }
  std::vector<bool>& selectedPyramids() { return selectedCells_[PYRAMID]; }
  std::vector<bool>& selectedTets() { return selectedCells_[TET]; }

private:
  const TetMeshForCombining& mesh_;
  std::vector<std::vector<HXTCombineCell>> cells_;
  std::vector<std::vector<bool>> selectedCells_;
};


/* Functions maintened for C API  */ 
void computePotentialHexes(
  TetMeshForCombining& tetMesh, double hexQualityThreshold, vector<HXTCombineCell>& hexes);

void computePotentialPrisms(
  TetMeshForCombining& tetMesh, double hexQualityThreshold, vector<HXTCombineCell>& prisms);

/** 
* Compute the quality for all provided cells 
* \pre The array holding qualities is allocated and owned by the client
*/
void computeCellQualityVector(
  const TetMeshWrapper& tets, 
  const std::vector<HXTCombineCell>& cells,
  double* qualities);

double cellQualityAPI(const TetMeshWrapper& tets, const HXTCombineCell& cell);

void countAndCompareCellsIMR(std::ofstream& out, TetMeshForCombining& tetMesh, const vector<double>& minQualities);

void incompatibilityGraph(
  const TetMeshWrapper& mesh, const std::vector<HXTCombineCell>& cells, HXTGraph& graph);


/**
 * 通过cell顶点规范化的顺序鉴定是否cell唯一
* \brief Returns the \param indices of unique instance of HXTCandidate cells
* \details Duplicates are detected using a normalized order on the cells vertices.
* The output size of indices gives the number of different cells in the \param cells vector.
*/
void uniqueCells(const vector<HXTCombineCell>& cells, vector<unsigned int>& indices);

void writeCandidateCellMESH(
  const HXTCombineCell& cell,
  const TetMeshWrapper& tetMesh,
  std::ofstream& out,
  bool withBoundaryTets);


} // namespace

#endif
