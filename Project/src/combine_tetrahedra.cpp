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

#include <fstream>
#include <string>
#include <chrono>

//#include "hxt_option.h"
//#include "hxt_message.h"

#include "algorithms.h"
#include "tet_mesh.h"
#include "hxt_combine_cpp_api.h"
#include "hxt_combine_cell.h"
using namespace HXTCombine;
/**
* \file combine_tetrahedra.cpp Executable to generate a hex-dominant mesh from a tet mesh
* \author Jeanne Pellerin
*/
//int my_process(const MatrixXf &V_, const MatrixXu &F_, const MatrixXu &T_){
int my_process(const MultiResolutionHierarchy& mRes, HXTCombineCellStore* combineRes) {


  using namespace HXTCombine;
  char* outputFile = "H://xxx.msh";
  double minQuality = 0.3; // necessary as long as I do not have the pyramid quality

  int hexFlag = 1;
  //int prismFlag = 1;
  //int pyramidFlag = 1;
  int prismFlag = 0;
  int pyramidFlag = 0;

  const std::string outputFileStr(outputFile);

  // Using my mesh
  MeshStore ioMesh;
  //readFileMESH(inputFileStr, ioMesh);
  myReadFileMESH(mRes.mV[0], mRes.mF, mRes.mT, ioMesh);
  auto start0 = std::chrono::high_resolution_clock::now();
  
  TetMeshForCombining tets(&ioMesh);
  
  auto finish0 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> t_mesh(finish0 - start0);
  std::cout << "Mesh structure built in "<< t_mesh.count() << " seconds" << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  HXTCombineCellStore TheResult(tets);

  if( hexFlag) {
    TheResult.computeHexes(minQuality);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t0(finish - start);
    std::cout << TheResult.hexes().size() << " potential hexes computed in " << t0.count() << " seconds" << std::endl;
  }

  if( prismFlag ) {
    auto start = std::chrono::high_resolution_clock::now();
    TheResult.computePrisms(minQuality);
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> tPrism(finish - start);
    std::cout << TheResult.prisms().size() << " potential prisms computed in " << tPrism.count() << " seconds" << std::endl;
  }

  if (pyramidFlag) {
    auto start = std::chrono::high_resolution_clock::now();
    TheResult.computePyramids(minQuality);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tPyramid(finish - start);
    std::cout << TheResult.pyramids().size() << " potential pyramids computed in " << tPyramid.count() << " seconds" << std::endl;
  }

  auto startSelect = std::chrono::high_resolution_clock::now();

  std::array<bool, 4> cellTypes { bool(hexFlag), bool(prismFlag), bool(pyramidFlag), true };
  // tf_怎样优化候选cell的选择策略
  TheResult.selectCellsGreedy(cellTypes);
 // TheResult.selectCellsGreedyLocal(cellTypes);
  auto endSelect = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> ts(endSelect - startSelect);

  if (hexFlag) std::cout << nbTrueValues(TheResult.selectedHexes())    << " selected hexes" << std::endl;
  // 计算下六面体体积占比
  //double hexsVolume = TheResult.computeHexesVolume();
//  std::cout << " hexsVolume: "   << hexsVolume << std::endl;

  if (prismFlag)   std::cout << nbTrueValues(TheResult.selectedPrisms())   << " selected prisms" << std::endl;
  if (pyramidFlag) std::cout << nbTrueValues(TheResult.selectedPyramids()) << " selected pyramids" << std::endl;
  std::cout << nbTrueValues(TheResult.selectedTets()) << " tetrahedra remain" << std::endl;
  std::cout << "Timings cell selection "<<  ts.count() << "seconds" << std::endl;
  TheResult.saveMSH(outputFileStr, cellTypes);

  return HXT_STATUS_OK;
}
