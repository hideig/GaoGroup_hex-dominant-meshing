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

#include<hxt_combine_cpp_api.h>

#include <cassert>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <functional>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <set>
#include <sstream>
#include <stack>
#include <vector>

#include <tet_mesh.h>
#include <basic_types.h>
#include <hxt_combine_cell.h>
#include <candidate_cell.h>
#include <algorithms.h>
#include <exact_cell_evaluation.h>
#include <compute_candidate_cells.h>

#include "hxt_omp.h"

/**
* \author Jeanne Pellerin
*/

namespace HXTCombine {

class StoreCandidateCells{
  /**
   * Amount of hexahedra stored in a buffer before trying.
   *
   * Increasing this value can decrease resource contention, but increases
   * memory consumption.
   */
  static constexpr std::size_t BUFFER_SIZE = 4096;

public:
  StoreCandidateCells(const TetMeshForCombining& tetMesh, double minQuality, vector<HXTCombineCell>& cells) :
    tets_(tetMesh), qualityThreshold_(minQuality),
    filler_(omp_get_max_threads(), cells)
  {
    // Reserving space definitely improves the performances on big models
    // But no real good guess. This depends a lot on the input mesh and on the min quality

    // It depends on the cell types, on the tet mesh, on the min quality !
    filler_.reserve(20*tetMesh.nbVertices());
  }
 // For these hexes..
  void operator()(const vector<VertexIndex>& vertices)
  {
    int id = omp_get_thread_num();

    double quality = cellQuality(tets_, vertices);
    if (quality > qualityThreshold_ && isCellReallyValid(tets_, vertices)){
      CandidateCell &cell = filler_.emplaceBack( id, vertices, tets_);

      // tf_ 输入网格的颜色属性有影响吗？随机化初始多种颜色/不考虑颜色
      //if(!cell.checkTetsTopology() ) {
      if(!cell.checkTetsTopology() || !cell.haveTetsSameColor() || !cell.checkQuadFacetColor()) {
        filler_.popBack(id);
      }
    }
  }
  // For these pyramids..
  void operator()(const vector<VertexIndex>& vertices, TetIndex t0)
  {
    int id = omp_get_thread_num();

    double quality = cellQuality(tets_, vertices);
    if (quality > qualityThreshold_ && isCellReallyValid(tets_, vertices)) {
      CandidateCell &cell = filler_.emplaceBack(id, vertices, t0, tets_);
      if (!cell.checkTetsTopology() || !cell.haveTetsSameColor() || !cell.checkQuadFacetColor()) {
        filler_.popBack(id);
      }
    }
  }

  /**
   * Flushes the buffer of all threads.
   *
   * Must always be called once all cells have been found.
   */
  void flush() {
    filler_.flush();
  }

private:
  const TetMeshForCombining& tets_;
  const double qualityThreshold_;
  ParallelArrayFiller<HXTCombineCell, CandidateCell, CreateCombineCell> filler_;
};


struct CandidateCellStatistics {
  enum CellType { CUBE, BOTELLA, YAMAKAWA, ALL, FALSE_VALID, FAILED, COLOR_TETS, COLOR_FACETS };

  const std::string typeNames_[8]{
    "Cube", "Botella", "Yamakawa", "All Hexes", "False valid",
    "Invalid nb tets", "Invalid tet color", "Invalid facet color" };

  unsigned int data_[8] = { 0,0,0,0,0,0,0,0 };

  void addCell(const TetMeshForCombining& tets, const vector<VertexIndex>& vertices)
  {
    // Min scaled jacobian is above the threshold but the hex is not really valid.
    if (!isCellReallyValid(tets, vertices)) data_[FALSE_VALID]++;
    else {
      CandidateCell hex(vertices, tets);
      if (!hex.checkTetsTopology())        data_[FAILED]++;
      else if (!hex.haveTetsSameColor())   data_[COLOR_TETS]++;
      else if (!hex.checkQuadFacetColor()) data_[COLOR_FACETS]++;
      else {
        if (hex.isCubeDecomposition())      data_[CUBE]++;
        if (hex.isBotellaDecomposition())   data_[BOTELLA]++;
        if (hex.isYamakawaDecomposition())  data_[YAMAKAWA]++;
        data_[ALL]++;
      }
    }
  }
  void addCell(const TetMeshForCombining& tets, const vector<VertexIndex>& vertices, TetIndex t0)
  {
    // Min scaled jacobian is above the threshold but the hex is not really valid.
    if (!isCellReallyValid(tets, vertices)) data_[FALSE_VALID]++;
    else {
      CandidateCell hex(vertices, t0, tets);
      if (!hex.checkTetsTopology())        data_[FAILED]++;
      else if (!hex.haveTetsSameColor())   data_[COLOR_TETS]++;
      else if (!hex.checkQuadFacetColor()) data_[COLOR_FACETS]++;
      else {
        if (hex.isCubeDecomposition())      data_[CUBE]++;
        if (hex.isBotellaDecomposition())   data_[BOTELLA]++;
        if (hex.isYamakawaDecomposition())  data_[YAMAKAWA]++;
        data_[ALL]++;
      }
    }
  }

  void printCategories(std::ofstream& out) const
  {
    for (unsigned int i = 0; i <8; ++i)
      out << std::setw(20) << std::left << typeNames_[i];
  }
  void printData(std::ofstream& out) const
  {
    for (unsigned int i = 0; i <8; ++i)
      out << std::setw(20) << std::left << data_[i];
  }
  void operator+=(const CandidateCellStatistics& rhs)
  {
    for (unsigned int i = 0; i <8; ++i)
      data_[i] += rhs.data_[i];
  }
  void reset() {
    for (unsigned int i = 0; i < 8; ++i)
      data_[i]=0;
  }
};


/**
 * Specifically written to produce results for IMR 2017 and compare our algorithm
 * to the pattern based method
 */
class ClassifyCandidateHexes {
public:
  ClassifyCandidateHexes(const TetMeshForCombining& tetMesh, double minQuality,
    vector<CandidateCellStatistics>& nbHexesPerThreadPerType)
    :tets_(tetMesh), qualityThreshold_(minQuality), stats_(nbHexesPerThreadPerType)
  {
    stats_.resize(omp_get_max_threads());
  }
  void operator()(const vector<VertexIndex>& hexVertices)
  {
    if (cellQuality(tets_, hexVertices) > qualityThreshold_) {
      unsigned int iThread = omp_get_thread_num();
      stats_[iThread].addCell(tets_, hexVertices);
    }
  }
  void operator()(const vector<VertexIndex>& vertices, TetIndex t0)
  {
    double quality = cellQuality(tets_, vertices);
    if (quality > qualityThreshold_ ) {
      unsigned int iThread = omp_get_thread_num();
      stats_[iThread].addCell(tets_, vertices, t0);
    }
  }

private:
  const TetMeshForCombining& tets_;
  const double qualityThreshold_;
  vector<CandidateCellStatistics>& stats_;
};


//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

void countAndCompareCellsIMR(std::ofstream& out, TetMeshForCombining& tetMesh, const vector<double>& minQualities)
{
  CandidateCellStatistics stats;

  out << std::setw(20) << std::left << "Q min"
      << std::setw(20) << std::left << "timing";
  stats.printCategories(out);
  out << std::endl;

  unsigned int nbThreads = omp_get_max_threads();

  for (unsigned int i = 0; i < minQualities.size(); ++i) {
    double minQuality = minQualities[i];

    auto t0 = std::chrono::high_resolution_clock::now();

    vector<CandidateCellStatistics> threadStats(nbThreads);
    ClassifyCandidateHexes counter(tetMesh, minQuality, threadStats);

    computeAllHex(counter, tetMesh, minQuality);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> countTiming(t1 - t0);

    stats.reset();
    for (unsigned int j = 0; j < nbThreads; ++j) {
      stats += threadStats[j];
    }

    out << std::setw(20) << std::left << minQuality
        << std::setw(20) << std::left << countTiming.count();
    stats.printData(out);
    out << std::endl;
  }
  // PRISMS
  for (unsigned int i = 0; i < minQualities.size(); ++i) {
    double minQuality = minQualities[i];

    auto t0 = std::chrono::high_resolution_clock::now();

    vector<CandidateCellStatistics> threadStats(nbThreads);
    ClassifyCandidateHexes counter(tetMesh, minQuality, threadStats);

    computeAllPrisms(counter, tetMesh, minQuality);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> countTiming(t1 - t0);

    stats.reset();
    for (unsigned int j = 0; j < nbThreads; ++j) {
      stats += threadStats[j];
    }
    out << std::setw(20) << std::left << minQuality
      << std::setw(20) << std::left << countTiming.count();
    stats.printData(out);
    out << std::endl;
  }
  //PYRAMIDS
  for (unsigned int i = 0; i < minQualities.size(); ++i) {
    double minQuality = minQualities[i];

    auto t0 = std::chrono::high_resolution_clock::now();

    vector<CandidateCellStatistics> threadStats(nbThreads);
    ClassifyCandidateHexes counter(tetMesh, minQuality, threadStats);

    computeAllPyramids(counter, tetMesh, 0);

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> countTiming(t1 - t0);

    stats.reset();
    for (unsigned int j = 0; j < nbThreads; ++j) {
      stats += threadStats[j];
    }

    out << std::setw(20) << std::left << minQuality
      << std::setw(20) << std::left << countTiming.count();
    stats.printData(out);
    out << std::endl;
  }

  out << std::endl;
}


//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
class hexVol{
public:
  hexVol(std::vector<vec3>& pts){
    pCenter=pts[0];
    LeftBackTop=pts[1];
    LeftFrontTop=pts[2];
    RightBackTop=pts[3];
    RightFrontTop=pts[4];
    LeftFrontBottom=pts[5];
    RightFrontBottom=pts[6];
    RightBackBottom=pts[7];
    LeftBackBottom=pts[8];
  }
  double distance(vec3& p1,vec3& p2)  {
      return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2) + pow(p1[2] - p2[2], 2));
  }

    double cos(double a,double b,double c)  {
      return (a * a + b * b - c * c) / (2 * a * b);
  }
    double Volume() {
      double fVolume = 0.0;
      //劈分成12个不规则四面体；

      //上
      double pPA = distance(LeftBackTop, pCenter);
      double pPB = distance(LeftFrontTop, pCenter);
      double pPC = distance(RightBackTop, pCenter);
      double pAB = distance(LeftBackTop, LeftFrontTop);
      double pAC = distance(LeftBackTop, RightBackTop);
      double pBC = distance(LeftFrontTop, RightBackTop);
      double x = cos(pPA, pPB, pAB);
      double y = cos(pPB, pPC, pBC);
      double z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z))/6;

      pPA = distance(LeftFrontTop, pCenter);
      pPB = distance(RightFrontTop, pCenter);
      pPC = distance(RightBackTop, pCenter);
      pAB = distance(RightFrontTop, LeftFrontTop);
      pAC = distance(LeftFrontTop, RightBackTop);
      pBC = distance(RightFrontTop, RightBackTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      //下
      pPA = distance(LeftFrontBottom, pCenter);
      pPB = distance(RightFrontBottom, pCenter);
      pPC = distance(RightBackBottom, pCenter);
      pAB = distance(RightFrontBottom, LeftFrontBottom);
      pAC = distance(LeftFrontBottom, RightBackBottom);
      pBC = distance(RightFrontBottom, RightBackBottom);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      pPA = distance(LeftFrontBottom, pCenter);
      pPB = distance(RightFrontBottom, pCenter);
      pPC = distance(RightBackBottom, pCenter);
      pAB = distance(RightFrontBottom, LeftFrontBottom);
      pAC = distance(LeftFrontBottom, RightBackBottom);
      pBC = distance(RightFrontBottom, RightBackBottom);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;
      //左
      pPA = distance(LeftFrontTop, pCenter);
      pPB = distance(LeftBackTop, pCenter);
      pPC = distance(LeftBackBottom, pCenter);
      pAB = distance(LeftFrontTop, LeftBackTop);
      pAC = distance(LeftFrontTop, LeftBackBottom);
      pBC = distance(LeftBackTop, LeftBackBottom);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      pPA = distance(LeftFrontBottom, pCenter);
      pPB = distance(LeftBackBottom, pCenter);
      pPC = distance(LeftFrontTop, pCenter);
      pAB = distance(LeftFrontBottom, LeftBackBottom);
      pAC = distance(LeftFrontBottom, LeftFrontTop);
      pBC = distance(LeftBackBottom, LeftFrontTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;
      //右
      pPA = distance(RightFrontTop, pCenter);
      pPB = distance(RightBackTop, pCenter);
      pPC = distance(RightBackBottom, pCenter);
      pAB = distance(RightFrontTop, RightBackTop);
      pAC = distance(RightFrontTop, RightBackBottom);
      pBC = distance(RightBackTop, RightBackBottom);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      pPA = distance(RightFrontBottom, pCenter);
      pPB = distance(RightBackBottom, pCenter);
      pPC = distance(RightFrontTop, pCenter);
      pAB = distance(RightFrontBottom, RightBackBottom);
      pAC = distance(RightFrontBottom, RightFrontTop);
      pBC = distance(RightBackBottom, RightFrontTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;
      //前
      pPA = distance(RightFrontTop, pCenter);
      pPB = distance(RightFrontBottom, pCenter);
      pPC = distance(LeftFrontTop, pCenter);
      pAB = distance(RightFrontBottom, RightFrontTop);
      pAC = distance(RightFrontTop, LeftFrontTop);
      pBC = distance(RightFrontBottom, LeftFrontTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      pPA = distance(LeftFrontBottom, pCenter);
      pPB = distance(RightFrontBottom, pCenter);
      pPC = distance(LeftFrontTop, pCenter);
      pAB = distance(RightFrontBottom, LeftFrontBottom);
      pAC = distance(LeftFrontBottom, LeftFrontTop);
      pBC = distance(RightFrontBottom, LeftFrontTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;
      //后
      pPA = distance(RightBackTop, pCenter);
      pPB = distance(RightBackBottom, pCenter);
      pPC = distance(LeftBackTop, pCenter);
      pAB = distance(RightBackBottom, RightBackTop);
      pAC = distance(RightBackTop, LeftBackTop);
      pBC = distance(RightBackBottom, LeftBackTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;

      pPA = distance(LeftBackBottom, pCenter);
      pPB = distance(RightBackBottom, pCenter);
      pPC = distance(LeftBackTop, pCenter);
      pAB = distance(RightBackBottom, LeftBackBottom);
      pAC = distance(LeftBackBottom, LeftBackTop);
      pBC = distance(RightBackBottom, LeftBackTop);
      x = cos(pPA, pPB, pAB);
      y = cos(pPB, pPC, pBC);
      z = cos(pPA, pPC, pAC);
      fVolume += (pPA * pPB * pPC * sqrt(1 + 2 * x * y * z - x * x - y * y - z * z)) / 6;
     // std::cout << "fVolume: " << fVolume << std:: endl;
      return fVolume;
  }

private:
    vec3 pCenter;
    vec3 LeftBackTop;
    vec3 LeftFrontTop;
    vec3 RightBackTop;
    vec3 RightFrontTop;
    vec3 LeftFrontBottom;
    vec3 RightFrontBottom;
    vec3 RightBackBottom;
    vec3 LeftBackBottom;
};

 double HXTCombineCellStore::computeHexesVolume(){
   double res = 0.0;
   int size = (this->hexes()).size() ;
   int ct = 0;
   for(int i = 0; i < size; i++){
      if(this->selectedHexes()[i]){
        std::vector<vec3> pts;
        ct++;
        HXTCombineCell& cell = this->hexes()[i];
        vec3 center;
        for(int j = 0; j < 8; j++){
          center+=this->mesh_.point(cell.vertexes[j]);
        }
        center/=8;
        pts.push_back(center);
        for(int j = 0; j < 8; j++){
          pts.push_back(this->mesh_.point(cell.vertexes[j]));
        }
        hexVol* temp = new hexVol(pts);
        double ddd = temp->Volume();
        if(ddd >0 && ddd < 200)res+=ddd;
        delete temp;
    }
   }
   std::cout <<  "res: " << res << std::endl; 
   std::cout << "hex number: " << ct << std::endl; 
 return res;
 }

void HXTCombineCellStore::computeHexes(double minQuality) {
  StoreCandidateCells store(mesh_, minQuality, hexes());

  computeAllHex(store, mesh_, minQuality);

  store.flush();
}
void HXTCombineCellStore::computePrisms(double minQuality) {
  StoreCandidateCells store(mesh_, minQuality, prisms());

  computeAllPrisms(store, mesh_, minQuality);

  store.flush();
}

void HXTCombineCellStore::computePyramids(double minQuality) {
  StoreCandidateCells store(mesh_, minQuality, pyramids());

  computeAllPyramids(store, mesh_, 0);

  store.flush();
}


void computePotentialHexes(
  TetMeshForCombining& tetMesh, double hexQualityThreshold, vector<HXTCombineCell>& hexes)
{
  StoreCandidateCells store(tetMesh, hexQualityThreshold, hexes);
  computeAllHex(store, tetMesh, hexQualityThreshold);
  store.flush();
}

void computePotentialPrisms(
  TetMeshForCombining& tetMesh, double prismQualityThreshold, vector<HXTCombineCell>& prisms)
{
  StoreCandidateCells store(tetMesh, prismQualityThreshold, prisms);
  computeAllPrisms(store, tetMesh, prismQualityThreshold);
  store.flush();
}



//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

double cellQualityAPI(const TetMeshWrapper& tets, const HXTCombineCell& cell) {
  if (cell.isHex()) {
    vec3 points[8];
    for (HexVertexIndex j = 0; j < 8; ++j) {
      points[j] = tets.point(cell.vertex(j));
    }
    return hexFullQuality(points);
  }
  else if (cell.isPrism()) {
    vector<vec3> points(6);
    for (PrismVertexIndex j = 0; j <6; ++j) {
      points[j] = tets.point(cell.vertex(j));
    }
    return cellApproximateQuality<Prism>(points.data());
  }
  else if (cell.isPyramid()) {
    vector<vec3> points(5);
    for (CellVertexIndex j = 0; j < 5; ++j) {
      points[j] = tets.point(cell.vertex(j));
    }
    return cellApproximateQuality<Pyramid>(points.data());
  }
  else return 0;
}

void computeCellQualityVector(
  const TetMeshWrapper& tets,
  const std::vector<HXTCombineCell>& cells,
  double* qualities)
{
#pragma omp parallel for
  for (int i = 0; i < (int)cells.size(); ++i) {
    qualities[i] = cellQualityAPI(tets, cells[i]);
  }
}

}
