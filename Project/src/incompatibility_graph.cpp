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

#include "hxt_combine_cpp_api.h"

#include <cassert>
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <thread>
#include <vector>

#include "algorithms.h"
#include "candidate_cell.h"
//#include <hxt_mwis.h>

#include "cell_types.h"
#include "incompatibility_graph.h"

#include "hxt_omp.h"


/**
* \author Kilian Verhetsel
*/

namespace {

using namespace HXTCombine;
using std::array;


struct VertexToCell {
  VertexToCell() {}
  VertexToCell(VertexIndex v, CellIndex c) :
    vertex(v), cell(c) {};

  bool operator<(const VertexToCell& rhs) const
  {
    return vertex < rhs.vertex;
  }
  VertexIndex vertex;
  CellIndex cell;
};

vector<VertexToCell>
fillVertexToCell(const TetMeshWrapper& tets,
                 const std::vector<HXTCombineCell>& cells)
{
  vector<VertexToCell> vertexToCell;

  const int threadNum = omp_get_max_threads();
  ParallelArrayFiller<VertexToCell, VertexToCell, Identity<VertexToCell>, 4096*32> filler(
    threadNum, vertexToCell);
  filler.reserve(cells.size() * 8);

  #pragma omp parallel for
  for (int i = 0; i < (int)cells.size(); ++i) {
    const HXTCombineCell& curCell = cells[i];
    const VertexIndex *curVertices = curCell.vertices();
    for (CellVertexIndex v = 0; v < curCell.nbVertices(); ++v) {
      filler.emplaceBack(omp_get_thread_num(), curVertices[v], i);
    }
  }

  filler.flush();
  return vertexToCell;
}

//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------

class IncompatibilityDetector {
  IncompatibilityAdder &adder_;

  const std::vector<VertexToCell> &vertexToCell_;
  const std::vector<HXTCombineCell> &cells_;
  std::size_t threadId_, threadNum_;

public:
  IncompatibilityDetector(IncompatibilityAdder &adder,
                          const std::vector<VertexToCell> &vertexToCell,
                          const std::vector<HXTCombineCell> &cells,
                          std::size_t threadId, std::size_t threadNum):
    adder_(adder),
    vertexToCell_(vertexToCell),
    cells_(cells),
    threadId_(threadId), threadNum_(threadNum)
    {}

  void operator()() {
    std::size_t produced = 0;

    static constexpr std::size_t BLOCK_SIZE = 8192;

    std::size_t step = threadNum_ * BLOCK_SIZE;

    const std::size_t n = vertexToCell_.size();
    for (std::size_t currentStart = threadId_ * BLOCK_SIZE;
         currentStart < n; currentStart += step) {
      produced += exploreRange(currentStart,
                               std::min(currentStart + BLOCK_SIZE, n));
    }
  }

private:
  std::size_t exploreRange(std::size_t start, std::size_t end) {
    std::size_t produced = 0;

    std::size_t from, to;
    for (from = moveToStart(start, end); from < end; from = to) {
      to = findVertexEnd(from);
      VertexIndex v = vertexToCell_[from].vertex;

      for (unsigned int i = from; i < to; ++i) {
        CellIndex h0 = vertexToCell_[i].cell;
        const HXTCombineCell& cell0 = cells_[h0];
        for (unsigned int j = i+1; j < to; ++j) {
          CellIndex h1 = vertexToCell_[j].cell;

          if (h0 == h1) continue;
          const HXTCombineCell& cell1 = cells_[h1];

          if ((h0 < h1 && firstSharedVertex(cell0, cell1) != v) ||
              (h1 < h0 && firstSharedVertex(cell1, cell0) != v)) {
            continue;
          }

          if (!areCellsCompatible(cell0, cell1)) {
            produced += 2;
            adder_.addPair(h0, h1);
            adder_.addPair(h1, h0);
          }
        }
      }
    }

    return produced;
  }

  std::size_t moveToStart(std::size_t i, std::size_t limit) {
    while (i != 0 &&
           vertexToCell_[i-1].cell == vertexToCell_[i].cell &&
           i < limit) {
      i++;
    }
    return i;
  }

  std::size_t findVertexEnd(std::size_t from) {
    VertexIndex v = vertexToCell_[from].vertex;
    while (from < vertexToCell_.size() && vertexToCell_[from].vertex == v)
      from++;
    return from;
  }
};

}


namespace HXTCombine {

//void incompatibilityGraph(const TetMeshWrapper& tets, const std::vector<HXTCombineCell>& cells, HXTGraph& graph)
//{
//  if (cells.empty()) return;
//
//  std::vector<VertexToCell> vertexToCell = fillVertexToCell(tets, cells);
//  std::size_t numCells = cells.size();
//  std::sort(vertexToCell.begin(), vertexToCell.end());
//
//  HXTIndex **adjacencyList;
//  if (hxtMalloc((void**)&adjacencyList, numCells*sizeof(*adjacencyList)) < 0) {
//    throw std::bad_alloc();
//  }
//
//  HXTIndex *degrees;
//  if (hxtMalloc((void**)&degrees, numCells*sizeof(*degrees)) < 0) {
//    hxtFree((void*)&adjacencyList);
//    throw std::bad_alloc();
//  }
//
//  HXTIndex *capacities;
//  if (hxtMalloc((void**)&capacities, numCells*sizeof(*capacities)) < 0) {
//    hxtFree((void*)&adjacencyList);
//    hxtFree((void*)&degrees);
//    throw std::bad_alloc();
//  }
//
//  for (HXTIndex i = 0; i < numCells; ++i) {
//    degrees[i] = 0;
//    capacities[i] = 16;
//    if (hxtMalloc((void*)&adjacencyList[i],
//                  capacities[i]*sizeof(*adjacencyList[i])) < 0) {
//      for (HXTIndex j = 0; j < i; j++)
//        hxtFree((void*)&adjacencyList[i]);
//      hxtFree((void*)&adjacencyList);
//      hxtFree((void*)&degrees);
//      throw std::bad_alloc();
//    }
//  }
//
//  addToGraph(adjacencyList, degrees, capacities,
//             [&](IncompatibilityAdder &adder, int i, int n) {
//               IncompatibilityDetector(adder, vertexToCell, cells, i, n)();
//             });
//
//  // Sort what we gathered and remove duplicatas.
//  // Is this efficient?
//  #pragma omp parallel for
//  for (int i = 0; i < (int)numCells; ++i) {
//    std::sort(adjacencyList[i], adjacencyList[i] + degrees[i]);
//    degrees[i] = std::unique(adjacencyList[i], adjacencyList[i] + degrees[i]) -
//      adjacencyList[i];
//  }
//
//  hxtFree((void*)&capacities);
//
//  if (hxtGraphCreate(&graph, numCells, adjacencyList, degrees) < 0)
//    throw std::bad_alloc();
//
//  for (HXTIndex i = 0; i < numCells; ++i) {
//    hxtFree((void*)&adjacencyList[i]);
//  }
//  hxtFree((void*)&adjacencyList);
//  hxtFree((void*)&degrees);
//  hxtFree((void*)&capacities);
//}

} // namespace
