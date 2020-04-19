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



//#include "basic_types.h"

#include <cassert>
#include <map>
#include <limits.h>
#include "cell_types.h"
/**
* \author Jeanne Pellerin
*/

namespace HXTCombine {
  // These declarations are necessary to compile with GCC
  constexpr CellVertexIndex Teta::facetVertex[][3];
  constexpr CellVertexIndex Teta::edgeVertex[][2];

  constexpr CellVertexIndex Hexa::facetVertex[][4];
  constexpr CellVertexIndex Hexa::edgeVertex[][2];
  constexpr CellVertexIndex Hexa::quadFacetTriangleVertex[][3];
  constexpr CellVertexIndex Hexa::quadFacetDiagonalVertex[][2];
  constexpr CellFacetIndex  Hexa::vertexAdjacentFacet[][3];
  constexpr CellVertexIndex Hexa::vertexAdjacentVertex[][3];
  constexpr CellVertexIndex Hexa::orientedFacetVertex[][4];
  constexpr unsigned int Hexa::vertexToBoundaryTriangle[][8][8];


  constexpr CellVertexIndex Pyramid::facetVertex[][4];
  constexpr CellVertexIndex Pyramid::edgeVertex[][2];
  constexpr CellVertexIndex Pyramid::quadFacetTriangleVertex[][3];
  constexpr CellVertexIndex Pyramid::quadFacetDiagonalVertex[][2];
  constexpr CellVertexIndex Pyramid::vertexAdjacentVertex[][4];
  constexpr CellVertexIndex Pyramid::vertexAdjacentFacet[][4];
  constexpr unsigned int Pyramid::vertexToBoundaryTriangle[][5][5];


  constexpr CellVertexIndex Prism::facetVertex[][4];
  constexpr CellVertexIndex Prism::edgeVertex[][2];
  constexpr CellVertexIndex Prism::quadFacetTriangleVertex[][3];
  constexpr CellVertexIndex Prism::quadFacetDiagonalVertex[][2];
  constexpr CellFacetIndex  Prism::vertexAdjacentFacet[][3];
  constexpr CellVertexIndex Prism::vertexAdjacentVertex[][3];
  constexpr CellVertexIndex Prism::orientedFacetVertex[][4];
  constexpr unsigned int Prism::vertexToBoundaryTriangle[][6][6];
}


namespace {

template<class T> 
CellVertexIndex minVertexIndex(const VertexIndex* vertices) 
{
  CellVertexIndex result = 0;
  VertexIndex min = INT_MAX;
  for (CellVertexIndex i = 0; i < T::nbVertices; ++i) {
    if (vertices[i] < min) {
      min = vertices[i];
      result = i;
    }
  }
  return result;
}

template<class T>
CellFacetVertexIndex vertexIndexInFacet(const VertexIndex* vertices, CellFacetIndex f, unsigned int v)
{
  unsigned int nbFacetVertices = f < T::nbTriFacets ? 3 : 4;
  for (CellFacetVertexIndex lv = 0; lv < nbFacetVertices; ++lv) {
    CellVertexIndex cur = T::facetVertex[f][lv];
    if (vertices[v] == vertices[cur]) return lv;
  }
  return NO_ID;
}


template<class T>
CellVertexIndex incidentVertex(trindex vertices) {
  for (CellVertexIndex v = 0; v < T::nbVertices; ++v) {
    trindex curAdjacent(
      T::vertexAdjacentVertex[v][0],
      T::vertexAdjacentVertex[v][1],
      T::vertexAdjacentVertex[v][2]);
    if (curAdjacent == vertices) return v;
  }
  assert(false);
  return NO_ID;
}

template<class T>
CellVertexIndex lastVertexInQuadFacet(trindex vertices)
{
  for (CellFacetIndex f = T::nbTriFacets; f < T::nbFacets; ++f) {
    trindex t0(T::facetVertex[f][0], T::facetVertex[f][1], T::facetVertex[f][3]);
    if (t0 == vertices) return T::facetVertex[f][2];

    trindex t1(T::facetVertex[f][0], T::facetVertex[f][1], T::facetVertex[f][2]);
    if (t1 == vertices) return T::facetVertex[f][3];

    trindex t2(T::facetVertex[f][0], T::facetVertex[f][2], T::facetVertex[f][3]);
    if (t2 == vertices) return T::facetVertex[f][1];

    trindex t3(T::facetVertex[f][1], T::facetVertex[f][2], T::facetVertex[f][3]);
    if (t3 == vertices) return T::facetVertex[f][0];
  }
  assert(false);
  return NO_ID;
}

template<class T>
std::map<VertexIndex, CellVertexIndex> diagonalOppositeVerticesIncidentQuadFacets(
  const VertexIndex* vertices, CellVertexIndex v0)
{
  std::map<VertexIndex, CellVertexIndex> result;

  for (unsigned int i = 0; i < 3; i++) {
    CellFacetIndex f = T::vertexAdjacentFacet[v0][i];
    if (f < T::nbTriFacets) continue;

    CellFacetVertexIndex v0InFacet = vertexIndexInFacet<T>(vertices, f, v0);
    assert(v0InFacet != NO_ID);

    CellFacetVertexIndex oppositeToV0 = (v0InFacet+2)%4;

    CellVertexIndex v = T::facetVertex[f][oppositeToV0];
    result.insert(std::make_pair(vertices[v], v));
  }
  assert( (T::nbVertices == 6 && result.size()==2) || (T::nbVertices == 8 && result.size()==3));
  return result;
}

} // anonymous namespace

namespace HXTCombine {

void computePrismVertexNormalizationPermutation(
  const VertexIndex vertices[6], vector<PrismVertexIndex>& permutation)
{
  permutation.resize(6, NO_ID);
  permutation[0] = minVertexIndex<Prism>(vertices);

  std::map<VertexIndex, PrismVertexIndex> v1v2 =
    diagonalOppositeVerticesIncidentQuadFacets<Prism>(vertices, permutation[0]);
  unsigned int i = 1;
  for (auto it = v1v2.begin(); it != v1v2.end(); ++it) {
    permutation[i] = it->second;
    ++i;
  }
  trindex v3Adj(permutation[0], permutation[1], permutation[2]);   // v3 is adjacent to v0, v1, v2
  permutation[3] = incidentVertex<Prism>(v3Adj);

  trindex v4Quad(permutation[0], permutation[1], permutation[3]);  // v4 is in quad facet v0, v3, v1
  trindex v5Quad(permutation[0], permutation[2], permutation[3]);  // v5 is in quad facet v0, v3, v2
  permutation[4] = lastVertexInQuadFacet<Prism>(v4Quad);
  permutation[5] = lastVertexInQuadFacet<Prism>(v5Quad);
}


void computeHexVertexNormalizationPermutation(
  const VertexIndex vertices[8], vector<HexVertexIndex>& permutation)
{
  permutation.resize(8, NO_ID);
  permutation[0] = minVertexIndex<Hexa>(vertices);

  std::map<VertexIndex, HexVertexIndex> v1v2v3 =
    diagonalOppositeVerticesIncidentQuadFacets<Hexa>(vertices, permutation[0]);
  unsigned int i = 1;
  for (auto it = v1v2v3.begin(); it != v1v2v3.end(); ++it) {
    permutation[i] = it->second;
    ++i;
  }
  trindex v4Adj(permutation[0], permutation[1], permutation[2]);  // v4 is adjacent to v0, v1, v2
  trindex v5Adj(permutation[0], permutation[1], permutation[3]);  // v5 is adjacent to v0, v1, v3
  trindex v6Adj(permutation[0], permutation[2], permutation[3]);  // v6 is adjacent to v0, v2, v3
  trindex v7Adj(permutation[1], permutation[2], permutation[3]);  // v7 is adjacent to v1, v2, v3
  permutation[4] = incidentVertex<Hexa>(v4Adj);
  permutation[5] = incidentVertex<Hexa>(v5Adj);
  permutation[6] = incidentVertex<Hexa>(v6Adj);
  permutation[7] = incidentVertex<Hexa>(v7Adj);
}

}

