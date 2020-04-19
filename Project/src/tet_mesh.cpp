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

#include "tet_mesh.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <utility>
#include <string>
#include <string.h>

#include "basic_types.h"

 /**
 * \author Jeanne Pellerin
 */

using namespace std;
namespace HXTCombine {

	trindex tetFacet(const MeshStore& tets, TetIndex t, TetFacetIndex lf) {
		VertexIndex v0 = tets.tetCorners[4 * t + Teta::facetVertex[lf][0]];
		VertexIndex v1 = tets.tetCorners[4 * t + Teta::facetVertex[lf][1]];
		VertexIndex v2 = tets.tetCorners[4 * t + Teta::facetVertex[lf][2]];
		return trindex(v0, v1, v2);
	}

	/**
	 * 计算体与面的邻接关系，耗时耗内存
	* This is expensive and requires quite a lot of memory
	*
	* Is it the fastest way to compute these adjacencies. How can we parallelize that ?
	*/
	void computeTetAdjacencies(MeshStore& store) {

		AdjacentTetIndex NO_ADJACENT = store.NO_ADJACENT_ID;

		unsigned int nbTets = store.tetCorners.size() / 4;
		unsigned int nbVertices = store.vertices.size();

		// 1. Chain corners around vertices
		vector<TetCornerIndex> nextTetCornerAroundVertex(store.tetCorners.size(), NO_ID);
		vector<TetCornerIndex> vertexToTetCorner(nbVertices, NO_ID);
		for (TetIndex t = 0; t < nbTets; ++t) {
			for (TetVertexIndex lv = 0; lv < 4; ++lv) {
				TetCornerIndex c = 4 * t + lv;
				VertexIndex v = store.tetCorners[c];
				nextTetCornerAroundVertex[c] = vertexToTetCorner[v];
				vertexToTetCorner[v] = c;
			}
		}

		// 2. Get tet adjacencies
		store.tetAdjacent.resize(4 * nbTets, NO_ADJACENT);
		for (TetIndex t = 0; t < nbTets; ++t) {
			for (TetFacetIndex f = 0; f < 4; ++f) {
				if (store.tetAdjacent[4 * t + f] == NO_ADJACENT) {
					trindex facet = tetFacet(store, t, f);
					VertexIndex v0 = facet.indices[0];

					for (TetCornerIndex c = vertexToTetCorner[v0]; c != NO_ID; c = nextTetCornerAroundVertex[c]) {
						TetIndex tetAroundV0 = c / 4;
						if (tetAroundV0 != t) {
							// Find the adjacent facet index in tetAroundV0
							TetFacetIndex adjacentF = NO_ADJACENT;
							if (tetFacet(store, tetAroundV0, 0) == facet) adjacentF = 0;
							else if (tetFacet(store, tetAroundV0, 1) == facet) adjacentF = 1;
							else if (tetFacet(store, tetAroundV0, 2) == facet) adjacentF = 2;
							else if (tetFacet(store, tetAroundV0, 3) == facet) adjacentF = 3;

							if (adjacentF != NO_ADJACENT) {
								store.tetAdjacent[4 * t + f] = 4 * tetAroundV0 + adjacentF;
								store.tetAdjacent[4 * tetAroundV0 + adjacentF] = 4 * t + f;
								break;
							}
						}
					}
				}
			}
		}
	}

	void buildhxtVertexFromCoordinates(const std::vector<double>& points, std::vector<double>& hxtvertices)
	{
		size_t nbVertices = points.size() / 3;
		hxtvertices.resize(4 * nbVertices);
		for (unsigned int i = 0; i < nbVertices; ++i) {
			hxtvertices[4 * i] = points[3 * i];
			hxtvertices[4 * i + 1] = points[3 * i + 1];
			hxtvertices[4 * i + 2] = points[3 * i + 2];
		}
	}

	/** Reading .mesh file - only the vertices - triangles - quads
	*/
	void myReadFileMESH(const MatrixXf &V_, const MatrixXu &F_, const MatrixXu &T_, HXTCombine::MeshStore& mesh, bool read_triangles)	{
		std::vector<double>& vertexVector = mesh.vertices;
		std::vector<VertexIndex>& triangles = mesh.triangleCorners;
		std::vector<VertexIndex>& tets = mesh.tetCorners;
		std::vector<VertexIndex>& hex = mesh.hexCorners;

		std::vector<ColorIndex>& triangleColors = mesh.triangleColors;
		std::vector<ColorIndex>& hexColors = mesh.hexColors;
		std::vector<ColorIndex>& tetColors = mesh.tetColors;

		int nbv = V_.cols();
		vertexVector.resize(4 * nbv);
		for (int i = 0; i < nbv; i++) {
			double x = V_(0, i), y = V_(1, i), z = V_(2, i);
			vertexVector[4 * i + 0] = x;
			vertexVector[4 * i + 1] = y;
			vertexVector[4 * i + 2] = z;
		}
		
		int nbf = F_.cols();
		triangles.reserve(3 * nbf);
		triangleColors.reserve(nbf);
		for (int i = 0; i < nbf; i++) {
			int n[3], cl = 1;
			for (int j = 0; j < 3; j++) {
				n[j] = F_(j, i);
				triangles.push_back(n[j]);
			}
			triangleColors.push_back(cl);
		}
		
		int nbt = T_.cols();
		tets.reserve(4 * nbt);
		tetColors.reserve(nbt);
		for (int i = 0; i < nbt; i++) {
			int n[4], cl = 1;
			for (int j = 0; j < 4; j++) {
				n[j] = T_(j, i);
				tets.push_back(n[j]);
			}
			tetColors.push_back(cl);
		}
		// Compute adjacencies
		computeTetAdjacencies(mesh);
	}

	void readFileMESH(const std::string& filename, HXTCombine::MeshStore& mesh, bool read_triangles)
	{
		FILE *fp = std::fopen(filename.c_str(), "r");
		if (!fp) {
			std::cout << "Unable to open file " << filename << std::endl;
			std::cout << strerror(errno) << std::endl;
			return;
		}

		char buffer[256];
		if (!fgets(buffer, sizeof(buffer), fp)) { fclose(fp); return; }

		char str[256];
		int format;
		sscanf(buffer, "%s %d", str, &format);

		std::vector<double>& vertexVector = mesh.vertices;

		std::vector<VertexIndex>& triangles = mesh.triangleCorners;
		std::vector<VertexIndex>& quads = mesh.quadCorners;
		std::vector<VertexIndex>& tets = mesh.tetCorners;
		std::vector<VertexIndex>& hex = mesh.hexCorners;

		std::vector<ColorIndex>& triangleColors = mesh.triangleColors;
		std::vector<ColorIndex>& quadColors = mesh.quadColors;
		std::vector<ColorIndex>& hexColors = mesh.hexColors;
		std::vector<ColorIndex>& tetColors = mesh.tetColors;


		while (!feof(fp)) {
			if (!fgets(buffer, 256, fp)) break;
			if (buffer[0] != '#') { // skip comments and empty lines
				str[0] = '\0';
				sscanf(buffer, "%s", str);
				if (!strncmp(buffer, "Dimension 3", 11)) {
				}else if (!strcmp(str, "Dimension")) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
				}else if (!strcmp(str, "Vertices")) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
					int nbv;
					sscanf(buffer, "%d", &nbv);
					std::cout << nbv << " vertices  ";
					vertexVector.resize(4 * nbv);
					for (int i = 0; i < nbv; i++) {
						if (!fgets(buffer, sizeof(buffer), fp)) break;
						int dum;
						double x, y, z;
						sscanf(buffer, "%lf %lf %lf %d", &x, &y, &z, &dum);
						vertexVector[4 * i + 0] = x;
						vertexVector[4 * i + 1] = y;
						vertexVector[4 * i + 2] = z;
					}
				}else if (!strcmp(str, "Triangles") && read_triangles) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
					int nbe;
					sscanf(buffer, "%d", &nbe);
					std::cout << nbe << " triangles  ";
					triangles.reserve(3 * nbe);
					triangleColors.reserve(nbe);
					for (int i = 0; i < nbe; i++) {
						if (!fgets(buffer, sizeof(buffer), fp)) break;
						int n[3], cl;
						sscanf(buffer, "%d %d %d %d", &n[0], &n[1], &n[2], &cl);
						for (int j = 0; j < 3; j++) {
							n[j]--;
							triangles.push_back(n[j]);
						}
						triangleColors.push_back(cl);
					}
				}else if (!strcmp(str, "Quadrilaterals")) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
					int nbe;
					sscanf(buffer, "%d", &nbe);
					std::cout << nbe << " quads  ";
					quads.reserve(4 * nbe);
					quadColors.reserve(nbe);
					for (int i = 0; i < nbe; i++) {
						if (!fgets(buffer, sizeof(buffer), fp)) break;
						int n[4], cl;
						sscanf(buffer, "%d %d %d %d %d", &n[0], &n[1], &n[2], &n[3], &cl);
						for (int j = 0; j < 4; j++) {
							n[j]--;
							quads.push_back(n[j]);
						}
						quadColors.push_back(cl);
					}
				}else if (!strcmp(str, "Tetrahedra")) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
					int nbe;
					sscanf(buffer, "%d", &nbe);
					std::cout << nbe << " tetrahedra  ";
					tets.reserve(4 * nbe);
					tetColors.reserve(nbe);
					for (int i = 0; i < nbe; i++) {
						if (!fgets(buffer, sizeof(buffer), fp)) break;
						int n[4], cl;
						sscanf(buffer, "%d %d %d %d %d", &n[0], &n[1], &n[2], &n[3], &cl);
						for (int j = 0; j < 4; j++) {
							n[j]--;
							tets.push_back(n[j]);
						}
						tetColors.push_back(cl);
					}
				}else if (!strcmp(str, "Hexahedra")) {
					if (!fgets(buffer, sizeof(buffer), fp)) break;
					int nbe;
					sscanf(buffer, "%d", &nbe);
					std::cout << nbe << " hexahedra ";
					hex.reserve(8 * nbe);
					hexColors.reserve(nbe);
					for (int i = 0; i < nbe; i++) {
						if (!fgets(buffer, sizeof(buffer), fp)) break;
						int n[8], cl;
						sscanf(buffer, "%d %d %d %d %d %d %d %d %d", &n[0], &n[1], &n[2], &n[3],
							&n[4], &n[5], &n[6], &n[7], &cl);
						for (int j = 0; j < 8; j++) {
							n[j]--;
							hex.push_back(n[j]);
						}
						hexColors.push_back(cl);
					}
				}
			}
		}
		std::cout << std::endl << std::endl;

		fclose(fp);

		// Compute adjacencies
		computeTetAdjacencies(mesh);
	}


	//void myreadFileMESH(const std::string& filename, HXTCombine::MeshStore& mesh, bool read_triangles)
	//{
	//	FILE *fp = std::fopen(filename.c_str(), "r");
	//	if (!fp) {
	//		std::cout << "Unable to open file " << filename << std::endl;
	//		return;
	//	}

	//	char buffer[256];
	//	if (!fgets(buffer, sizeof(buffer), fp)) { fclose(fp); return; }

	//	char str[256];
	//	int format;
	//	sscanf(buffer, "%s %d", str, &format);

	//	std::vector<double>& vertexVector = mesh.vertices;

	//	std::vector<VertexIndex>& triangles = mesh.triangleCorners;
	//	std::vector<VertexIndex>& quads = mesh.quadCorners;
	//	std::vector<VertexIndex>& tets = mesh.tetCorners;
	//	std::vector<VertexIndex>& hex = mesh.hexCorners;

	//	std::vector<ColorIndex>& triangleColors = mesh.triangleColors;
	//	std::vector<ColorIndex>& quadColors = mesh.quadColors;
	//	std::vector<ColorIndex>& hexColors = mesh.hexColors;
	//	std::vector<ColorIndex>& tetColors = mesh.tetColors;


	//	while (!feof(fp)) {
	//		if (!fgets(buffer, 256, fp)) break;
	//		if (buffer[0] != '#') { // skip comments and empty lines
	//			str[0] = '\0';
	//			sscanf(buffer, "%s", str);
	//			if (!strcmp(str, "Vertices")) {
	//				if (!fgets(buffer, sizeof(buffer), fp)) break;
	//				int nbv;
	//				sscanf(buffer, "%d", &nbv);
	//				std::cout << nbv << " vertices  ";
	//				vertexVector.resize(4 * nbv);
	//				for (int i = 0; i < nbv; i++) {
	//					if (!fgets(buffer, sizeof(buffer), fp)) break;
	//					double x, y, z;
	//					sscanf(buffer, "%lf %lf %lf ", &x, &y, &z);
	//					vertexVector[4 * i + 0] = x;
	//					vertexVector[4 * i + 1] = y;
	//					vertexVector[4 * i + 2] = z;
	//				}
	//			}
	//			else if (!strcmp(str, "Triangles") && read_triangles) {
	//				if (!fgets(buffer, sizeof(buffer), fp)) break;
	//				int nbe;
	//				sscanf(buffer, "%d", &nbe);
	//				std::cout << nbe << " triangles  ";
	//				triangles.reserve(3 * nbe);
	//				triangleColors.reserve(nbe);
	//				for (int i = 0; i < nbe; i++) {
	//					if (!fgets(buffer, sizeof(buffer), fp)) break;
	//					int n[3], cl;
	//					sscanf(buffer, "%d %d %d %d", &cl, &n[0], &n[1], &n[2]);
	//					for (int j = 0; j < 3; j++) {
	//						n[j]--;
	//						triangles.push_back(n[j]);
	//					}
	//					triangleColors.push_back(cl);
	//				}
	//			}
	//			else if (!strcmp(str, "Tetrahedra")) {
	//				if (!fgets(buffer, sizeof(buffer), fp)) break;
	//				int nbe;
	//				sscanf(buffer, "%d", &nbe);
	//				std::cout << nbe << " tetrahedra  ";
	//				tets.reserve(4 * nbe);
	//				tetColors.reserve(nbe);
	//				for (int i = 0; i < nbe; i++) {
	//					if (!fgets(buffer, sizeof(buffer), fp)) break;
	//					int n[4], cl;
	//					sscanf(buffer, "%d %d %d %d %d", &cl, &n[0], &n[1], &n[2], &n[3]);
	//					for (int j = 0; j < 4; j++) {
	//						n[j]--;
	//						tets.push_back(n[j]);
	//					}
	//					tetColors.push_back(cl);
	//				}
	//			}
	//		}
	//	}
	//	std::cout << std::endl << std::endl;

	//	fclose(fp);

	//	// Compute adjacencies
	//	computeTetAdjacencies(mesh);
	//}
	void TetMeshForCombining::initConnectivityTables()	{
		fillFacets();
		chainTetCornersAroundVertices();
		fillVertexToVertices();
		if (nbTriangles() > 0) connectTetsToTriangles();
	}

	/*
	* This is not bad and works well, but how can I parallelize it ??
	*
	* \todo See if I can use Celestin's version for this
	*/
	void TetMeshForCombining::chainTetCornersAroundVertices() {
		vertexToTetCorner_.resize(nbVertices(), NO_ID);
		nextTetCornerAroundVertex_.resize(4 * nbTets(), NO_ID);

		for (TetIndex t = 0; t < nbTets(); ++t) {
			for (TetVertexIndex lv = 0; lv < 4; ++lv) {
				TetCornerIndex c = 4 * t + lv;
				VertexIndex v = vertex(t, lv);
				nextTetCornerAroundVertex_[c] = vertexToTetCorner_[v];
				vertexToTetCorner_[v] = c;
			}
		}
	}


	void TetMeshForCombining::fillVertexToVertices(){
		vertexToVertices_.resize(nbVertices());
		vector<VertexIndex> temp;
		for (VertexIndex v = 0; v < nbVertices(); ++v) {
			temp.resize(0);
			TetCornerIndex c = vertexToTetCorner_[v];
			while (c != NO_ID) {
				TetIndex t = c / 4;
				TetVertexIndex lv = c % 4;
				temp.push_back(vertex(t, Teta::facetVertex[lv][0]));
				temp.push_back(vertex(t, Teta::facetVertex[lv][1]));
				temp.push_back(vertex(t, Teta::facetVertex[lv][2]));
				c = nextTetCornerAroundVertex_[c];
			}
			std::sort(temp.begin(), temp.end(), std::greater<VertexIndex>());
			auto it = std::unique(temp.begin(), temp.end());
			temp.resize(std::distance(temp.begin(), it));
			std::swap(temp, vertexToVertices_[v]);
		}
	}

	void TetMeshForCombining::fillFacets()	{
		facets_.resize(nbVertices());

		unsigned int nb = 0;
		for (TetIndex t = 0; t < nbTets(); ++t) {
			for (TetFacetIndex f = 0; f < 4; ++f) {
				TetIndex adjacent = adjacentTet(t, f);
				if (t < adjacent) {
					trindex facet = tetFacet(t, f);
					facets_[facet.indices[0]].push_back(std::make_pair(facet.indices[1], facet.indices[2]));
					nb++;
				}
			}
		}
		for (auto &v : facets_) {
			std::sort(v.begin(), v.end());
		}
		std::cout << "There are " << nb << " facets " << std::endl;
	}

	/** Maybe could be done at the same time than tet connection to avoid
	*  one loop over the tets.
	*/
	void TetMeshForCombining::connectTetsToTriangles(){
		if (nbTriangles() == 0) return;

		// 0. Allocate storage space
		adjacentTriangles_.resize(4 * nbTets(), NO_ADJACENT);

		// 1. Chain triangle corners around triangle vertices
		vector<TetCornerIndex> nextTriangleCornerAroundVertex(3 * nbTriangles(), NO_ID);
		vector<TetCornerIndex> vertexToTriangleCorner(nbVertices(), NO_ID);
		for (TriangleIndex t = 0; t < nbTriangles(); ++t) {
			for (TriangleVertexIndex lv = 0; lv < 3; ++lv) {
				TriangleCornerIndex c = 3 * t + lv;
				VertexIndex v = triangleVertex(t, lv);
				nextTriangleCornerAroundVertex[c] = vertexToTriangleCorner[v];
				vertexToTriangleCorner[v] = c;
			}
		}

		// 2. Set adjacencies between tets and triangles
		for (TetIndex tet = 0; tet < nbTets(); ++tet) {
			for (TetFacetIndex f = 0; f < 4; ++f) {
				TriangleIndex adjacentTriangle = adjacentTriangles_[4 * tet + f];
				if (adjacentTriangle == NO_ADJACENT) {
					trindex facet = tetFacet(tet, f);
					VertexIndex v0 = facet.indices[0];
					for (TriangleCornerIndex c = vertexToTriangleCorner[v0]; c != NO_ID; c = nextTriangleCornerAroundVertex[c]) {
						TriangleIndex tri = c / 3;
						trindex triFacet = triangleFacet(tri);
						if (facet == triFacet) {
							adjacentTriangles_[4 * tet + f] = tri;
							break;
						}
					}
				}
			}
		}
	}
}
