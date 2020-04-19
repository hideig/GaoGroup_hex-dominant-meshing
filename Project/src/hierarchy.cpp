#include "hierarchy.h"
#include "meshio.h"
#include "timer.h"
#include "quat.h"
#include "bvh.h"
#include "orient_triangle_mesh.h"
#include "dedge.h"
#include "subdivide.h"
#include "tri_tri_intersection.h"
#include "tet_mesh.h"

MultiResolutionHierarchy::MultiResolutionHierarchy() {
    mV = { MatrixXf::Zero(3, 1) };
    mN = { MatrixXf::Zero(3, 1) };
    mO = { MatrixXf::Zero(3, 1) };
	my_mO = { MatrixXf::Zero(3, 1) };
    mQ = { MatrixXf::Zero(4, 1) };
    mBVH = nullptr;
	ratio_scale = 3.0;
	tElen_ratio = 1.0;
	tet_elen = 0.6;
	re_color = true;
	Qquadric = splitting = decomposes = doublets = triangles = false;
}
//bool MultiResolutionHierarchy::myLoad(const TetMeshForCombining &tets){
//	std::lock_guard<ordered_lock> lock(mMutex);
//
//	mV.resize(1);
//	mV[0] = MatrixXf::Zero(3, 1);
//	mF = MatrixXu::Zero(3, 1);
//
//	myLoadTetMesh(mV[0], mF, mT, tets);
//
//	mV.resize(1);
//	mAABB = AABB(
//		mV[0].rowwise().minCoeff(),
//		mV[0].rowwise().maxCoeff()
//	);
//
//	ms = compute_mesh_stats(mF, mV[0]);
//	diagonalLen = 3 * (mAABB.max - mAABB.min).norm() / 100;
//	ratio_scale = ms.mAverageEdgeLength * 3.5 / diagonalLen;
//	tet_elen = tElen_ratio * ratio_scale * diagonalLen * 0.3;
//
//	return true;
//}

bool MultiResolutionHierarchy::load(const std::string &filename) {
    std::lock_guard<ordered_lock> lock(mMutex);

    mV.resize(1);
	mV[0] = MatrixXf::Zero(3, 1);
	mF = MatrixXu::Zero(3, 1);
	
	try {
		load_obj(filename, mF, mV[0]);
		//myLoadTetMesh(filename, mV[0], mF, mT);
	}
	catch (const std::exception &e) {
		std::cout << "failed loading obj file." << std::endl;
#ifdef T_VTAG
		std::vector<std::vector<uint32_t>> mFs2D_O;
		load_off(filename, mFs2D_O, mV[0]);
		std::vector<tuple_E> mEs_O;
		std::vector<std::vector<uint32_t>> mFes_O;
		construct_tEs_tFEs(mFs2D_O, mFes_O, mEs_O);
		tagging_singularities_T_nodes(mV[0], mEs_O, mFs2D_O);
		char path[300];
		sprintf(path, "%s%s", filename.c_str(), "_V_flag.txt");
		write_Vertex_Types_TXT(V_flag, path);

#endif
		return false;
	}

	mV.resize(1);
	mAABB = AABB(
		mV[0].rowwise().minCoeff(),
		mV[0].rowwise().maxCoeff()
	);

	ms = compute_mesh_stats(mF, mV[0]); // �����mAverageEdgeLength
	diagonalLen = 3 * (mAABB.max - mAABB.min).norm() / 100;

	std::cout << "mAverageEdgeLength: " << mAverageEdgeLength << std::endl;
	ratio_scale = ms.mAverageEdgeLength * 3.5 / diagonalLen;
	std::cout << "ratio_scale: " << ratio_scale << std::endl;
	tet_elen = tElen_ratio * ratio_scale * diagonalLen * 0.3;
	std::cout << "tet_elen: " << tet_elen << std::endl;

    return true;
}
bool MultiResolutionHierarchy::my_load(const std::string &filename) {
	std::lock_guard<ordered_lock> lock(mMutex);

	mV.resize(1);
	mV[0] = MatrixXf::Zero(3, 1);
	mF = MatrixXu::Zero(3, 1);

	try {
		//load_obj(filename, mF, mV[0]);
		myLoadTetMesh(filename, mV[0], mF, mT);
	}
	catch (const std::exception &e) {
		std::cout << "failed loading obj file." << std::endl;
#ifdef T_VTAG
		std::vector<std::vector<uint32_t>> mFs2D_O;
		load_off(filename, mFs2D_O, mV[0]);
		std::vector<tuple_E> mEs_O;
		std::vector<std::vector<uint32_t>> mFes_O;
		construct_tEs_tFEs(mFs2D_O, mFes_O, mEs_O);
		tagging_singularities_T_nodes(mV[0], mEs_O, mFs2D_O);
		char path[300];
		sprintf(path, "%s%s", filename.c_str(), "_V_flag.txt");
		write_Vertex_Types_TXT(V_flag, path);

#endif
		return false;
	}

	mV.resize(1);
	mAABB = AABB(
		mV[0].rowwise().minCoeff(),
		mV[0].rowwise().maxCoeff()
	);

	ms = compute_mesh_stats(mF, mV[0]); // �����mAverageEdgeLength
	diagonalLen = 3 * (mAABB.max - mAABB.min).norm() / 100;

	std::cout << "mAverageEdgeLength: " << mAverageEdgeLength << std::endl;
	ratio_scale = ms.mAverageEdgeLength * 3.5 / diagonalLen;
	std::cout << "ratio_scale: " << ratio_scale << std::endl;
	tet_elen = tElen_ratio * ratio_scale * diagonalLen * 0.3;
	std::cout << "tet_elen: " << tet_elen << std::endl;

	return true;
}
MeshStats MultiResolutionHierarchy::compute_mesh_stats(const MatrixXu &F_, const MatrixXf &V_, bool deterministic)
{
	MeshStats stats;
	cout << "Computing mesh statistics .. ";
	cout.flush();
	auto map = [&](const tbb::blocked_range<uint32_t> &range, MeshStats stats) -> MeshStats {
		for (uint32_t f = range.begin(); f != range.end(); ++f) {
			Vector3f v[3] = { V_.col(F_(0, f)), V_.col(F_(1, f)), V_.col(F_(2, f)) };
			Vector3f face_center = Vector3f::Zero();

			for (int i = 0; i<3; ++i) {
				Float edge_length = (v[i] - v[i == 2 ? 0 : (i + 1)]).norm();
				stats.mAverageEdgeLength += edge_length;
				stats.mMaximumEdgeLength = std::max(stats.mMaximumEdgeLength, (double)edge_length);
				stats.mAABB.expandBy(v[i]);
				face_center += v[i];
			}
			face_center *= 1.0f / 3.0f;

			Float face_area = 0.5f * (v[1] - v[0]).cross(v[2] - v[0]).norm();
			stats.mSurfaceArea += face_area;
			stats.mWeightedCenter += face_area * face_center;
		}
		return stats;
	};

	auto reduce = [](MeshStats s0, MeshStats s1) -> MeshStats {
		MeshStats result;
		result.mSurfaceArea = s0.mSurfaceArea + s1.mSurfaceArea;
		result.mWeightedCenter = s0.mWeightedCenter + s1.mWeightedCenter;
		result.mAverageEdgeLength =
			s0.mAverageEdgeLength + s1.mAverageEdgeLength;
		result.mMaximumEdgeLength =
			std::max(s0.mMaximumEdgeLength, s1.mMaximumEdgeLength);
		result.mAABB = AABB::merge(s0.mAABB, s1.mAABB);
		return result;
	};

	tbb::blocked_range<uint32_t> range(0u, (uint32_t)F_.cols(), GRAIN_SIZE);

	if (deterministic)
		stats = tbb::parallel_deterministic_reduce(range, MeshStats(), map, reduce);
	else
		stats = tbb::parallel_reduce(range, MeshStats(), map, reduce);

	stats.mAverageEdgeLength /= F_.cols() * 3;
	stats.mWeightedCenter /= stats.mSurfaceArea;

	return stats;
}
bool MultiResolutionHierarchy::tet_meshing()
{
	MatrixXf V;
	MatrixXu F;
	if(!tetMesh()){
		V = mV[0];
		F = mF;
	}else {
		V = mV[0];
		F = mF;
	}
	Vector3f minV = mV[0].rowwise().minCoeff() + Vector3f(-0.1, -0.1, -0.1),
		maxV = mV[0].rowwise().maxCoeff() + Vector3f(0.1, 0.1, 0.1);

	MatrixXf Vs(3, 8);
	for(uint32_t i=0;i<8;i++)
		for (uint32_t j = 0; j < 3; j++) {
			short bit = ((1 << j) & i) >> j;
			Vs(j, i) = (bit *minV[j] + (1 - bit)*maxV[j]);
		}
	MatrixXu tris_Cube(12, 3);
	tris_Cube <<
		0, 2, 3,
		0, 3, 1,
		3, 2, 7,
		7, 2, 6,
		3, 7, 5,
		3, 5, 1,
		1, 5, 0,
		0, 5, 4,
		4, 5, 7,
		4, 7, 6,
		4, 6, 2,
		4, 2, 0;
	tris_Cube.transposeInPlace();
	orient_triangle_mesh_index(Vs, tris_Cube);

	tetgenio in_bg_, in, addin, in_bg, out_, out;

	in_bg_.numberofpoints = 8;
	in_bg_.pointlist = new REAL[8 * 3];
	for (uint32_t i = 0; i < 8; i++)
		for (uint32_t j = 0; j < 3; j++)
			in_bg_.pointlist[3 * i + j] = Vs(j, i);
	in_bg_.numberoffacets = 12;
	in_bg_.facetlist = new tetgenio::facet[12];
	in_bg_.facetmarkerlist = new int[in_bg_.numberoffacets];
	tetgenio::facet *f0;
	tetgenio::polygon *p0;
	for (uint32_t i = 0; i < in_bg_.numberoffacets; i++) {
		f0 = &in_bg_.facetlist[i];
		f0->numberofpolygons = 1;
		f0->polygonlist = new tetgenio::polygon[f0->numberofpolygons];
		f0->numberofholes = 0;
		f0->holelist = NULL;
		p0 = &f0->polygonlist[0];
		p0->numberofvertices = 3;
		p0->vertexlist = new int[p0->numberofvertices];
		p0->vertexlist[0] = tris_Cube(0, i);
		p0->vertexlist[1] = tris_Cube(1, i);
		p0->vertexlist[2] = tris_Cube(2, i);
		in_bg_.facetmarkerlist[i] = 1;
	}
	tetrahedralize("pq", &in_bg_, &out_);

	in_bg.numberofpoints = out_.numberofpoints;
	in_bg.pointlist = new REAL[in_bg.numberofpoints * 3];
	for (uint32_t i = 0; i < 3 * out_.numberofpoints; i++)
		in_bg.pointlist[i] = out_.pointlist[i];
	in_bg.numberoftetrahedra = out_.numberoftetrahedra;
	in_bg.tetrahedronlist = new int[out_.numberoftetrahedra * 4];
	for (uint32_t i = 0; i < 4 * out_.numberoftetrahedra; i++)
		in_bg.tetrahedronlist[i] = out_.tetrahedronlist[i];
	in_bg.pointmtrlist = new double[in_bg.numberofpoints];
	in_bg.numberofpointmtrs = 1;
	std::cout << "target tet edge len: " << tet_elen << endl;
	for (int i = 0; i<in_bg.numberofpoints; i++)
		in_bg.pointmtrlist[i] = tet_elen;
	
	tetgenio::facet *f;
	tetgenio::polygon *p;

	in.firstnumber = 0;

	in.numberofpoints = V.cols();
	in.pointlist = new REAL[in.numberofpoints * 3];
	in.pointmarkerlist = new int[in.numberofpoints];
	in.pointmtrlist = new double[in.numberofpoints];
	for (int i = 0; i<in.numberofpoints; i++){
		in.pointlist[3 * i + 0] = V(0,i);
		in.pointlist[3 * i + 1] = V(1,i);
		in.pointlist[3 * i + 2] = V(2,i);
		in.pointmarkerlist[i] = 1;
	}

	in.numberoffacets = F.cols();
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for (int i = 0; i<in.numberoffacets; i++)
	{
		f = &in.facetlist[i];
		f->numberofpolygons = 1;
		f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
		f->numberofholes = 0;
		f->holelist = NULL;
		p = &f->polygonlist[0];
		p->numberofvertices = 3;
		p->vertexlist = new int[p->numberofvertices];
		p->vertexlist[0] = F(0,i);
		p->vertexlist[1] = F(1,i);
		p->vertexlist[2] = F(2,i);

		in.facetmarkerlist[i] = 1;
	}
	tetrahedralize("pqm", &in, &out, &addin, &in_bg);
	mV[0].setZero(); mV[0].resize(3, out.numberofpoints);
	for (uint32_t i = 0; i < out.numberofpoints; i++) {
		for (uint32_t j = 0; j < 3; j++)
			mV[0](j, i) = out.pointlist[3 * i + j];
	}
	mT.setZero(); 
	mT.resize(4, out.numberoftetrahedra);
	for (int i = 0; i < out.numberoftetrahedra; i++){
		for (uint32_t j = 0; j < 4; j++)
			mT(j, i) = out.tetrahedronlist[4 * i + j];
	}
	//Fs
	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, bool>> tempF;
	tempF.reserve(mT.cols() * 4);
	std::vector<Vector3u> Fs;
	for (uint32_t t = 0; t < mT.cols(); ++t) {
		for (uint32_t f = 0; f < 4; ++f) {
			uint32_t v0 = mT(tet_faces[f][0], t), v1 = mT(tet_faces[f][1], t), v2 = mT(tet_faces[f][2], t);
			if (v0 > v1) std::swap(v0, v1);
			if (v1 > v2) std::swap(v2, v1);
			if (v0 > v1) std::swap(v0, v1);
			tempF.push_back(std::make_tuple(v0, v1, v2, t, true));
		}
	}
	std::sort(tempF.begin(), tempF.end());
	Fs.clear();
	Fs.reserve(tempF.size() / 3);
	int F_num = -1, f_b = 0; std::vector<bool> f_boundary; f_boundary.reserve(tempF.size() / 2);
	for (uint32_t i = 0; i < tempF.size(); ++i) {
		if (i == 0 || (i != 0 &&
			(std::get<0>(tempF[i]) != std::get<0>(tempF[i - 1]) ||
				std::get<1>(tempF[i]) != std::get<1>(tempF[i - 1]) ||
				std::get<2>(tempF[i]) != std::get<2>(tempF[i - 1])))) {
			F_num++;
			Vector3u v(std::get<0>(tempF[i]), std::get<1>(tempF[i]), std::get<2>(tempF[i]));
			Fs.push_back(v);
			f_boundary.push_back(true);
			f_b++;
		}else {
			f_boundary[F_num] = false;
			f_b--;
		}
	}
	vector<vector<uint32_t>> HF(f_b);
	mF.resize(3, f_b); f_b = 0;
	for (uint32_t f = 0; f < F_num+1; f++)
		if (f_boundary[f]) {
			mF.col(f_b) = Fs[f];
			HF[f_b].push_back(Fs[f][0]);
			HF[f_b].push_back(Fs[f][1]);
			HF[f_b++].push_back(Fs[f][2]);
		}

	orient_polygon_mesh(mV[0], HF);
	for (uint32_t i = 0; i < HF.size(); i++)
		for (uint32_t j = 0; j < 3; j++)mF(j, i) = HF[i][j];
	//orient_triangle_mesh_index(mV[0], mF);

	std::cout << "V, F, T: " << mV[0].cols() << " " << mF.cols() << " " << mT.cols() << endl;
	return true;
}
void MultiResolutionHierarchy::build() {
	Timer<> timer;
	mV.resize(1);

	if (mBVH)
		delete mBVH;

	timer.beginStage("Computing face and vertex normals");
	mN.resize(1);
	mC.resize(1);
	nV_boundary_flag.resize(1);
	mN[0].setZero(3, mV[0].cols());
	mC[0].setZero(3, mV[0].cols());
	mNF.resize(3, mF.cols()); mCF.resize(3, mF.cols());
	VectorXi count(mV[0].cols());
	count.setZero();
	for (uint32_t i = 0; i<mF.cols(); ++i) {
		uint32_t i0 = mF(0, i), i1 = mF(1, i), i2 = mF(2, i);
		Vector3f v0 = mV[0].col(i0), v1 = mV[0].col(i1), v2 = mV[0].col(i2);
		Vector3f n = (v1 - v0).cross(v2 - v0).normalized();
		mNF.col(i) = n;
		mCF.col(i) += (v0 + v1 + v2) / 3;
		mN[0].col(i0) += n; mN[0].col(i1) += n; mN[0].col(i2) += n;
		count[i0]++; count[i1]++; count[i2]++;
	}

	for (uint32_t i = 0; i<mN[0].cols(); ++i) {
		if (mN[0].col(i) != Vector3f::Zero()) {
			Vector3f d1 = mN[0].col(i) / count[i],
				d2 = mN[0].col(i).normalized();
			if (tetMesh() && d1.dot(d2) < 0.85f)
				d2 = Vector3f::Zero();;//
			mN[0].col(i) = d2;
			if (d2 != Vector3f::Zero())
				mC[0].col(i) = mV[0].col(i);
		}
	}

	vnfs.clear();// ���������
	vnfs.resize(mV[0].cols());
	for (uint32_t i = 0; i < mF.cols(); ++i) for (uint32_t j = 0; j < 3; j++) vnfs[mF(j, i)].push_back(i);

	timer.endStage();

	timer.beginStage("Computing adjacency data structure");
	mL.clear(); mP.clear();
	nV_boundary_flag[0].clear(); nV_boundary_flag[0].resize(mV[0].cols(), false);

	if (tetMesh()) {
		for (uint32_t i = 0; i < mF.cols(); ++i) for (uint32_t j = 0; j < 3; j++) nV_boundary_flag[0][mF(j, i)] = true;

		std::vector<std::pair<uint32_t, uint32_t>> adj;
		adj.reserve(mT.cols() * 12);
		for (uint32_t t = 0; t < mT.cols(); ++t) {
			const int tet_edges[6][2] = {
				{ 0, 1 },{ 0, 2 },{ 0, 3 },{ 1, 2 },{ 1, 3 },{ 2, 3 } };
			for (int i = 0; i < 6; ++i) {
				uint32_t v0 = mT(tet_edges[i][0], t);
				uint32_t v1 = mT(tet_edges[i][1], t);
				adj.push_back(std::make_pair(v0, v1));
				adj.push_back(std::make_pair(v1, v0));
			}
		}
		std::sort(adj.begin(), adj.end());
		adj.erase(std::unique(adj.begin(), adj.end()), adj.end());

		std::vector<Triplet> triplets;
		for (auto item : adj)
			triplets.push_back(Triplet(item.first, item.second, 1.f));
		mL.resize(1);
		mL[0].resize(mV[0].cols(), mV[0].cols());
		mL[0].setFromTriplets(triplets.begin(), triplets.end());
	}else {
		construct_tEs_tFEs(mF, nFes, nEs);
		//nV_nes, tag boundary V
		nV_nes.clear(); nV_nes.resize(mV[0].cols());
		for (uint32_t i = 0; i < nEs.size(); i++) {
			uint32_t v0 = std::get<0>(nEs[i]);
			uint32_t v1 = std::get<1>(nEs[i]);
			nV_nes[v0].push_back(i);
			nV_nes[v1].push_back(i);
			if (std::get<2>(nEs[i])) {
				nV_boundary_flag[0][v0] = nV_boundary_flag[0][v1] = true;
			}
		}

		std::vector<std::pair<uint32_t, uint32_t>> adj;
		adj.reserve(mF.cols() * 6);
		for (uint32_t f = 0; f < mF.cols(); ++f) {
			for (int i = 0; i < 3; ++i) {
				uint32_t v0 = mF(i, f);
				uint32_t v1 = mF((i + 1) % 3, f);
				adj.push_back(std::make_pair(v0, v1));
				adj.push_back(std::make_pair(v1, v0));
			}
		}
		std::sort(adj.begin(), adj.end());
		adj.erase(std::unique(adj.begin(), adj.end()), adj.end());

		std::vector<Triplet> triplets;
		for (auto item : adj)
			triplets.push_back(Triplet(item.first, item.second, 1.f));
		mL.resize(1);
		mL[0].resize(mV[0].cols(), mV[0].cols());
		mL[0].setFromTriplets(triplets.begin(), triplets.end());
	}
	for (uint32_t i = 0; i < (uint32_t)mL[0].rows(); ++i) {
		Float sum = 1 / mL[0].row(i).sum();
		mL[0].row(i) *= sum;
		mL[0].coeffRef(i, i) = -sum;
	}
	mL[0].makeCompressed();
	timer.endStage();

	struct WeightedEdge {
		WeightedEdge(uint32_t _i0, uint32_t _i1, Float weight)
			: weight(weight), i0(_i0), i1(_i1) {
			if (i0 > i1)
				std::swap(i0, i1);
		}

		bool operator<(const WeightedEdge &e) const {
			return std::tie(weight, i0, i1) < std::tie(e.weight, e.i0, e.i1);
		}

		Float weight;
		uint32_t i0, i1;
	};

	timer.beginStage("Building hierarchy");
	while (mL[mL.size() - 1].cols() > 1) {
		const MatrixXf &V = mV[mV.size() - 1];
		const MatrixXf &N = mN[mN.size() - 1];
		const MatrixXf &C = mC[mC.size() - 1];
		const vector<bool> &VB = nV_boundary_flag[nV_boundary_flag.size() - 1];
		const SMatrix &L = mL[mL.size() - 1];
		std::vector<bool> collapsed(L.cols(), false);
		std::vector<bool> visited(L.cols(), false);
		std::set<WeightedEdge> edges;

		double edgeSum = 0;
		size_t edgeCount = 0;
		for (int k = 0; k < L.outerSize(); ++k) {
			for (SMatrix::InnerIterator it(L, k); it; ++it) {
				if (it.col() == it.row())
					continue;
				Float length = (V.col(it.row()) - V.col(it.col())).norm();
				edgeSum += length;
				edgeCount += 1;
				edges.insert(WeightedEdge(it.row(), it.col(), length));
			}
		}
		if (mL.size() == 1)
			mAverageEdgeLength = edgeSum / edgeCount;

		std::vector<Triplet> P_triplets, R_triplets;
		std::vector<Vector3f> V_next, N_next, C_next;
		std::map<uint32_t, uint32_t> vertex_map;

		uint32_t nVertices = 0; vector<bool> vb_flag(V.cols(), false);
		for (auto const &e : edges) {
			visited[e.i0] = visited[e.i1] = true;
			if (collapsed[e.i0] || collapsed[e.i1])
				continue;
			collapsed[e.i0] = true;
			collapsed[e.i1] = true;
			P_triplets.push_back(Triplet(e.i0, nVertices, 1.0f));
			P_triplets.push_back(Triplet(e.i1, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, e.i0, 0.5f));
			R_triplets.push_back(Triplet(nVertices, e.i1, 0.5f));
			V_next.push_back(0.5f * (V.col(e.i0) + V.col(e.i1)));

			if (VB[e.i0] || VB[e.i1]) vb_flag[nVertices] = true;

			Vector3f n = N.col(e.i0) + N.col(e.i1);
			Vector3f c = C.col(e.i0) + C.col(e.i1);
			if (N.col(e.i0) != Vector3f::Zero() &&
				N.col(e.i1) != Vector3f::Zero()) {
				if (tetMesh()) {
					n = N.col(e.i0);
					c = C.col(e.i0);
				}else {
					n.normalize();
					c *= 0.5f;
				}
			}

			N_next.push_back(n);
			C_next.push_back(c);

			vertex_map[e.i0] = nVertices;
			vertex_map[e.i1] = nVertices;
			nVertices++;
		}

		for (uint32_t i = 0; i<V.cols(); ++i) {
			if (collapsed[i] || !visited[i])
				continue;
			P_triplets.push_back(Triplet(i, nVertices, 1.0f));
			R_triplets.push_back(Triplet(nVertices, i, 1.0f));
			V_next.push_back(V.col(i));
			N_next.push_back(N.col(i));
			C_next.push_back(C.col(i));
			vertex_map[i] = nVertices;

			if (VB[i]) vb_flag[nVertices] = true;

			nVertices++;
		}
		vb_flag.resize(nVertices);

		if (mL.size() != 1)
			std::cout << ", ";
		std::cout << nVertices;
		std::cout.flush();

		SMatrix P(V.cols(), nVertices), R(nVertices, V.cols());

		P.setFromTriplets(P_triplets.begin(), P_triplets.end());
		R.setFromTriplets(R_triplets.begin(), R_triplets.end());

		SMatrix L2 = R*L*P;
		MatrixXf V2(3, nVertices), N2(3, nVertices), C2(3, nVertices), Q2(4, nVertices);
		for (uint32_t i = 0; i<nVertices; ++i) {
			V2.col(i) = V_next[i];
			N2.col(i) = N_next[i];
			C2.col(i) = C_next[i];
		}

		nV_boundary_flag.push_back(vb_flag);
		mP.push_back(std::move(P));
		mN.push_back(std::move(N2));
		mV.push_back(std::move(V2));
		mC.push_back(std::move(C2));
		mL.push_back(L2);
	}
	std::cout << " ";
	timer.endStage();

	mQ.resize(mL.size());
	mO.resize(mL.size());
	my_mO.resize(mL.size());

	pcg32 rng;
	if (tetMesh()) {
		for (uint32_t i = 0; i < mL.size(); ++i) {
			mQ[i].resize(4, mV[i].cols());
			mO[i].resize(3, mV[i].cols());
			my_mO[i].resize(3, mV[i].cols());
			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				mQ[i].col(j) = Quaternion::Random(rng);

				mO[i].col(j) = aabbRand(mAABB, rng);
			}
		}
	}else {
		for (uint32_t i = 0; i < mL.size(); ++i) {
			mQ[i].resize(3, mV[i].cols());
			mO[i].resize(3, mV[i].cols());
			my_mO[i].resize(3, mV[i].cols());

			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				if (i == 0 && nV_boundary_flag[i][j]) {
					std::vector<uint32_t> vs;
					for (auto eid : nV_nes[j]) {
						if (!std::get<2>(nEs[eid])) continue;

						uint32_t v0 = std::get<0>(nEs[eid]);
						uint32_t v1 = std::get<1>(nEs[eid]);
						if (v0 == j) vs.push_back(v1);
						else vs.push_back(v0);
					}
					Vector3f direct0, direct1; direct1.setZero();
					direct0 = (mV[0].col(j) - mV[0].col(vs[0])).normalized();
					for (uint32_t k = 1; k < vs.size(); k++) {
						Vector3f direct_ = (mV[0].col(vs[k]) - mV[0].col(j)).normalized();
						direct1 += direct_;
					}

					if (std::abs(direct0.dot(direct1)) < 0.5)
						mQ[i].col(j) = direct0;
					else if (direct0 + direct1 == Vector3f::Zero())
						mQ[i].col(j) = direct0;
					else
						mQ[i].col(j) = (direct0 + direct1).normalized();
					continue;
				}

				Vector3f n = mN[i].col(j), v = mV[i].col(j);
				Vector3f s, t;
				coordinate_system(n, s, t);
				float angle = rng.nextFloat() * 2 * M_PI;
				mQ[i].col(j) = s * std::cos(angle) + t * std::sin(angle);
			}


			for (uint32_t j = 0; j < mV[i].cols(); ++j) {
				Vector3f n = mN[i].col(j), v = mV[i].col(j);
				rng.nextFloat();
				Vector3f o = aabbRand(mAABB, rng);
				o -= n.dot(o - v) * n;
				mO[i].col(j) = o;
			}
		}
		//propagate up
		for (uint32_t i = 1; i < mL.size(); ++i) {
			for (int k = 0; k < mP[i - 1].outerSize(); ++k) {
				SMatrix::InnerIterator it(mP[i - 1], k);
				for (; it; ++it) {
					if (nV_boundary_flag[i - 1][it.row()])
						mQ[i].col(it.col()) = mQ[i - 1].col(it.row());
				}
			}
		}
	}
	mOrientationIterations = 0;
	mPositionIterations = 0;

	mBVH = new BVH(&mF, &mV[0], mAABB);
	mBVH->build();
	mScale = tet_elen;
	//mScale = diagonalLen * ratio_scale;
	mInvScale = 1.f / mScale;

	std::cout << "diagonalLen, ratio_scale, mScale: " << diagonalLen << ", " << ratio_scale << ", " << mScale << std::endl;

	sta.tN = mF.cols();
	sta.tetN = mT.cols();
}
void MultiResolutionHierarchy::construct_tEs_tFEs(MatrixXu & F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs) {
	mtFes.clear(); mtEs.clear();

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(F.cols() * 3);
	mtFes.resize(F.cols());
	for (uint32_t f = 0; f < F.cols(); ++f) {
		for (uint32_t e = 0; e < 3; ++e) {
			uint32_t v0 = F(e, f), v1 = F((e + 1) % 3, f);
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(3);
		mtFes[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mtEs.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) || std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mtEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(mtEs[E_num]) = false;

		mtFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
}
void MultiResolutionHierarchy::construct_tEs_tFEs(std::vector<std::vector<uint32_t>> &F, std::vector<std::vector<uint32_t>> &mtFes, std::vector<tuple_E> &mtEs) {
	mtFes.clear(); mtEs.clear();

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(F.size() * 3);
	mtFes.resize(F.size());
	for (uint32_t f = 0; f < F.size(); ++f) {
		for (uint32_t e = 0; e < F[f].size(); ++e) {
			uint32_t v0 = F[f][e], v1 = F[f][(e+1)%F[f].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(F[f].size());
		mtFes[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	mtEs.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) || std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			mtEs.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(mtEs[E_num]) = false;

		mtFes[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
}
void MultiResolutionHierarchy::orient_polygon_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF, vector<vector<uint32_t>> &HFE, vector<tuple_E> &Es) {
	//orient surface
	if (!HF.size())return;
	uint32_t start_f = 0;
	vector<bool> flag(HF.size(), true);
	vector<vector<uint32_t>> Efs(Es.size());
	for (uint32_t i = 0; i < HFE.size();i++)for (auto e : HFE[i])Efs[e].push_back(i);
	flag[start_f] = false;
	std::queue<uint32_t> pf_temp; pf_temp.push(start_f);
	while (!pf_temp.empty()) {
		uint32_t fid = pf_temp.front(); pf_temp.pop();
		for (auto eid : HFE[fid]) for (auto nfid : Efs[eid]) {
			if (!flag[nfid]) continue;
			uint32_t v0 = std::get<0>(Es[eid]), v1 = std::get<1>(Es[eid]);
			int32_t v0_pos = std::find(HF[fid].begin(), HF[fid].end(), v0) - HF[fid].begin();
			int32_t v1_pos = std::find(HF[fid].begin(), HF[fid].end(), v1) - HF[fid].begin();

			if ((v0_pos + 1) % HF[fid].size() != v1_pos) swap(v0, v1);

			int32_t v0_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v0) - HF[nfid].begin();
			int32_t v1_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v1) - HF[nfid].begin();

			if ((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) std::reverse(HF[nfid].begin(), HF[nfid].end());

			pf_temp.push(nfid); flag[nfid] = false;
		}
		if (pf_temp.empty()) {
			bool found = false;
			for (uint32_t i = 0; i < flag.size(); i++)if (flag[i]) {
				start_f = i;
				flag[start_f] = false; pf_temp.push(start_f);
				found = true;
			}
			if (!found) break;
		}
	}

	Float res = 0;
	Vector3f ori; ori.setZero();
	for (uint32_t i = 0; i < HF.size(); i++) {
		auto &fvs = HF[i];
		Vector3f center; center.setZero(); for (auto vid : fvs) center += HV.col(vid); center /= fvs.size();

		for (uint32_t j = 0; j < fvs.size(); j++) {
			Vector3f x = HV.col(fvs[j]) - ori, y = HV.col(fvs[(j + 1) % fvs.size()]) - ori, z = center - ori;
			res += -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
		}
	}
	if (res > 0) {
		for (uint32_t i = 0; i < HF.size(); i++) std::reverse(HF[i].begin(), HF[i].end());
	}
}
void MultiResolutionHierarchy::orient_polygon_mesh(MatrixXf &HV, vector<vector<uint32_t>> &HF) {
	//orient surface
	if (!HF.size())return;

	vector<vector<uint32_t>> HFE(HF.size());
	vector<tuple_E> Es;

	std::vector<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t, int>> temp;
	temp.reserve(HF.size() * 3);
	HFE.resize(HF.size());
	for (uint32_t f = 0; f < HF.size(); ++f) {
		for (uint32_t e = 0; e < HF[f].size(); ++e) {
			uint32_t v0 = HF[f][e], v1 = HF[f][(e + 1) % HF[f].size()];
			if (v0 > v1) std::swap(v0, v1);
			temp.push_back(std::make_tuple(v0, v1, f, e, Edge_tag::B));
		}
		std::vector<uint32_t> fes(HF[f].size());
		HFE[f] = fes;
	}
	std::sort(temp.begin(), temp.end());
	Es.reserve(temp.size() / 2);
	int E_num = -1;
	for (uint32_t i = 0; i < temp.size(); ++i) {
		if (i == 0 || (i != 0 && (std::get<0>(temp[i]) != std::get<0>(temp[i - 1]) ||
			std::get<1>(temp[i]) != std::get<1>(temp[i - 1])))) {
			E_num++;
			Es.push_back(std::make_tuple(std::get<0>(temp[i]), std::get<1>(temp[i]), true, 0, std::get<4>(temp[i]), E_num, -1, 0));
		}
		else if (i != 0 && (std::get<0>(temp[i]) == std::get<0>(temp[i - 1]) &&
			std::get<1>(temp[i]) == std::get<1>(temp[i - 1])))
			std::get<2>(Es[E_num]) = false;
		HFE[std::get<2>(temp[i])][std::get<3>(temp[i])] = E_num;
	}
	vector<vector<uint32_t>> Efs(Es.size());
	for (uint32_t i = 0; i < HFE.size(); i++)for (auto e : HFE[i])Efs[e].push_back(i);

	uint32_t start_f = 0;
	vector<bool> flag(HF.size(), true);
	flag[start_f] = false;
	std::queue<uint32_t> pf_temp; pf_temp.push(start_f);
	while (!pf_temp.empty()) {
		uint32_t fid = pf_temp.front(); pf_temp.pop();
		for (auto eid : HFE[fid]) for (auto nfid : Efs[eid]) {
			if (!flag[nfid]) continue;
			uint32_t v0 = std::get<0>(Es[eid]), v1 = std::get<1>(Es[eid]);
			int32_t v0_pos = std::find(HF[fid].begin(), HF[fid].end(), v0) - HF[fid].begin();
			int32_t v1_pos = std::find(HF[fid].begin(), HF[fid].end(), v1) - HF[fid].begin();

			if ((v0_pos + 1) % HF[fid].size() != v1_pos) swap(v0, v1);

			int32_t v0_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v0) - HF[nfid].begin();
			int32_t v1_pos_ = std::find(HF[nfid].begin(), HF[nfid].end(), v1) - HF[nfid].begin();

			if ((v0_pos_ + 1) % HF[nfid].size() == v1_pos_) std::reverse(HF[nfid].begin(), HF[nfid].end());

			pf_temp.push(nfid); flag[nfid] = false;
		}
		if (pf_temp.empty()) {
			bool found = false;
			for (uint32_t i = 0; i < flag.size(); i++)if (flag[i]) {
				start_f = i;
				flag[start_f] = false; pf_temp.push(start_f);
				found = true;
			}
			if (!found) break;
		}
	}

	Float res = 0;
	Vector3f ori; ori.setZero();
	for (uint32_t i = 0; i < HF.size(); i++) {
		auto &fvs = HF[i];
		Vector3f center; center.setZero(); for (auto vid : fvs) center += HV.col(vid); center /= fvs.size();

		for (uint32_t j = 0; j < fvs.size(); j++) {
			Vector3f x = HV.col(fvs[j]) - ori, y = HV.col(fvs[(j + 1) % fvs.size()]) - ori, z = center - ori;
			res += -((x[0] * y[1] * z[2] + x[1] * y[2] * z[0] + x[2] * y[0] * z[1]) - (x[2] * y[1] * z[0] + x[1] * y[0] * z[2] + x[0] * y[2] * z[1]));
		}
	}
	if (res > 0) {
		for (uint32_t i = 0; i < HF.size(); i++) std::reverse(HF[i].begin(), HF[i].end());
	}
}