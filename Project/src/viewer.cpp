﻿#include "viewer.h"
#include "im_resources.h"
#include "timer.h"
#include "bvh.h"
#include <nanogui/serializer/opengl.h>
#include "hxt_combine_cpp_api.h"
#include "tet_mesh.h"
using namespace HXTCombine;
//int my_process(const MatrixXf &V_, const MatrixXu &F_, const MatrixXu &T_);
int my_process(const MultiResolutionHierarchy &mRes, HXTCombineCellStore* combineRes);
void my_show_edges(MatrixXf &Result_Vs, std::vector<tuple_E> &Es_to_render, MatrixXf &Result_edges)
{
	Result_edges.resize(6, 2 * Es_to_render.size());
	//for rendering edges
	for (uint32_t i = 0; i < Es_to_render.size(); ++i) {
		Vector3f color;
		if (std::get<4>(Es_to_render[i]) == Edge_tag::R)
			color = Vector3f(1, 0, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::B)
			color = Vector3f(0, 0, 1);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::D)
			color = Vector3f(0, 1, 0);
		else if (std::get<4>(Es_to_render[i]) == Edge_tag::H)
			color = Vector3f(1, 1, 1);

		uint32_t i0 = std::get<0>(Es_to_render[i]), i1 = std::get<1>(Es_to_render[i]);

		Result_edges.col(i * 2 + 0) << Result_Vs.col(i0), color;
		Result_edges.col(i * 2 + 1) << Result_Vs.col(i1), color;
	}
}

Viewer::Viewer(std::string &filename, bool fullscreen)
    : Screen(Vector2i(1280, 960), "Robust Quad/Hex-dominant Meshes", true),
      mOptimizer(nullptr) {

    mOptimizer = new Optimizer(mRes);

    /* Initialize shaders for rendering geometry and fields */
    mMeshShader.init("mesh_shader",
        (const char *) shader_mesh_vert,
        (const char *) shader_mesh_frag,
        (const char *) shader_mesh_geo);

    mTetShader.init("tet_shader",
        (const char *) shader_tet_vert,
        (const char *) shader_mesh_frag,
        (const char *) shader_tet_geo);

    mOrientationFieldShaderTet.init("orientation_field_shader_tet",
        (const char *) shader_orientation_field_tet_vert,
        (const char *) shader_orientation_field_tet_frag,
        (const char *) shader_orientation_field_tet_geo);

    mOrientationFieldShaderTri.init("orientation_field_shader_tri",
        (const char *) shader_orientation_field_tri_vert,
        (const char *) shader_orientation_field_tri_frag,
        (const char *) shader_orientation_field_tri_geo);

    mPositionFieldShader.init("position_field_shader",
        (const char *) shader_position_field_vert,
        (const char *) shader_position_field_frag);
	//mOtherEdge.init("other_positions",
	//	(const char *)shader_position_field_vert,
	//	(const char *)shader_position_field_frag);

    mOrientationSingularityShaderTet.init("orientation_singularity_shader_tet",
        (const char *) shader_singularity_tet_vert,
        (const char *) shader_singularity_tet_frag);

    mOrientationSingularityShaderTri.init("orientation_singularity_shader_tri",
        (const char *) shader_singularity_tri_vert,
        (const char *) shader_singularity_tri_frag,
        (const char *) shader_singularity_tri_geo);

    mPositionSingularityShaderTet.init("position_singularity_shader_tet",
        (const char *) shader_singularity_tet_vert,
        (const char *) shader_singularity_tet_frag);

    mPositionSingularityShaderTri.init("position_singularity_shader_tri",
        (const char *) shader_singularity_tri_vert,
        (const char *) shader_singularity_tri_frag,
        (const char *) shader_singularity_tri_geo);

    mExtractionResultShader.init("extraction_result",
        (const char *) shader_singularity_tet_vert,
        (const char *) shader_singularity_tet_frag);
	mExtractionResultShader2.init("extraction_result2",
		(const char *)shader_singularity_tet_vert,
		(const char *)shader_singularity_tet_frag);

	
	mEdge_color_morph2.init("mEdge_color_morph2",
		(const char *)shader_singularity_tet_vert,
			(const char *)shader_singularity_tet_frag);


	mExtractionResultShader_E_done.init("Shader_E_local",
		(const char *)shader_singularity_tet_vert,
		(const char *)shader_singularity_tet_frag);



	mExtractionResultShader_F_done.init("Shader_F_local",
		(const char *)shader_mesh_vert,
		(const char *)shader_mesh_frag,
		(const char *)shader_mesh_geo);

    /* Default view setup */
    mCamera.arcball = Arcball();
    mCamera.arcball.setSize(mSize);
    mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
    mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
    mCamera.zoom = 1.0f;
    mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
    mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
    mBaseColorBoundary = mRes.tetMesh() ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
    mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
    mSpecularColorBoundary = mRes.tetMesh() ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
    mTranslate = false;

	/* Scan over example files in the 'datasets' directory */
	auto ctx = nvgContext();
	try {
		mExampleImages = nanogui::loadImageDirectory(ctx, "resources");
	}catch (const std::runtime_error &e) {
	}
	mExampleImages.insert(mExampleImages.begin(),
		std::make_pair(nvgImageIcon(ctx, loadmesh), ""));

	/* Initialize user interface */
    Window *window = new Window(this, "");
    window->setPosition(Vector2i(15, 15));
    window->setLayout(new GroupLayout());

	PopupButton *openBtn0 = new PopupButton(window, "Open mesh");
	openBtn0->setBackgroundColor(Color(0, 255, 0, 25));
	openBtn0->setIcon(ENTYPO_ICON_FOLDER);
	Popup *popup0 = openBtn0->popup();
	VScrollPanel *vscroll = new VScrollPanel(popup0);
	ImagePanel *panel0 = new ImagePanel(vscroll);
	panel0->setImages(mExampleImages);
	panel0->setCallback([&, openBtn0](int i) {
		openBtn0->setPushed(false);

		std::string filename2 = mExampleImages[i].second;
		std::string extension;
		if (filename2.size() > 4)
			extension = str_tolower(filename2.substr(filename2.size() - 4));

		if (filename2.empty()) {
#ifdef T_VTAG
			filename2 = nanogui::file_dialog({
				{ "obj", "Wavefront OBJ" }, { "off", "OFF" }
			}, false);


#else
			filename2 = nanogui::file_dialog({
				{ "obj", "Wavefront obj" }
			}, false);
			//filename2 = nanogui::file_dialog({
			//	{ "mesh", "Wavefront mesh" }
			//}, false);
#endif
			if (filename2 == "")
				return;
		}

		if(!mRes.load(filename2)) return;
		//if(!mRes.my_load(filename2)) return;
		mScaleBox->setValue(mRes.scale());

		filename = filename2;
		mRes.outpath = filename2;
		/* Default view setup */
		mCamera.arcball = Arcball();
		mCamera.arcball.setSize(mSize);
		mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
		mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
		mCamera.zoom = 1.0f;
		mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
		mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
		//mBaseColor = Vector4f(1.0f, 1.0f, 1.0f, 1.f);
		mBaseColorBoundary = mRes.tetMesh() ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
		mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
		mSpecularColorBoundary = mRes.tetMesh() ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
		mTranslate = false;

	});

	if (mRes.mV[0].cols()) {
		/* Default view setup */
		mCamera.arcball = Arcball();
		mCamera.arcball.setSize(mSize);
		mCamera.modelTranslation = -mRes.aabb().center().cast<float>();
		mCamera.modelZoom = 3.0f / mRes.aabb().extents().cwiseAbs().maxCoeff();
		mCamera.zoom = 1.0f;
		mLightPosition = Vector3f(0.0f, 0.3f, 5.0f);
		mBaseColor = Vector4f(0.4f, 0.5f, 0.7f, 1.f);
		//mBaseColor = Vector4f(1.0f, 1.0f, 1.0f, 1.f);
		mBaseColorBoundary = mRes.tetMesh() ? Vector4f(0.0f, 0.0f, 1.0f, .2f) : mBaseColor;
		mSpecularColor = Vector4f(1.f, 1.f, 1.f, 1.f);
		mSpecularColorBoundary = mRes.tetMesh() ? Vector4f(1.f, 1.f, 1.f, .2f) : mSpecularColor;
		mTranslate = false;
	}

//Parameters
	new Label(window, "Output scale", "sans-bold");
	mScaleBox = new FloatBox<Float>(window);
	mScaleBox->setValue(mRes.scale());
	mScaleBox->setEditable(true);
	mScaleBox->setAlignment(TextBox::Alignment::Right);
	mScaleBox->setId("outputscale");

//2D&3D
	Widget *statePanel = new Widget(window);
	statePanel->setLayout(new BoxLayout(Orientation::Horizontal, nanogui::Alignment::Middle, 0, 5));
	mSolveDatastructureBtn = new Button(statePanel, "Surface", ENTYPO_ICON_FLASH);
	mSolveDatastructureBtn->setBackgroundColor(Color(0, 0, 255, 25));
	mSolveDatastructureBtn->setFlags(Button::Flags::ToggleButton);
	mSolveDatastructureBtn->setChangeCallback([&](bool value) {
		mRes.build();

		mTetShader.bind();

		MatrixXf vertexColors = MatrixXf::Zero(4, mRes.vertexCount());
		mTetShader.uploadAttrib("position", mRes.V());
		mTetShader.uploadIndices(mRes.T());
		mTetShader.uploadAttrib("color", vertexColors);

		mMeshShader.bind();
		mMeshShader.shareAttrib(mTetShader, "position");
		mMeshShader.uploadIndices(mRes.F());
		mMeshShader.uploadAttrib("normal", mRes.N());


		mOrientationFieldShaderTri.bind();
		mOrientationFieldShaderTri.shareAttrib(mTetShader, "position");
		mOrientationFieldShaderTri.uploadAttrib("q", mRes.Q());
		mOrientationFieldShaderTri.uploadAttrib("n", mRes.N());

		mPositionFieldShader.bind();
		mPositionFieldShader.uploadAttrib("o", mRes.O());
		mOtherEdge.bind();
		mOtherEdge.uploadAttrib("o", mRes.my_O());

		mLayers[Layers::Boundary]->setChecked(true);
		mLayers[Layers::PositionField]->setChecked(true);

	//	mLayers[Layers::OtherEdge]->setChecked(true);
	});

	mTmeshingBtn = new Button(statePanel, "Volume", ENTYPO_ICON_FLASH);
	mTmeshingBtn->setBackgroundColor(Color(0, 0, 255, 25));
	mTmeshingBtn->setFlags(Button::Flags::ToggleButton);
	mTmeshingBtn->setChangeCallback([&](bool value) {
		mRes.tet_meshing();
		mRes.build();
		mTetShader.bind();
		MatrixXf vertexColors = MatrixXf::Zero(4, mRes.vertexCount());
		mTetShader.uploadAttrib("position", mRes.V());
		mTetShader.uploadIndices(mRes.T());
		mTetShader.uploadAttrib("color", vertexColors);

		mMeshShader.bind();
		mMeshShader.shareAttrib(mTetShader, "position");
		mMeshShader.uploadIndices(mRes.F());
		mMeshShader.uploadAttrib("normal", mRes.N());

		mOrientationFieldShaderTet.bind();
		mOrientationFieldShaderTet.shareAttrib(mTetShader, "position");
		mOrientationFieldShaderTet.uploadAttrib("q", mRes.Q());

		mPositionFieldShader.bind();
		mPositionFieldShader.uploadAttrib("o", mRes.O());
	    mOtherEdge.bind();
		mOtherEdge.uploadAttrib("o", mRes.my_O());

		mLayers[Layers::Boundary]->setChecked(true);
		mLayers[Layers::PositionField]->setChecked(true);
		//mLayers[Layers::OtherEdge]->setChecked(true);
		/*
		std::cout << "------------------- my_process ----------------" << mRes.mScale << std::endl;
		
		MeshStore ioMesh;
		myReadFileMESH(mRes.mV[0], mRes.mF, mRes.mT, ioMesh);
		auto start0 = std::chrono::high_resolution_clock::now();

		TetMeshForCombining tets(&ioMesh);

		auto finish0 = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> t_mesh(finish0 - start0);
		std::cout << "Mesh structure built in " << t_mesh.count() << " seconds" << std::endl;

		auto start = std::chrono::high_resolution_clock::now();

		HXTCombineCellStore TheResult(tets);

		int hexFlag = 1;
		int prismFlag = 0, pyramidFlag = 0;
		double minQuality = 0.;
		if (hexFlag) {
			TheResult.computeHexes(minQuality);
			auto finish = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> t0(finish - start);
			std::cout << TheResult.hexes().size() << " potential hexes computed in " << t0.count() << " seconds" << std::endl;
		}

		if (prismFlag) {
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

		std::array<bool, 4> cellTypes{ bool(hexFlag), bool(prismFlag), bool(pyramidFlag), true };
		TheResult.selectCellsGreedyLocal(cellTypes);
		//TheResult.selectCellsGreedy(cellTypes);
		auto endSelect = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> ts(endSelect - startSelect);
		if (hexFlag) std::cout << nbTrueValues(TheResult.selectedHexes()) << " selected hexes" << std::endl;
		if (prismFlag)   std::cout << nbTrueValues(TheResult.selectedPrisms()) << " selected prisms" << std::endl;
		if (pyramidFlag) std::cout << nbTrueValues(TheResult.selectedPyramids()) << " selected pyramids" << std::endl;
		std::cout << nbTrueValues(TheResult.selectedTets()) << " tetrahedra remain" << std::endl;
		std::cout << "Timings cell selection " << ts.count() << "seconds" << std::endl;

		 //初始化 mpEs， mV_tag， F_tag;
		unsigned int num = TheResult.mesh_.nbVertices();
		mRes.mVv_tag.resize(3, num);
		for (uint32_t i = 0; i < num; ++i) {
			mRes.mVv_tag(0, i) = TheResult.mesh_.point(i)[0];
			mRes.mVv_tag(1, i) = TheResult.mesh_.point(i)[1];
			mRes.mVv_tag(2, i) = TheResult.mesh_.point(i)[2];
		}

		const TetMeshForCombining& mesh_= TheResult.mesh_;
		unsigned int edgeIndex = 0;
		//// TETS 
		
		for (unsigned int t = 0; t < mesh_.nbTets(); ++t) {
			unsigned int v0, v1, v2, v3;
			if (!((TheResult.selectedCells_[3])[t])) continue;
			else {
				v0 = mesh_.vertex(t, 0);
				v1 = mesh_.vertex(t, 1);
				v2 = mesh_.vertex(t, 2);
				v3 = mesh_.vertex(t, 3);
			}
			//v0, v1, boundary, energy, color, edge index, xy/yz/xz plane, timestamp
			mRes.mpEes.push_back(std::make_tuple(v0, v1, false, 0, Edge_tag::H, edgeIndex++, -1, 0));
			mRes.mpEes.push_back(std::make_tuple(v1, v2, false, 0, Edge_tag::H, edgeIndex++, -1, 0));
			mRes.mpEes.push_back(std::make_tuple(v2, v3, false, 0, Edge_tag::H, edgeIndex++, -1, 0));
			mRes.mpEes.push_back(std::make_tuple(v3, v0, false, 0, Edge_tag::H, edgeIndex++, -1, 0));
		}
		// OTHER CELLS
		for (unsigned int type = 0; type + 1 < cellTypes.size(); ++type) {
			if (!cellTypes[type]) continue;
			const std::vector<HXTCombineCell>& cells = TheResult.cells_[type];
			const std::vector<bool>& selected = TheResult.selectedCells_[type];
			for (unsigned int i = 0; i < cells.size(); ++i) {
				if (selected[i]) {
					unsigned int v[8];
					if (cells[i].isHex()) {
						for (unsigned int j = 0; j < 8; j++) {
							v[j] = cells[i].vertex(j);
						}
						for (unsigned int j = 0; j < 7; ++j) {
							for (unsigned int k = j + 1; k < 8; ++k) {
								bindex edge(v[j], v[k]);
								if (isEdge(cells[i], edge)) {
									mRes.mpEes.push_back(std::make_tuple(v[j], v[k], false, 0, Edge_tag::B, edgeIndex++, -1, 0));
									//mRes.mpEes.push_back(std::make_tuple(v[j], v[k], false, 0, edgeIndex % 4, edgeIndex++, -1, 0));
								}
							}
						}

						for (unsigned int j = 0; j < 16; ++j) {
							for (unsigned int j1 = j + 1; j1 < 16; ++j1) {
								for (unsigned int j2 = j1 + 1; j2 < 16; ++j2) {
									for (unsigned int j3 = j2 + 1; j3 < 16; ++j3) {
										quadindex facet(v[j%8], v[j1%8], v[j2%8], v[j3%8]);
										if (isQuadFacet(cells[i], facet)) {
											mRes.Ff_tag.push_back(std::vector<uint32_t>{v[j%8], v[j1%8], v[j2%8], v[j3%8]});
										}
									}
								}
							}
						}

					}else if (cells[i].isPrism()) {
					}else if (cells[i].isPyramid()) {
					}
				}
			}
		}
		*/

	});


//Rosy
    mSolveOrientationBtn = new Button(window, "Rosy", ENTYPO_ICON_FLASH);
    mSolveOrientationBtn->setBackgroundColor(Color(0, 0, 255, 25));
    mSolveOrientationBtn->setFlags(Button::Flags::ToggleButton);
    mSolveOrientationBtn->setChangeCallback([&](bool value) {
        mOptimizer->setOptimizeOrientations(value);
        mOptimizer->notify();

        mLayers[Layers::OrientationField]->setChecked(true);
        mLayers[Layers::OrientationSingularities]->setChecked(true);
        mLayers[Layers::PositionField]->setChecked(false);
		//mLayers[Layers::OtherEdge]->setChecked(false);
        mLayers[Layers::PositionSingularities]->setChecked(false);
        if (value == false)
            updateOrientationSingularities();

		std::cout << "_Rosy_______________mscale: " << mRes.mScale << std::endl;
    });
//Posy
	mSolvePositionBtn = new Button(window, "Posy", ENTYPO_ICON_FLASH);
    mSolvePositionBtn->setBackgroundColor(Color(0, 0, 255, 25));
    mSolvePositionBtn->setFlags(Button::Flags::ToggleButton);
    mSolvePositionBtn->setChangeCallback([&](bool value) {

        mOptimizer->setOptimizePositions(value);
        mOptimizer->notify();

        mLayers[Layers::OrientationField]->setChecked(false);
        mLayers[Layers::OrientationSingularities]->setChecked(false);
        mLayers[Layers::PositionField]->setChecked(true);
		//mLayers[Layers::OtherEdge]->setChecked(true);
        mLayers[Layers::PositionSingularities]->setChecked(true);
        if (value == false)
            updatePositionSingularities();

		std::cout << "_Posy_______________mscale: " << mRes.mScale << std::endl;
    });
//Extraction
    mExtractBtn = new Button(window, "Extract", ENTYPO_ICON_FLASH);
    mExtractBtn->setBackgroundColor(Color(0, 255, 0, 25));
    mExtractBtn->setCallback([&]() {
		if (!mRes.tetMesh()) {
			mRes.re_color = true;
			mRes.doublets = true;
			mRes.splitting = true;
			mRes.triangles = false;//true;
			mRes.decomposes = true;
			mRes.meshExtraction2D();
		}else{
			mRes.re_color = true;
			mRes.splitting = true;
			mRes.doublets = false;
			mRes.triangles = false;
			mRes.decomposes = false;
			//mRes.meshExtraction3D();
			mRes.meshExtraction3D(mRes);
		}
		mLayers[PositionSingularities]->setChecked(false);
		mLayers[Layers::PositionField]->setChecked(false);
		mLayers[Layers::Boundary]->setChecked(false);
		mLayers[OrientationSingularities]->setChecked(false);
		mLayers[Layers::OrientationField]->setChecked(false);
		//mLayers[Layers::OtherEdge]->setChecked(false);

		mExtractionResultShader_F_done.bind();
		mExtractionResultShader_F_done.uploadAttrib("position", mRes.mV_final_rend);
		mExtractionResultShader_F_done.uploadIndices(mRes.F_final_rend);
		mShow_F_done->setChecked(true);

		auto const &R4 = mRes.E_final_rend;
		mExtractionResultShader_E_done.bind();
		mExtractionResultShader_E_done.uploadAttrib("position", MatrixXf(R4.block(0, 0, 3, R4.cols())));
		mExtractionResultShader_E_done.uploadAttrib("color", MatrixXf(R4.block(3, 0, 3, R4.cols())));
		mShow_E_done->setChecked(true);

	});

	// my Combine Extraction
	mCombineBtn = new Button(window, "Combine", ENTYPO_ICON_FLASH);
	mCombineBtn->setBackgroundColor(Color(0, 255, 0, 25));
	mCombineBtn->setCallback([&]() {
		// 显示combine的结果
		std::cout << "----------------- combineShow ---------------- "  << std::endl;

		mRes.E_final_rend.setZero();
		mRes.E_final_rend.resize(6, 2 * mRes.mpEes.size());
		mRes.composit_edges_colors(mRes.mVv_tag, mRes.mpEes, mRes.E_final_rend);
	//	mRes.composit_edges_colors(mRes.mVv_tag, mRes.mpEes, mRes.E_final_rend);
	//	mRes.composit_edges_centernodes_triangles(mRes.Ff_tag, mRes.mVv_tag, mRes.E_final_rend, mRes.mV_final_rend, mRes.F_final_rend);
		
		mLayers[PositionSingularities]->setChecked(false);
		mLayers[Layers::PositionField]->setChecked(false);
		mLayers[Layers::Boundary]->setChecked(false);
		mLayers[OrientationSingularities]->setChecked(false);
		mLayers[Layers::OrientationField]->setChecked(false);
		mExtractionResultShader_F_done.bind();
		mExtractionResultShader_F_done.uploadAttrib("position", mRes.mV_final_rend);
		mExtractionResultShader_F_done.uploadIndices(mRes.F_final_rend);
		mShow_F_done->setChecked(true);
		auto const &R4 = mRes.E_final_rend;
		mExtractionResultShader_E_done.bind();
		mExtractionResultShader_E_done.uploadAttrib("position", MatrixXf(R4.block(0, 0, 3, R4.cols())));
		mExtractionResultShader_E_done.uploadAttrib("color", MatrixXf(R4.block(3, 0, 3, R4.cols())));
		mShow_E_done->setChecked(true);
	});


 	//Config Layers
	PopupButton *openBtn3 = new PopupButton(window, "Config Layers");
	auto popup3 = openBtn3->popup();
	popup3->setLayout(new GroupLayout());
	Configlayers[Config_Layers::Alignment] = new CheckBox(popup3, "Boundary alignment");
	Configlayers[Config_Layers::Extrinsic] = new CheckBox(popup3, "Extrinsic smoothing");
	Configlayers[Config_Layers::Randomization] = new CheckBox(popup3, "Randomization");
	Configlayers[Config_Layers::Hierarchy] = new CheckBox(popup3, "Hierarchy");
	int ctr = 0;
	for (auto l : Configlayers) {
		l->setChecked(true);
		l->setId("configlayer." + std::to_string(ctr++));
	}
	//Render Layers
	PopupButton *openBtn = new PopupButton(window, "Render Layers");
	auto popup = openBtn->popup();
	popup->setLayout(new GroupLayout());

	mLayers[Layers::Tetrahedra] = new CheckBox(popup, "Tetrahedra");
	mLayers[Layers::OrientationField] = new CheckBox(popup, "Orientation field");
	mLayers[Layers::OrientationSingularities] = new CheckBox(popup, "Orientation singularities");
	mLayers[Layers::PositionField] = new CheckBox(popup, "Position field");
	mLayers[Layers::PositionSingularities] = new CheckBox(popup, "Position singularities");
	mLayers[Layers::Boundary] = new CheckBox(popup, "Boundary");
	mLayers[Layers::BoundaryWireframe] = new CheckBox(popup, "Boundary wireframe");
	mLayers[Layers::OtherEdge] = new CheckBox(popup, "Other positions");

	ctr = 0;
	for (auto l : mLayers) {
		l->setChecked(false);
		l->setId("layer." + std::to_string(ctr++));
	}
//morphing
	PopupButton *MorphingBtn = new PopupButton(window, "Morphing");
	Popup *morphPopup = MorphingBtn->popup();
	morphPopup->setAnchorHeight(61);

	morphPopup->setLayout(new GroupLayout());

	mEdgeTagging = new CheckBox(morphPopup, "Coloring");
	mEdgeTagging->setId("showedgetags");
	mEdgeTagging->setChecked(false);
	mEdgeTagging->setCallback([&](bool value) {
		if (!mRes.tetMesh())
			mRes.init_edge_tagging2D();
		else{
			mRes.init_edge_tagging3D();
		}
		auto const &R = mRes.E_rend;
		mExtractionResultShader.bind();
		mExtractionResultShader.uploadAttrib("position", MatrixXf(R.block(0, 0, 3, R.cols())));
		mExtractionResultShader.uploadAttrib("color", MatrixXf(R.block(3, 0, 3, R.cols())));

		mLayers[Layers::PositionField]->setChecked(false);
		mLayers[Layers::PositionSingularities]->setChecked(false);
		mLayers[Layers::Boundary]->setChecked(false);
	});
	Slider *slider2 = new Slider(morphPopup);
	slider2->setValue(0.0);
	auto cb = [&](Float value) {
		mRes.E_I_rend = (1 - value) * mRes.E_rend + value*mRes.E_O_rend;
	};
	cb(0.0f);
	slider2->setCallback(cb);
	slider2->setId("slider2");
//output
	PopupButton *ConstraintsBtn = new PopupButton(window, "Output");
	ConstraintsBtn->setIcon(ENTYPO_ICON_ROCKET);
	ConstraintsBtn->setBackgroundColor(Color(100, 0, 0, 25));
	Popup *ConstraintsPopup = ConstraintsBtn->popup();
	ConstraintsPopup->setAnchorHeight(61);

	ConstraintsPopup->setLayout(new GroupLayout());

	mShow_F_done = new CheckBox(ConstraintsPopup, "Face");
	mShow_F_done->setId("showdoneF");
	mShow_F_done->setChecked(false);
	mShow_E_done = new CheckBox(ConstraintsPopup, "Edge");
	mShow_E_done->setId("showdoneE");
	mShow_E_done->setChecked(false);

	mOutputBtn = new Button(ConstraintsPopup, "Output", ENTYPO_ICON_FLASH);
	mOutputBtn->setBackgroundColor(Color(0, 255, 0, 25));
	mOutputBtn->setCallback([&]() {
		char patho[300];
		if (!mRes.tetMesh()) {
			sprintf(patho, "%s%s", filename.c_str(), "_surout.obj");
			write_surface_mesh_OBJ(mRes.mV_tag, mRes.F_tag, patho);

#ifdef T_VTAG
			sprintf(patho, "%s%s", filename.c_str(), "_V_flag.txt");
			write_Vertex_Types_TXT(mRes.V_flag, patho);
#endif		
		}else {
			sprintf(patho, "%s%s", filename.c_str(), ".HYBRID");
			write_volume_mesh_HYBRID(mRes.mV_tag, mRes.F_tag, mRes.P_tag, mRes.Hex_flag, mRes.PF_flag, patho);
		}
	});
//slicing
	new Label(window, "Slicing plane", "sans-bold");
	Slider *slider = new Slider(window);
	mSplit = Vector4f(1.f, 0.f, 0.f, mRes.aabb().center().x());
	slider->setValue(0.0);
	auto cb2 = [&](Float value) {
		float offset = (mRes.aabb().max.x() - mRes.aabb().min.x()) * 0.5;
		mSplit.w() = -((1 - value) * (mRes.aabb().min.x() - offset) + value * (mRes.aabb().max.x() + offset));

		if (mShow_F_done->checked() || mShow_E_done->checked()) {
			if (!mRes.ECs.size()) return;
			if (!mRes.P_tag.size() && !mRes.F_tag.size()) return;
			//compute renderable faces
			std::vector<tuple_E> Es;
			std::vector<std::vector<uint32_t>> Fs, Fes;

			vector<bool> flag(mRes.ECs.size(), true), flag_F(mRes.F_tag.size(), false);
			for (uint32_t i = 0; i < mRes.ECs.size(); i++)
				if (mSplit.dot(mRes.ECs[i]) < 0) flag[i] = false;
			if (mRes.P_tag.size()) {
				for (uint32_t i = 0; i < flag.size(); i++)if (flag[i]) {
					for (auto f : mRes.P_tag[i]) flag_F[f] = !flag_F[f];
				}
			}
			else flag_F = flag;
			for (uint32_t i = 0; i < flag_F.size(); i++)if (flag_F[i]) {
				Fs.push_back(mRes.F_tag[i]);
				Fs.push_back(mRes.F_tag[i]); reverse(Fs[Fs.size() - 1].begin(), Fs[Fs.size() - 1].end());
			}
			if (!Fs.size()) return;

			mRes.construct_tEs_tFEs(Fs, Fes, Es);
			mRes.orient_polygon_mesh(mRes.mV_tag, Fs, Fes, Es);
			mRes.E_final_rend.setZero();
			mRes.E_final_rend.resize(6, 2 * Es.size());
			mRes.composit_edges_colors(mRes.mV_tag, Es, mRes.E_final_rend);
			mRes.composit_edges_centernodes_triangles(Fs, mRes.mV_tag, mRes.E_final_rend, mRes.mV_final_rend, mRes.F_final_rend);


			mExtractionResultShader_F_done.bind();
			mExtractionResultShader_F_done.uploadAttrib("position", mRes.mV_final_rend);
			mExtractionResultShader_F_done.uploadIndices(mRes.F_final_rend);

			auto const &R4 = mRes.E_final_rend;
			mExtractionResultShader_E_done.bind();
			mExtractionResultShader_E_done.uploadAttrib("position", MatrixXf(R4.block(0, 0, 3, R4.cols())));
			mExtractionResultShader_E_done.uploadAttrib("color", MatrixXf(R4.block(3, 0, 3, R4.cols())));
		}
	};
	cb2(0.0f);
	slider->setCallback(cb2);
	slider->setId("slider1");
//layout
    performLayout();
}

Viewer::~Viewer() {
    mOptimizer->shutdown();
    delete mOptimizer;
}

bool Viewer::mouseMotionEvent(const Vector2i &p, const Vector2i &rel,
                              int button, int modifiers) {
    if (!Screen::mouseMotionEvent(p, rel, button, modifiers)) {
        if (mCamera.arcball.motion(p)) {
            //repaint();
        } else if (mTranslate) {
            Eigen::Matrix4f model, view, proj;
            computeCameraMatrices(model, view, proj);
            float zval = project(mRes.aabb().center().cast<float>(), view * model, proj, mSize).z();
            Eigen::Vector3f pos1 =
                unproject(Eigen::Vector3f(p.x(), mSize.y() - p.y(), zval),
                          view * model, proj, mSize);
            Eigen::Vector3f pos0 = unproject(
                Eigen::Vector3f(mTranslateStart.x(),
                                mSize.y() - mTranslateStart.y(), zval),
                view * model, proj, mSize);
            mCamera.modelTranslation =
                mCamera.modelTranslation_start + (pos1 - pos0);
            //repaint();
        }
    }
    return true;
}

bool Viewer::mouseButtonEvent(const Vector2i &p, int button, bool down, int modifiers) {
    if (!Screen::mouseButtonEvent(p, button, down, modifiers)) {
        if (button == GLFW_MOUSE_BUTTON_1 && modifiers == 0) {
            mCamera.arcball.button(p, down);
        } else if (button == GLFW_MOUSE_BUTTON_2 ||
                   (button == GLFW_MOUSE_BUTTON_1 && modifiers == GLFW_MOD_SHIFT)) {
            mCamera.modelTranslation_start = mCamera.modelTranslation;
            mTranslate = true;
            mTranslateStart = p;
        }
    }
    if (button == GLFW_MOUSE_BUTTON_1 && !down)
        mCamera.arcball.button(p, false);
    if (!down) {
        mTranslate = false;
    }
    return true;
}

bool Viewer::resizeEvent(const Vector2i &size) {
    mCamera.arcball.setSize(mSize);
    return true;
}

bool Viewer::scrollEvent(const Vector2i &p, const Eigen::Vector2f &rel) {
    if (!Screen::scrollEvent(p, rel)) {
        mCamera.zoom = std::max(0.1, mCamera.zoom * (rel.y() > 0 ? 1.1 : 0.9));
        //repaint();
    }
    return true;
}

void Viewer::updatePositionSingularities() {
    if (mRes.tetMesh()) {
        mRes.detectPositionSingularitiesTet();
        auto const &S = mRes.positionSingularities();
        mPositionSingularityShaderTet.bind();
        mPositionSingularityShaderTet.uploadAttrib("position", MatrixXf(S.block(0, 0, 3, S.cols())));
        mPositionSingularityShaderTet.uploadAttrib("color",    MatrixXf(S.block(3, 0, 3, S.cols())));
    } else {
        mRes.detectPositionSingularitiesTri();
        auto const &S = mRes.positionSingularities();
        mPositionSingularityShaderTri.bind();
        mPositionSingularityShaderTri.uploadAttrib("position", MatrixXf(S.block(0, 0, 3, S.cols())));
        mPositionSingularityShaderTri.uploadAttrib("normal",   MatrixXf(S.block(3, 0, 3, S.cols())));
        mPositionSingularityShaderTri.uploadAttrib("color",    MatrixXf(S.block(6, 0, 3, S.cols())));
    }
}

void Viewer::updateOrientationSingularities() {
    if (mRes.tetMesh()) {
        mRes.detectOrientationSingularitiesTet();
        auto const &S = mRes.orientationSingularities();
        mOrientationSingularityShaderTet.bind();
        mOrientationSingularityShaderTet.uploadAttrib("position", MatrixXf(S.block(0, 0, 3, S.cols())));
        mOrientationSingularityShaderTet.uploadAttrib("color",    MatrixXf(S.block(3, 0, 3, S.cols())));
    } else {
        mRes.detectOrientationSingularitiesTri();
        auto const &S = mRes.orientationSingularities();
        mOrientationSingularityShaderTri.bind();
        mOrientationSingularityShaderTri.uploadAttrib("position", MatrixXf(S.block(0, 0, 3, S.cols())));
        mOrientationSingularityShaderTri.uploadAttrib("normal",   MatrixXf(S.block(3, 0, 3, S.cols())));
        mOrientationSingularityShaderTri.uploadAttrib("color",    MatrixXf(S.block(6, 0, 3, S.cols())));
    }
}

void Viewer::drawContents() {
	glClearColor(.5, .5, .5, 1.0f);
	//glClearColor(1, 1, 1, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

	mOptimizer->setAlignment(Configlayers[Config_Layers::Alignment]->checked());
	mOptimizer->setRandomization(Configlayers[Config_Layers::Randomization]->checked());
	mOptimizer->setHierarchy(Configlayers[Config_Layers::Hierarchy]->checked());

	mRes.setScale(mScaleBox->value());

	if (!mOptimizer->active()) {
		if (mSolveOrientationBtn->pushed()) {
			mSolveOrientationBtn->setPushed(false);
			updateOrientationSingularities();
		}
		if (mSolvePositionBtn->pushed()) {
			mSolvePositionBtn->setPushed(false);
			updatePositionSingularities();
		}
	}
	else if (!mOptimizer->hierarchy()) {
		if (mSolveOrientationBtn->pushed())
			updateOrientationSingularities();
	}
	if (mTmeshingBtn->pushed()) {
		mTmeshingBtn->setPushed(false);
	}
	if (mSolveDatastructureBtn->pushed()) {
		mSolveDatastructureBtn->setPushed(false);
	}
	Eigen::Matrix4f model, view, proj;
	computeCameraMatrices(model, view, proj);
	Eigen::Matrix4f mvp = proj * view * model;
	Eigen::Vector4f civ =
		(view * model).inverse() * Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);

	if (mRes.tetMesh()) {
		mOrientationFieldShaderTet.bind();
		mOrientationFieldShaderTet.uploadAttrib("q", mRes.Q());
	}
	else {
		mOrientationFieldShaderTri.bind();
		mOrientationFieldShaderTri.uploadAttrib("q", mRes.Q());
	}

	mPositionFieldShader.bind();
	mPositionFieldShader.uploadAttrib("o", mRes.O());
	//mPositionFieldShader.uploadAttrib("o", mRes.my_O());

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDisable(GL_BLEND);

	if (mLayers[Tetrahedra]->checked()) {
		mTetShader.bind();
		mTetShader.setUniform("light_position", mLightPosition);
		mTetShader.setUniform("model", model);
		mTetShader.setUniform("view", view);
		mTetShader.setUniform("proj", proj);
		mTetShader.setUniform("base_color", mBaseColor);
		//mTetShader.setUniform("specular_color", mSpecularColor);
		mTetShader.setUniform("split", mSplit);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		mTetShader.drawIndexed(GL_LINES_ADJACENCY, 0, mRes.tetCount() * 4);
		//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	glPointSize(5);
	// This must be enabled, otherwise glLineWidth has no effect
	//glEnable(GL_LINE_SMOOTH);
	//glLineWidth(5);
	if (mLayers[OrientationField]->checked()) {
		auto &shader = mRes.tetMesh() ? mOrientationFieldShaderTet : mOrientationFieldShaderTri;
		shader.bind();
		shader.setUniform("mvp", mvp);
		shader.setUniform("split", mSplit, false);
		shader.setUniform("scale", mRes.averageEdgeLength() / 3);
		shader.drawArray(GL_POINTS, 0, mRes.tetMesh() ? mRes.tetCount() : mRes.faceCount());
	}

	if (mLayers[PositionField]->checked()) {
		mPositionFieldShader.bind();
		mPositionFieldShader.setUniform("mvp", mvp);
		mPositionFieldShader.setUniform("split", mSplit, false);
		mPositionFieldShader.drawArray(GL_POINTS, 0, mRes.vertexCount());
	}
	//if (mLayers[OtherEdge]->checked()) {
	//	mOtherEdge.bind();
	//	mOtherEdge.setUniform("mvp", mvp);
	//	mOtherEdge.setUniform("split", mSplit, false);
	//	mOtherEdge.drawArray(GL_POINTS, 0, mRes.vertexCount());
	//}
	if (mLayers[OrientationSingularities]->checked()) {
		auto &shader = mRes.tetMesh() ? mOrientationSingularityShaderTet : mOrientationSingularityShaderTri;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.setUniform("scale", mRes.averageEdgeLength(), false);
		shader.drawArray(mRes.tetMesh() ? GL_LINES : GL_POINTS, 0, mRes.orientationSingularities().cols());
	}

	if (mLayers[PositionSingularities]->checked()) {
		auto &shader = mRes.tetMesh() ? mPositionSingularityShaderTet : mPositionSingularityShaderTri;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.setUniform("scale", mRes.averageEdgeLength(), false);
		shader.drawArray(mRes.tetMesh() ? GL_LINES : GL_POINTS, 0, mRes.positionSingularities().cols());
	}

	if (mLayers[Boundary]->checked()) {
		if (mRes.tetMesh()) {
			glEnable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}
		mMeshShader.bind();
		mMeshShader.setUniform("light_position", mLightPosition);
		mMeshShader.setUniform("model", model);
		mMeshShader.setUniform("view", view);
		mMeshShader.setUniform("proj", proj);
		mMeshShader.setUniform("base_color", mBaseColorBoundary);
		//mMeshShader.setUniform("specular_color", mSpecularColorBoundary);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		if (mRes.tetMesh()) {
			glDepthFunc(GL_LEQUAL);
			glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
			mMeshShader.drawIndexed(GL_TRIANGLES, 0, mRes.faceCount());
			glDepthFunc(GL_EQUAL);
			glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
			mMeshShader.drawIndexed(GL_TRIANGLES, 0, mRes.faceCount());
			glDepthFunc(GL_LEQUAL);
		}
		else {
			mMeshShader.drawIndexed(GL_TRIANGLES, 0, mRes.faceCount());
		}
		glDisable(GL_POLYGON_OFFSET_FILL);
	}

	if (mLayers[BoundaryWireframe]->checked()) {
		mMeshShader.bind();
		mMeshShader.setUniform("light_position", mLightPosition);
		mMeshShader.setUniform("model", model);
		mMeshShader.setUniform("view", view);
		mMeshShader.setUniform("proj", proj);
		mMeshShader.setUniform("base_color", Vector4f(Vector4f::Constant(0.f)));
		mMeshShader.setUniform("specular_color", Vector4f(Vector4f::Constant(0.f)));
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		mMeshShader.drawIndexed(GL_TRIANGLES, 0, mRes.faceCount());
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	if (mEdgeTagging->checked()) {
		mExtractionResultShader.bind();
		mExtractionResultShader.uploadAttrib("position", MatrixXf(mRes.E_I_rend.block(0, 0, 3, mRes.E_I_rend.cols())));
		mExtractionResultShader.uploadAttrib("color", MatrixXf(mRes.E_I_rend.block(3, 0, 3, mRes.E_I_rend.cols())));
		auto &shader = mExtractionResultShader;
		shader.bind();
		shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.drawArray(GL_LINES, 0, mRes.E_I_rend.cols());
	}


	if (mShow_F_done->checked()) {
		mExtractionResultShader_F_done.bind();
		mExtractionResultShader_F_done.setUniform("light_position", mLightPosition);
		mExtractionResultShader_F_done.setUniform("model", model);
		mExtractionResultShader_F_done.setUniform("view", view);
		mExtractionResultShader_F_done.setUniform("proj", proj);
		//mExtractionResultShader_F_done.setUniform("split", mSplit, false);
		mExtractionResultShader_F_done.setUniform("base_color", mBaseColorBoundary);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		mExtractionResultShader_F_done.drawIndexed(GL_TRIANGLES, 0, mRes.F_final_rend.cols());
		glDisable(GL_POLYGON_OFFSET_FILL);
	}
	if (mShow_E_done->checked())
	{
		auto &shader = mExtractionResultShader_E_done;
		shader.bind();
		//shader.setUniform("split", mSplit, false);
		shader.setUniform("mvp", mvp);
		shader.drawArray(GL_LINES, 0, mRes.E_final_rend.cols());
	}

}

bool Viewer::keyboardEvent(int key, int scancode, int action, int modifiers) {
    if (Screen::keyboardEvent(key, scancode, action, modifiers))
        return true;
    if (action != GLFW_PRESS)
        return false;
     return false;
}

void Viewer::computeCameraMatrices(Eigen::Matrix4f &model,
                                   Eigen::Matrix4f &view,
                                   Eigen::Matrix4f &proj) {
    view = lookAt(mCamera.eye, mCamera.center, mCamera.up);

    float fH = std::tan(mCamera.viewAngle / 360.0f * M_PI) * mCamera.dnear;
    float fW = fH * (float) mSize.x() / (float) mSize.y();

    proj = frustum(-fW, fW, -fH, fH, mCamera.dnear, mCamera.dfar);
    model = mCamera.arcball.matrix();

	model = model * scale(Eigen::Vector3f::Constant(mCamera.zoom * mCamera.modelZoom));
	model = model * translate(mCamera.modelTranslation);
}


