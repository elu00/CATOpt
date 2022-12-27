#include "args/args.hxx"
#include <string>
#include <iostream>
#include <memory>


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
using std::string;
using std::shared_ptr;

#include "Common.h"


#include "IntrinsicFlattening.h"
#include "EmbeddingOptimization.h"


shared_ptr<ManifoldSurfaceMesh> mesh;
shared_ptr<VertexPositionGeometry> geometry;
polyscope::SurfaceMesh *psMesh;

// TODO: rewrite this
void planarMapping(int N) {
    for (int m = 0; m <= N; m++){
        geometry->requireEdgeLengths();
        IntrinsicFlattening flattener(mesh, geometry->edgeLengths);
        auto [thetas, betas] = flattener.solveFromPlane((double)m/(double)N);

        /*
        CircleWrapper patterns(mesh, betas, psMesh);
        patterns.solve("fin" + std::string(3 - std::to_string(m).length(), '0') +  std::to_string(m) );
        */
    }
}

EmbeddingOptimization* E;
int subdivisions = 3;

// intrinsically flatten the mesh, then try to embed it in the plane with a N * N subdivision on each triangle.
void embedding(int N) {
    geometry->requireEdgeLengths();
    IntrinsicFlattening flattener(mesh, geometry->edgeLengths);
    cout << "solving intrinsic..." << endl;
    CornerData<double> beta = flattener.solveIntrinsicOnly();
    cout << "solved" << endl;
    E = (new EmbeddingOptimization(mesh, geometry, beta));
    auto [submesh, subgeometry] = E->initializeSubdivision(N);
    cout << "EmbeddingOptimization initialized" << endl;
    E->initializeLM();
    cout << "LM initialized" << endl;
}
// TODO: rewrite this
void surfaceToPlane() {
    geometry->requireEdgeLengths();
    IntrinsicFlattening flattener(mesh, geometry->edgeLengths);
    /*
    EdgeData<double> intersectionAngles = flattener.solveKSS();
    CircleWrapper patterns(mesh, intersectionAngles, psMesh);
    //patterns.solveKSS();
    */
}
bool LMInitialized = false;
int MAX_ITERS = 100;
double fairnessNormalization = 1.;
void myCallback() {
    ImGui::SliderInt("Number of subdivisions", &subdivisions, 2, 5);  // set a float variable
    if (ImGui::Button("Create subdivisions")) {
        embedding(subdivisions);
        LMInitialized = true;
    }
    if (LMInitialized) {
        ImGui::InputInt("Max iteration count", &MAX_ITERS, 1, 2000);  // set a float variable
        ImGui::InputDouble("Fairness Weight", &fairnessNormalization, 0.01, 100.0);
        if (ImGui::Button("Run LM")) {
            //cout << E->fairnessNormalization << endl;
            //E->fairnessNormalization = fairnessNormalization;
            E->optimizeOneStep(MAX_ITERS);
        }
    }
}

int main(int argc, char **argv) {
    // Configure the argument parser
    string inputMeshPath;
    args::ArgumentParser parser("mesh file name");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a mesh name was given
    if (!inputFilename) {
        inputMeshPath = "../meshes/beanhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/tetrahedron.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    // polyscope sanity checks
    polyscope::init();
    psMesh = polyscope::registerSurfaceMesh(
            "original geometry",
            geometry->inputVertexPositions, mesh->getFaceVertexList(),
            polyscopePermutations(*mesh));

    mesh->compress();
    polyscope::state::userCallback = myCallback;

    polyscope::show();

    return EXIT_SUCCESS;
}
