// parsing and STL stuff
#include "args/args.hxx"
#include <iostream>
#include <memory>
#include <string>

// polyscope dependencies
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
using std::shared_ptr;
using std::string;

// geometry-central includes
#include "Common.h"

#include "EmbeddingOptimization.h"
#include "IntrinsicFlattening.h"
#include "PrescribeCurvature.h"

// geometry-central data
shared_ptr<ManifoldSurfaceMesh> mesh;
shared_ptr<VertexPositionGeometry> geometry;
// polyscope handle
polyscope::SurfaceMesh *psMesh;

// test to see if optimization stuff works
void planarMapping() {
    geometry->requireEdgeLengths();

    // calculate desired betas - in this case, we're just using the flat angles
    geometry->requireCornerAngles();
    CornerData<double> beta = geometry->cornerAngles;

    // calculate desired vertex curvatures
    VertexData<double> VertexCurvatures(*mesh);
    for (Vertex v : mesh->vertices()) {
        double accum = v.isBoundary() ? M_PI : 2 * M_PI;
        for (Corner c : v.adjacentCorners()) {
            accum -= beta[c];
        }
    }

    // calculate edge curvatures
    EdgeData<double> EdgeCurvatures(*mesh);
    for (Edge e : mesh->edges()) {
        EdgeCurvatures[e] = 0;
    }
    auto [l, betaHat] = PrescribeCurvature(mesh, geometry->edgeLengths, beta,
                                           VertexCurvatures, EdgeCurvatures);
    // Seam edges
    EdgeData<bool> S(*mesh, false);
    CornerData<Vector2> positions = LayoutMesh(mesh, l, S);
    CATToSVG(mesh, positions, betaHat, "test.svg");
}

EmbeddingOptimization *E;
int subdivisions = 3;

// intrinsically flatten the mesh, then perform an NxN subdivision on each
// triangle.
void embeddingInitialization(int N) {
    geometry->requireEdgeLengths();
    IntrinsicFlattening flattener(mesh, geometry->edgeLengths);
    cout << "Solving for intrinstic angles..." << endl;
    // corner angles from our solve
    CornerData<double> beta = flattener.SolveIntrinsicOnly();
    cout << "Solved, checking validity of solution..." << endl;
    flattener.CheckConstraintsIntrinsicOnly(beta);
    cout << "Solution is valid" << endl;
    E = (new EmbeddingOptimization(mesh, geometry, beta));
    auto [submesh, subgeometry] = E->initializeSubdivision(N);
    cout << "EmbeddingOptimization initialized" << endl;
    E->initializeLM();
    cout << "LM initialized" << endl;
}
// TODO: write a function to map a surface to the plane
void surfaceToPlane() {
    geometry->requireEdgeLengths();
    IntrinsicFlattening flattener(mesh, geometry->edgeLengths);
    /*
  EdgeData<double> intersectionAngles = flattener.solveKSS();
  CircleWrapper patterns(mesh, intersectionAngles, psMesh);
  //patterns.solveKSS();
  */
}
// Imgui parameters
bool LMInitialized = false;
int MAX_ITERS = 100;
double fairnessNormalization = 1;
void imguiCallback() {
    ImGui::SliderInt("Number of subdivisions", &subdivisions, 2, 6);
    if (ImGui::Button("Create subdivisions")) {
        embeddingInitialization(subdivisions);
        LMInitialized = true;
    }
    // run LM optimization on subdivided mesh
    if (LMInitialized) {
        ImGui::InputInt("Max iteration count", &MAX_ITERS, 1, 2000);
        ImGui::InputDouble("Fairness Weight", &fairnessNormalization, 1e-10,
                           100.0);
        if (ImGui::Button("Run LM")) {
            E->fairnessNormalization = fairnessNormalization;
            E->optimizeOneStep(MAX_ITERS);
        }
    }
}

int main(int argc, char **argv) {
    //===========================Initialization=================================
    // Configure the argument parser
    string inputMeshPath;
    args::ArgumentParser parser("mesh file name");
    args::Positional<std::string> inputFilename(parser, "mesh", "path to mesh");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // choose a default mesh if none is provided
    if (!inputFilename) {
        inputMeshPath = "/home/elu/repos/catopt/meshes/beanhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/tetrahedron.obj";
        //  inputMeshPath = "/home/elu/repos/catopt/meshes/plane.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    mesh->compress();
    //===========================Polyscope setup================================
    polyscope::init();
    // initialize view of original mesh
    /*
        psMesh = polyscope::registerSurfaceMesh(
      "original geometry",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));
  */

    // planarMapping();
    polyscope::state::userCallback = imguiCallback;
    polyscope::show();

    return EXIT_SUCCESS;
}
