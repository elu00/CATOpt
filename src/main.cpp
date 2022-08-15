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
#include "CircleWrapper.h"


shared_ptr<ManifoldSurfaceMesh> mesh;
shared_ptr<VertexPositionGeometry> geometry;
polyscope::SurfaceMesh *psMesh;

void planarMapping(int N) {
    for (int m = 0; m <= N; m++){
        IntrinsicFlattening flattener(mesh, geometry);
        auto [thetas, betas] = flattener.solveFromPlane((double)m/(double)N);

        CircleWrapper patterns(mesh, betas, psMesh);
        patterns.solve("fin" + std::string(3 - std::to_string(m).length(), '0') +  std::to_string(m) );
    }
}

EmbeddingOptimization* E;
float t = 1e-4;
void embedding(int N) {
    IntrinsicFlattening flattener(mesh, geometry);
    cout << "solving intrinsic..." << endl;
    CornerData<double> beta = flattener.solveIntrinsicOnly();
    cout << "solved" << endl;
    E = (new EmbeddingOptimization(mesh, geometry, beta));
    auto [submesh, subgeometry] = E->solve(N);

}
void surfaceToPlane() {
    IntrinsicFlattening flattener(mesh, geometry);
    /*
    EdgeData<double> intersectionAngles = flattener.solveKSS();
    CircleWrapper patterns(mesh, intersectionAngles, psMesh);
    //patterns.solveKSS();
    */
}
void myCallback() {
    ImGui::InputFloat("param value", &t, 0.01f, 1.0f);  // set a float variable
    if (ImGui::Button("run subroutine")) {
        E->optimize(t);
    }
}

int main(int argc, char **argv) {
    // Configure the argument parser
    string inputMeshPath;
    args::ArgumentParser parser("");
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
        inputMeshPath = "/home/elu/repos/catopt/meshes/spotwithhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/SmallDisk.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/BumpyTorusPatch.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/triangle.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/beanhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/square.obj";
        inputMeshPath = "../meshes/plane.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/test.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/patch.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    // polyscope sanity checks
    polyscope::init();
    /*
       psMesh = polyscope::registerSurfaceMesh(
       "waaah",
       geometry->inputVertexPositions, mesh->getFaceVertexList(),
       polyscopePermutations(*mesh));
       */

    mesh->compress();
    polyscope::state::userCallback = myCallback;
    //embedding(2);
    planarMapping(5);
    //planarMapping(100);

    //polyscope::show();

    return EXIT_SUCCESS;
}
