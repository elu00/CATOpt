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
#include "CircleWrapper.h"

shared_ptr<ManifoldSurfaceMesh> mesh;
shared_ptr<VertexPositionGeometry> geometry;

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
        //inputMeshPath = "/home/elu/repos/catopt/meshes/cube.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/spotwithhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/plane.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/beanhole.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/test.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/patch.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/BumpyTorusPatch.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    // polyscope sanity checks
    polyscope::init();
    polyscope::SurfaceMesh *psMesh = polyscope::registerSurfaceMesh(
      "waaah",
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

    mesh->compress();
    IntrinsicFlattening flattener(mesh, geometry);
/*
    EdgeData<double> intersectionAngles = flattener.solveKSS();
    CircleWrapper patterns(mesh, intersectionAngles, psMesh);
    patterns.solveKSS();
    return EXIT_SUCCESS;
    */
    SolutionData sol = flattener.solveFromPlane();
    VertexData<double> vertexSum(*mesh);
    for (Vertex v: mesh->vertices()) {
        double accum = 0;
        for (Corner c: v.adjacentCorners()) {
            accum += sol.betas[c];
        }
        vertexSum[v] = accum;
    }
    psMesh->addVertexScalarQuantity("vertex sums", vertexSum);
    //psMesh->addHalfedgeScalarQuantity("solved", sol.betas);
    psMesh->addEdgeScalarQuantity("thetas", sol.thetas);
    //psMesh->addVertexScalarQuantity("infinite vertex", sol.infVertex);
    psMesh->addFaceScalarQuantity("faces", sol.fMask);

    CircleWrapper patterns(mesh, sol, psMesh);
    patterns.solve();
    


    polyscope::show();
    
    return EXIT_SUCCESS;
}