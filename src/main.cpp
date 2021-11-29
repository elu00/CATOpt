#include "args/args.hxx"
#include <string>
#include <iostream>
#include <memory>
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
        //inputMeshPath = "/home/elu/repos/catopt/meshes/plane.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/beanhole.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/test.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/patch.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/BumpyTorusPatch.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    mesh->compress();
    IntrinsicFlattening flattener(mesh, geometry);
    SolutionData sol = flattener.solve();
    CircleWrapper patterns(mesh, sol);
    patterns.solve();
    
    return EXIT_SUCCESS;
}