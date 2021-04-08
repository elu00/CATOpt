#include "CatOpt.h"
CatOpt::CatOpt(string s) {
    inputMeshPath = s;
    cout << "Initialized" << endl;

    //polyscope::init("openGL_mock");
    polyscope::init();

    //polyscope::state::userCallback = myCallback;
    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        "main",
        geometry->inputVertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));
    cout << "starting optimization" << endl;
    initializeQuantities();
    generateConstraints();

    subdivision();
    buildNewMesh();

    //fin = descent();

    //////////////////////////////////////////
    // BFF stuff
    conformalFlatten();
    //testSVG();

    //generateVisualization();
    // Give control to the polyscope gui
    //polyscope::show();

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
        //inputMeshPath = "/home/elu/repos/catopt/meshes/spotwithhole.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/cube.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    CatOpt c (inputMeshPath);
    
    return EXIT_SUCCESS;
}