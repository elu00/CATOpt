#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "optimo/solvers/sparse_primal_dual_qp.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;
using std::unique_ptr;

using optimo::solvers::SparsePrimalDualQP;

using Eigen::Matrix;
using Eigen::Dynamic;

// == Geometry-central data
unique_ptr<HalfedgeMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

// Mesh data
size_t nVertices;
size_t nEdges;
size_t nCorners;
EdgeData<size_t> eInd;
VertexData<size_t> vInd;
CornerData<size_t> cInd;
CornerData<double> cornerAngles;
VertexData<double> angleDefects; 

// Optimization Stuff
SparseMatrix<double> constraints;
Vector<double> x_init;
Vector<double> x;
Vector<double> rhs;
Vector<double> ineqRHS;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

void generateConstraints()
{
    // Initialization
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    cInd = mesh->getCornerIndices();
    geometry->requireVertexGaussianCurvatures();
    angleDefects = geometry->vertexGaussianCurvatures;

    // Equality constraint initialization
    constraints = SparseMatrix<double>(nVertices, nEdges);
    std::vector<Eigen::Triplet<double>> tripletList;
    rhs = Vector<double>(nVertices);
    for (size_t i = 0; i < nVertices; i++)
    {
        rhs[i] = angleDefects[mesh->vertex(i)] / 2.;
    }
    for (Edge e : mesh->edges())
    {
        tripletList.emplace_back(vInd[e.halfedge().vertex()], eInd[e], 1);
        tripletList.emplace_back(vInd[e.halfedge().twin().vertex()], eInd[e], 1);
    }
    constraints.setFromTriplets(tripletList.begin(), tripletList.end());

    return;
}

bool checkInequalityConstraints()
{
    cornerAngles = CornerData<double>(*mesh);
    CornerData<double> netAngles = CornerData<double>(*mesh);
    for (Corner c : mesh->corners())
    {
        cornerAngles[c] = geometry->cornerAngle(c);
        netAngles[c] = cornerAngles[c];
    }
    for (Corner c : mesh->corners())
    {
        Halfedge h = c.halfedge();
        netAngles[c] += x_init[eInd[h.edge()]];
        netAngles[h.next().corner()] += x_init[eInd[h.edge()]];
    }
    for (Corner c : mesh->corners())
    {
        if (netAngles[c] < 0 || netAngles[c] > 2 * PI)
        {
            return false;
        }
    }
    return true;
}

void generateVisualization()
{
    // Visualization
    /*
    EdgeData<double> initialGuess(*mesh);
    EdgeData<double> finalSolution(*mesh);
    for (Edge e : mesh->edges())
    {
        initialGuess[e] = x_init[eInd[e]];
        finalSolution[e] = x[eInd[e]];
        //cout << x[eInd[e]] << endl;
        //cout << x_init[eInd[e]] << endl;
    }
    */
    psMesh->addEdgeScalarQuantity("Initial Guess", x_init);
    psMesh->addEdgeScalarQuantity("Final Solution", x);
    psMesh->addVertexScalarQuantity("curvature",
            geometry->vertexGaussianCurvatures);
}

void optimizationStep()
{ 
    // Initialization
    // parameters are number of unknowns, number of inequality contraints, number of eq constraints
    const uint n = nEdges;
    const uint m = 2*nCorners;
    const uint p = nVertices;
    const uint l = n + m + p;
    SparsePrimalDualQP<double>::Params params(n, m, p);
    double min_value;
    params.reserve();
    params.d = Vector<double>::Zero(nEdges,1);
    // Equality constraints
    params.beq = rhs;
    for (Edge e : mesh->edges())
    {
        params.setAeqElement(vInd[e.halfedge().vertex()], eInd[e], 1);
        params.setAeqElement(vInd[e.halfedge().twin().vertex()], eInd[e], 1);
    }

    // Inequality constraints 
    Matrix<double, Dynamic, 1>& bin = params.bin;   
    for (Corner c : mesh->corners())
    {
        bin[cInd[c]] = 2*PI - geometry->cornerAngle(c);
        bin[cInd[c] + nCorners] = -geometry->cornerAngle(c);
        Halfedge h = c.halfedge();
        params.setAinElement(cInd[c], eInd[h.edge()] , 1);
        params.setAinElement(cInd[c] + nCorners, eInd[h.edge()] , -1);
        params.setAinElement(cInd[h.next().corner()], eInd[h.edge()] , 1);
        params.setAinElement(cInd[h.next().corner()] + nCorners, eInd[h.edge()] , -1);

    }

    // Optimization variable lol
    for (uint i = 0; i < n; i++) 
    {     
        params.setQElement(i, i, 1);
    }

    // Plug in initial guess, I think?
    Matrix<double, Dynamic, 1> y(l);
    for (uint i = 0; i < n; i++) y(0, i) = x_init[i];
    SparsePrimalDualQP<double> qp_solver;   
    auto res = qp_solver(&params, &y, &min_value);
    x = y.block(0, 0, n, 1);
    cout << "RESULT: " << res << endl;
    cout << "VALUE: " << min_value << endl;
    cout << "Optimization Done" << endl;
}



int main(int argc, char **argv)
{

    // Configure the argument parser
    args::ArgumentParser parser("geometry-central & Polyscope example project");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    // Make sure a mesh name was given
    if (!inputFilename)
    {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }
    // Initialize polyscope
    polyscope::init();

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));
    generateConstraints();
    cout << "Constraints generated" << endl;

    // Initialize with simple linear solve
    // This will screw up if SuiteSparse doesn't exist
    x_init = solve(constraints, rhs);
    if (checkInequalityConstraints())
    {
        cout << "No optimization required!" << endl;
        x = x_init;
    }
    else
    {
        cout << "Doing optimization now..." << endl;
        optimizationStep();
    }

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
            polyscope::guessNiceNameFromPath(args::get(inputFilename)),
            geometry->inputVertexPositions, mesh->getFaceVertexList(),
            polyscopePermutations(*mesh));

    generateVisualization();

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
