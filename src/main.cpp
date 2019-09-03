#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;
using std::unique_ptr;

// == Geometry-central data
unique_ptr<HalfedgeMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

// Mesh data
size_t nVertices;
size_t nEdges;
EdgeData<size_t> eInd;
VertexData<size_t> vInd;
CornerData<double> cornerAngles;

// Optimization Parameters
double t = 2.;
double mew = 1.5;
double epsilon = 0.001;

// Optimization Stuff
SparseMatrix<double> constraints;
Vector<double> x;
Vector<double> rhs;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

void generateConstraints()
{
    // Initialization
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    geometry->requireVertexGaussianCurvatures();
    VertexData<double> angleDefects = geometry->vertexGaussianCurvatures;

    // Matrix initialization
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
        netAngles[c] += x[eInd[h.edge()]];
        netAngles[h.next().corner()] += x[eInd[h.edge()]];
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

void barrierOpt()
{
}

inline double logBarrier(const double xi, const double upperBound)
{
    return -(1 / t) * log(-(upperBound - xi));
}

inline void vecProj(SparseMatrix<double> A, Vector<double> xi, Vector<double> b)
{
}
double objFunction()
{
    return 0.;
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
    // Initialize with simple linear solve
    // This will screw up if SuiteSparse doesn't exist
    x = solve(constraints, rhs);
    if (checkInequalityConstraints())
    {
        cout << "No log-barrier required!" << endl;
    }
    else
    {
        cout << "Doing log-barrier now..." << endl;
        barrierOpt();
    }

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(args::get(inputFilename)),
        geometry->inputVertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));

    // Visualization
    psMesh->addEdgeScalarQuantity("Alphas", x);
    psMesh->addVertexScalarQuantity("curvature",
                                  geometry->vertexGaussianCurvatures);

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
