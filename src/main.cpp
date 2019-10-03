#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "fusion.h"
#include <algorithm>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;
using std::unique_ptr;

using namespace mosek::fusion;
using namespace monty;
using std::vector;

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
//SparseMatrix<double> constraints;
//Vector<double> x_init;
vector<double> *sol;
vector<double> rhs;
vector<double> ineqRHS0;
vector<double> ineqRHS1;

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

    // Model initialization
    Model::t M = new Model(); auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", nEdges, Domain::inRange(-2*PI, 2 * PI));

   	vector<int> rows;
	vector<int> cols;
	vector<double> values;
    // Equality constraint initialization
    rhs = vector<double>(nVertices);
    for (size_t i = 0; i < nVertices; i++)
    {
        rhs[i] = angleDefects[mesh->vertex(i)] / 2.;
    }
    for (Edge e : mesh->edges())
    {
        rows.emplace_back(vInd[e.halfedge().vertex()]);
		cols.emplace_back(eInd[e]);
		values.emplace_back(1.);
        rows.emplace_back(vInd[e.halfedge().twin().vertex()]);
		cols.emplace_back(eInd[e]);
		values.emplace_back(1.);
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
	auto Meq = Matrix::sparse(nVertices, nEdges, r, c, v); 
    auto eqRHS = new_array_ptr(rhs);
    M->constraint("eq constraints", Expr::mul(Meq, x), Domain::equalsTo(eqRHS));
    cout << "eq generated" << endl;
    // inequality constraints
    ineqRHS0 = vector<double>(nCorners);
    ineqRHS1 = vector<double>(nCorners);
    rows.clear();
    cols.clear();
    values.clear();
    for (Corner c : mesh->corners())
    {
        ineqRHS0[cInd[c]] = - geometry->cornerAngle(c);
        ineqRHS1[cInd[c]] = 2*PI - geometry->cornerAngle(c);
        Halfedge h = c.halfedge();
        rows.emplace_back(cInd[c]);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(1.);
        rows.emplace_back(cInd[h.next().corner()]);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(1.);
    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
	auto Mineq = Matrix::sparse(nCorners, nEdges, r, c, v); 
    auto inRHS0 = new_array_ptr(ineqRHS0);
    auto inRHS1 = new_array_ptr(ineqRHS1);
    // sum is greater than 0
    M->constraint("ineq0 constraints", Expr::mul(Mineq, x), Domain::greaterThan(inRHS0));
    // sum is less thatn 2pi
    M->constraint("ineq1 constraints", Expr::mul(Mineq, x), Domain::lessThan(inRHS1));

    auto ones =  std::make_shared<ndarray<double,1>>(shape(nEdges),1.);
    cout << "ineq generated" << endl;
    M->objective(ObjectiveSense::Minimize, Expr::dot(ones,x));
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << x->level() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = x->getSize();
    auto xVal = x->level();
    std::cout << "\nOptimal primal objective: " << M->primalObjValue() <<"\n";
    sol = new vector<double>(nEdges);
    for(int i = 0; i < xsize; ++i)
    {
       (*sol)[i] =  (*xVal)[i];
    }
    return;
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
    //psMesh->addEdgeScalarQuantity("Initial Guess", x_init);
    //psMesh->addEdgeScalarQuantity("Final Solution", *sol);
    //psMesh->addVertexScalarQuantity("curvature",
    //        geometry->vertexGaussianCurvatures);
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
        //polyscope::init();

        // Load mesh
        std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));
        // Do optimization
        cout << "starting optimization";
        generateConstraints();

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
            polyscope::guessNiceNameFromPath(args::get(inputFilename)),
            geometry->inputVertexPositions, mesh->getFaceVertexList(),
            polyscopePermutations(*mesh));

    generateVisualization();

    cout << "registering mesh" << endl; 
    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
