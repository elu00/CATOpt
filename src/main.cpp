#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "fusion.h"

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
Vector<double> x_init;
Vector<double> x;
vector<double> rhs;
vector<double> ineqRHS;

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
	auto Meq = Matrix::sparse(r->size(), c->size(), r, c, v); 
    auto eqRHS = new_array_ptr(rhs);

    M->constraint("eq constraints", Expr::mul(Meq, x), Domain::equalsTo(eqRHS));
    cout << "eq generated" << endl;
    // inequality constraints
    ineqRHS = vector<double>(nCorners);
    rows.clear();
    cols.clear();
    values.clear();
    for (Corner c : mesh->corners())
    {
        ineqRHS[cInd[c]] = 2*PI - geometry->cornerAngle(c);
        ineqRHS[cInd[c] + nCorners] = -geometry->cornerAngle(c);
        Halfedge h = c.halfedge();
        rows.emplace_back(cInd[c]);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(1.);
        rows.emplace_back(cInd[c] + nCorners);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(-1.);
        rows.emplace_back(cInd[h.next().corner()]);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(1.);
        rows.emplace_back(cInd[h.next().corner()] + nCorners);
		cols.emplace_back(eInd[h.edge()]);
		values.emplace_back(-1.);

    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
	auto Mineq = Matrix::sparse(r->size(), c->size(), r, c, v); 
    auto inRHS = new_array_ptr(ineqRHS);
    M->constraint("ineq constraints", Expr::mul(Mineq, x), Domain::lessThan(inRHS));

    cout << "ineq generated" << endl;
    M->objective(ObjectiveSense::Minimize, x);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << "RESULT: " << "hehe" << endl;
    cout << "VALUE: " << "xd" << endl;
    cout << "Optimization Done" << endl;

    return;
}

void generateVisualization()
{

}
void optimizationStep()
{
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
        // Do optimization
        cout << "starting optimization";
        generateConstraints();

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
