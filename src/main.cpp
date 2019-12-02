#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

#include "args/args.hxx"
#include "imgui.h"

#include "fusion.h"
#include <algorithm>
#include <tuple>

#include "svg_gen.hpp"

#include "glm/vec3.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;
using std::unique_ptr;

using namespace mosek::fusion;
using namespace monty;
using std::vector;
using std::tuple;
using std::make_tuple;

// == Geometry-central data
unique_ptr<HalfedgeMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

// Mesh data
size_t nVertices;
size_t nFaces;
size_t nEdges;
size_t nCorners;
EdgeData<size_t> eInd;
VertexData<size_t> vInd;
CornerData<size_t> cInd;
FaceData<size_t> fInd;
CornerData<double> cornerAngles;
VertexData<double> angleDefects; 

// Optimization Stuff
// SparseMatrix<double> constraints;
// Vector<double> x_init;
vector<double>* sol;
vector<double> rhs;
vector<double> ineqRHS0;
vector<double> ineqRHS1;

size_t subdiv_level = 8;
size_t int_count;
size_t edge_count;
size_t vertex_count;
vector<Vector3> subdiv_points;
vector<tuple<size_t, size_t, double>> correct_dist;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

void generateConstraints()
{
    // Initialization
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    nFaces = mesh->nFaces();
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    cInd = mesh->getCornerIndices();
    fInd = mesh->getFaceIndices();
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
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint( Expr::vstack(t, x), Domain::inQCone() );
    M->objective( ObjectiveSense::Minimize, t );
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

  
void generateSVGs()
{
    int n = 0;
    for (Face f: mesh->faces())
    {
        double ij, jk, ki, a_ij, a_jk, a_ki;
        Halfedge h = f.halfedge();
        ij = geometry->edgeLength(h.edge()) * 100;
        a_ij = (*sol)[eInd[h.edge()]];
        h = h.next();
        jk = geometry->edgeLength(h.edge()) * 100;
        a_jk = (*sol)[eInd[h.edge()]];
        h = h.next();
        ki = geometry->edgeLength(h.edge()) * 100;
        a_ki = (*sol)[eInd[h.edge()]];
        CAT c = CAT(ij, jk, ki, a_ij, a_jk, a_ki);
        cout << ij << endl;
        cout << c.a_ij << endl;
        c.to_svg("current" + std::to_string(n) + ".svg");
        n++;
    }
}
inline Vector3 bary (Face f, double a, double b, double c)
{
    //cout << a + b + c << endl;
    auto it = f.halfedge();
    Vector3 i = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 j = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 k = geometry->inputVertexPositions[it.vertex()];
    return a*i + b*j + c*k;
}
inline double l2_dist(Face f, int j1, int k1, int j2, int k2)
{
    return 5.;
}


inline void interior_neighbors(Face f, int i, int j)
{

}
/* INDEXING CONVENTIONS:
 * Interior points - edge points - vertices
 * Interior points and edge points always have index in (0, subdiv_level) (note the strict inequality)
 */
inline size_t get_interior_index(Face f, int i, int j)
{
    return 6;
}
inline size_t get_edge_index(Face e, int i)
{
    return 7;
}
inline size_t get_vertex_index(Vertex v)
{
    return 5;
}
void subdivision()
{
    int_count = nFaces * (subdiv_level-1) * (subdiv_level-2)/2;
    edge_count = nEdges * (subdiv_level-1);
    vertex_count = nVertices;
    subdiv_points = vector<Vector3> (int_count + edge_count + vertex_count);
    // initialization of initial guesses
    for (Face f: mesh->faces())
    {
        size_t index = 0;
        for (size_t j = 1; j < subdiv_level - 1; j++)
        {
            for (size_t k = 1; k < subdiv_level - j - 1; k++)
            {
                Vector3 coordinate = bary(f, 1.-((double)(j+k))/subdiv_level, ((double)j)/subdiv_level, ((double)k)/subdiv_level);
                subdiv_points[fInd[f] * (subdiv_level-1) * (subdiv_level-2)/2 + index] = coordinate;
                index++;
                //cout << coordinate[0] << "," << coordinate[1] << "," << coordinate[2] << endl;
                //subdiv_points.push_back(glm::vec3{coordinate[0], coordinate[1], coordinate[2]});
            }
        }
       
    }
    for (Edge e: mesh->edges())
    {
        Vector3 v1 = geometry->inputVertexPositions[e.halfedge().vertex()];
        Vector3 v2 = geometry->inputVertexPositions[e.halfedge().twin().vertex()];
        for (size_t i = 1; i < subdiv_level - 1; i++)
        {
            subdiv_points[int_count + edge_count + eInd[e]] = ((double)i)/subdiv_level * v1 + (1. - ((double)i)/subdiv_level) * v2;
        }
    }
    for (Vertex v: mesh->vertices())
    {
        subdiv_points[int_count + edge_count + vInd[v]] = geometry->inputVertexPositions[v];
    }
    // initialization of "correct" lengths
    for (Face f: mesh->faces())
    {

    }
   
}
void generateVisualization()
{
    // Visualization
    //EdgeData<double> initialGuess(*mesh);
    //EdgeData<double> finalSolution(*mesh);
    //for (Edge e : mesh->edges())
    //{
    //    initialGuess[e] = x_init[eInd[e]];
    //    finalSolution[e] = x[eInd[e]];
        //cout << x[eInd[e]] << endl;
        //cout << x_init[eInd[e]] << endl;
    //}
    //psMesh->addEdgeScalarQuantity("Initial Guess", x_init);
    psMesh->addEdgeScalarQuantity("Final Solution", *sol);
    psMesh->addVertexScalarQuantity("curvature",
            geometry->vertexGaussianCurvatures);
    cout << "No of subdivision points:" << subdiv_points.size();

    polyscope::registerPointCloud("Initial Subdivision", subdiv_points);
}
 

int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("Optimization");
  args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

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

  // Make sure a mesh name was given
  if (!inputFilename) {
    std::cerr << "Please specify a mesh file as argument" << std::endl;
    return EXIT_FAILURE;
  }

  // Initialize polyscope
  cout << "Initialized" << endl;
  polyscope::init();

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  cout << "starting optimization";
  generateConstraints();
  subdivision();
  generateVisualization();
  //generateSVGs();
  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
