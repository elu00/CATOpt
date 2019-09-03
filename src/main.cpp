#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "Eigen/Sparse"
#include "Eigen/Dense"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using namespace Eigen;
using std::cout;
using std::endl;
using std::unique_ptr;

// == Geometry-central data
unique_ptr<HalfedgeMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

// Optimization Parameters
double t = 2.;
double mew = 1.5;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

SparseMatrix<double> generateConstraints()
{
    // Initialization
    size_t nVertices = mesh->nVertices;
    size_t nEdges = mesh->nEdges;
    EdgeData<size_t> eInd = mesh->getEdgeIndices();
    VertexData<size_t> vInd = mesh->getVertexIndices();
    geometry->requireVertexGaussianCurvatures();
    VertexData<double> angleDefects = geometry->vertexGaussianCurvatures;

    // Matrix initialization
    SparseMatrix<double> d0 = Eigen::SparseMatrix<double>(nVertices, nEdges);
    std::vector<Triplet<double>> tripletList;
    VectorXd rhs = VectorXd(nVertices);
    for (size_t i = 0; i < nVertices; i++) {
      rhs[i] = angleDefects[mesh->vertex(i)]/2;
    }
    for (Edge e : mesh->edges()) {
      tripletList.emplace_back(vInd[e.halfedge.vertex()], eInd[e], 1);
      tripletList.emplace_back(vInd[e.halfedge.twin.vertex()], eInd[e], 1);
    }
    d0.setFromTriplets(tripletList.begin(), tripletList.end());

    return d0;
    /*
    SparseLU<SparseMatrix<double>> solver;
    solver.compute(d0);
    if (solver.info() != Eigen::Success) {
      cout << "solving failed" << endl;
    }
    Vector<double> solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
      cout << "solving failed";
    }
    */

}

inline double logBarrier(const double x, const double upperBound)
{
  return -(1/t) * log(-(upperBound - x));
}

inline void vecProj(Eigen::SparseMatrix<double> A, Eigen::VectorXd x, Eigen::VectorXd b)
{
  
}
void objFunction()
{
    
}


int main(int argc, char **argv) {

  // Configure the argument parser
  args::ArgumentParser parser("geometry-central & Polyscope example project");
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
  polyscope::init();

  // Load mesh
  std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));
  SparseMatrix<double> constraints = generateConstraints();

  // Register the mesh with polyscope
  psMesh = polyscope::registerSurfaceMesh(
      polyscope::guessNiceNameFromPath(args::get(inputFilename)),
      geometry->inputVertexPositions, mesh->getFaceVertexList(),
      polyscopePermutations(*mesh));

  // Give control to the polyscope gui
  polyscope::show();

  return EXIT_SUCCESS;
}
