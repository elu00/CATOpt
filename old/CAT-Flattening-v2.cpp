#include "polyscope/polyscope.h"

#include <iostream>

#include "geometrycentral/geometry.h"
#include "geometrycentral/halfedge_mesh.h"
#include "geometrycentral/linear_solvers.h"
#include "geometrycentral/polygon_soup_mesh.h"

#include <Eigen/SparseLU>

#include "args/args.hxx"
#include "json/json.hpp"

#define GLM_ENABLE_EXPERIMENTAL
#include "glm/gtx/string_cast.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


using namespace geometrycentral;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

class CatData {
  // Initialized stuff
  Geometry<Euclidean>* geom;
  HalfedgeMesh* mesh;
  std::string niceName;

  // Original mesh information
  VertexData<size_t> vInd;
  size_t nVerts;
  size_t nHalfedges;
  size_t dim;
  VertexData<double> angleDefects;
  HalfedgeData<size_t> hInd;
  CornerData<size_t> cInd;

public:
  // Derived Information
  HalfedgeData<double> finalCurvature;
  VertexData<double> multiplier;

  HalfedgeData<double> d;
  HalfedgeData<double> alpha;

  EdgeData<double> netEdgeCurvature;
  CornerData<double> cornerAngle;
  CornerData<char> badCorners;

  // Solves the optimization problem
  void solveOptMatrix() {
    Eigen::SparseMatrix<double> d0 = Eigen::SparseMatrix<double>(dim, dim);
    std::vector<Eigen::Triplet<double>> tripletList;
    Vector<double> rhs = Vector<double>(dim);
    // cout << dim << endl;
    for (size_t i = 0; i < nHalfedges; i++) {
      tripletList.emplace_back(i, i, 1.);
      rhs[i] = 0.;
    }

    for (size_t i = nHalfedges; i < dim; i++) {
      rhs[i] = angleDefects[mesh->vertex(i - nHalfedges)];
    }
    for (EdgePtr e : mesh->edges()) {
      HalfedgePtr h1 = e.halfedge();
      HalfedgePtr h2 = h1.twin();
      size_t v1 = vInd[h1.vertex()];
      size_t v2 = vInd[h2.vertex()];
      tripletList.emplace_back(nHalfedges + v1, hInd[h1], 1);
      tripletList.emplace_back(nHalfedges + v1, hInd[h2], 1);
      tripletList.emplace_back(nHalfedges + v2, hInd[h1], 1);
      tripletList.emplace_back(nHalfedges + v2, hInd[h2], 1);

      tripletList.emplace_back(hInd[h1], nHalfedges + v1, 1);
      tripletList.emplace_back(hInd[h2], nHalfedges + v1, 1);
      tripletList.emplace_back(hInd[h1], nHalfedges + v2, 1);
      tripletList.emplace_back(hInd[h2], nHalfedges + v2, 1);
    }
    d0.setFromTriplets(tripletList.begin(), tripletList.end());
    cout << "Matrix built" << endl;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(d0);
    if (solver.info() != Eigen::Success) {
      cout << "solving failed" << endl;
    }
    Vector<double> solution = solver.solve(rhs);
    if (solver.info() != Eigen::Success) {
      cout << "solving failed";
    }
    // cout << "Matrix solved";
    for (size_t i = 0; i < nHalfedges; i++) {
      finalCurvature[i] = solution[i];
      //cout << solution[i] << endl;
    }
    for (size_t i = nHalfedges; i < dim; i++) {
      multiplier[i - nHalfedges] = solution[i];
    }
    polyscope::getSurfaceMesh(niceName)->addQuantity("Original Curvature", angleDefects);
    polyscope::getSurfaceMesh(niceName)->addQuantity("Curvature change", finalCurvature);
    polyscope::getSurfaceMesh(niceName)->addQuantity("Lagrange Multiplier", multiplier);
    return;
  }
  void checkAngles()
  {
    size_t bad = 0;
    for (CornerPtr c: mesh->corners())
    {
      cornerAngle[c] =  geom->angle(c);
    }
    for (CornerPtr c: mesh->corners())
    {
      HalfedgePtr h = c.halfedge();
      cornerAngle[c] += finalCurvature[h];
      cornerAngle[c.next()] += finalCurvature[h];
    }
    cout << bad << endl;
    for (CornerPtr c: mesh->corners())
    {
      if (cornerAngle[c] < 0 || cornerAngle[c] > 2 * M_PI)
      {
        bad++;
        badCorners[c] = true;
      }
    }
    cout << bad << endl;
  }
  CatData(std::string filename) {
    niceName = polyscope::utilities::guessNiceNameFromPath(filename);
    mesh = new HalfedgeMesh(PolygonSoupMesh(filename), geom);
    polyscope::registerSurfaceMesh(niceName, geom);

    vInd = mesh->getVertexIndices();
    nVerts = mesh->nVertices();
    hInd = mesh->getHalfedgeIndices();
    nHalfedges = mesh->nHalfedges();
    cInd = mesh->getCornerIndices();
    dim = nVerts + nHalfedges;

    finalCurvature = HalfedgeData<double>(mesh);
    multiplier = VertexData<double>(mesh);
    d = HalfedgeData<double>(mesh);
    alpha = HalfedgeData<double>(mesh);
    netEdgeCurvature = EdgeData<double>(mesh);
    cornerAngle = CornerData<double>(mesh);
    badCorners = CornerData<char>(mesh, false);

    geom->getVertexAngleDefects(angleDefects);
    solveOptMatrix();
    checkAngles();
    polyscope::getSurfaceMesh(niceName)->addQuantity("Central angles", alpha);
    delete geom;
    delete mesh;
  }
};


int main(int argc, char** argv) {
  // Configure the argument parser
  /*args::ArgumentParser parser("A simple demo of Polyscope.\nBy "
                              "Nick Sharp (nsharp@cs.cmu.edu)",
                              "");
  args::PositionalList<string> files(parser, "files", "One or more files to visualize");
  */
  // Options
  polyscope::options::autocenterStructures = true;
  // Initialize polyscope
  polyscope::init();
  CatData* c = new CatData("C:/spot1.obj");
  delete c;
  // Show the gui
  polyscope::show();

  return 0;
}
