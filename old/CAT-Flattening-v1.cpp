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

#include <nlopt.hpp>
#include <math.h>

using namespace geometrycentral;
using std::cerr;
using std::cout;
using std::endl;
using std::string;

double objWrapper(unsigned n, const double* x, double* grad, void* f_data);
void angleDefectWrapper(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data);
void validAngleWrapper(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data);
class CatData {
  // Initialized stuff
  Geometry<Euclidean>* geom;
  HalfedgeMesh* mesh;
  std::string niceName;

  // Original mesh information
  VertexData<size_t> vInd;
  size_t nVerts;
  size_t nHalfedges;
  size_t nCorners;
  size_t dim;
  VertexData<double> angleDefects;
  HalfedgeData<size_t> hInd;
  CornerData<size_t> cInd;
  //EdgeData<double> lengths;

public:
  // Derived Information
  HalfedgeData<double> theta;
  HalfedgeData<double> alpha;
  HalfedgeData<double> beta;

  double obj(unsigned n, const double* x, double* grad, void* f_data)
  {
    double accum = 0;
    for (size_t i = 0; i < n; i++)
    {
      accum += pow(x[i],2);
      if (grad) 
      {
        grad[i] = 2*x[i];
      }
    }
    return accum;
  }
  
  void angleDefectCalc(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data) 
  {
    // initialize angle defect
    for (size_t i = 0; i < m; i++)
    {
      result[i] = - angleDefects[i];
    }
    if (gradient) 
    {
      for (size_t i = 0; i < n*m ; i++)
      {
        gradient[i] = 0.;
      }
    }
    //go through halfedges
    for (EdgePtr e : mesh->edges()) 
    {
      HalfedgePtr he1 = e.halfedge();
      HalfedgePtr he2 = he1.twin();
      size_t h1 = hInd[he1];
      size_t h2 = hInd[he2];
      size_t v1 = vInd[he1.vertex()];
      size_t v2 = vInd[he2.vertex()];
      result[v1] += x[h1] + x[h2];
      result[v2] += x[h1] + x[h2];

      if (gradient) 
      {
        gradient[v1 * n + h1] = 1;
        gradient[v2 * n + h1] = 1;
        gradient[v1 * n + h2] = 1;
        gradient[v2 * n + h2] = 1;
      }
    }
    return;
  }
  void validAngleCalc(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data) 
  {
    // first coord specifies angle + alphas - 2pi < 0, second specifies -angle -alphas < 0 
    for (CornerPtr c: mesh->corners())
    {
      double curAngle =  geom->angle(c);
      result[2 * cInd[c]] = curAngle - 2 * M_PI;
      result[2 * cInd[c] + 1] = - curAngle;
    }
    for (CornerPtr c: mesh->corners())
    {
      size_t h = hInd[c.halfedge()];
      result[2 * cInd[c]] += x[h];
      result[2 * cInd[c] + 1] -= x[h];
      result[2 * cInd[c.next()]] += x[h];
      result[2 * cInd[c.next()] + 1] -= x[h];
      if (gradient) 
      {
        gradient[2 * cInd[c] * n + h] = 1;
        gradient[(2 * cInd[c] + 1) * n + h] = -1;
        gradient[2 * cInd[c.next()] * n + h] = 1;
        gradient[(2 * cInd[c.next()] + 1) * n + h] = -1;
      }
    }
    return;
  }

  void optimize() {
    nlopt::opt opt(nlopt::LN_COBYLA, nHalfedges);
    opt.set_min_objective(&objWrapper, this);
    opt.set_lower_bounds(-2 * M_PI);
    opt.set_upper_bounds(2 * M_PI);

    std::vector<double> tol1(nVerts, 1e-8);
    std::vector<double> tol2(2 * nCorners, 1e-8);

    opt.add_equality_mconstraint(&angleDefectWrapper, this, tol1);
    opt.add_inequality_mconstraint(&validAngleWrapper, this, tol2);

    //opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
    //opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
    opt.set_xtol_rel(1e-4);
	std::vector<double> x(nHalfedges, 0);
    double minf;
    try {
      nlopt::result result = opt.optimize(x, minf);
      std::cout << "found minimum" << std::setprecision(10) << minf << std::endl;
    } catch (std::exception& e) {
      std::cout << "nlopt failed: " << e.what() << std::endl;
    }
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
    nCorners = mesh->nCorners();
    dim = nVerts + nHalfedges;

    theta = HalfedgeData<double>(mesh);
    alpha = HalfedgeData<double>(mesh);
    beta = HalfedgeData<double>(mesh);

    geom->getVertexAngleDefects(angleDefects);
    optimize();
    delete geom;
    delete mesh;
  }
};

double objWrapper(unsigned n, const double* x, double* grad, void* f_data)
{
  return static_cast<CatData*>(f_data)->obj(n, x, grad, NULL);
}

void angleDefectWrapper(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data) 
{
  static_cast<CatData*>(func_data)->angleDefectCalc(m, result, n, x, gradient, NULL);
}

void validAngleWrapper(unsigned m, double* result, unsigned n, const double* x, double* gradient, void* func_data)
{
  static_cast<CatData*>(func_data)->validAngleCalc(m, result, n, x, gradient, NULL);
}

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



/* class CatDataOld {
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
  EdgeData<double> lengths;

  public:
    // Derived Information
    HalfedgeData<double> finalCurvature;
    VertexData<double> multiplier;

    HalfedgeData<double> theta;
    HalfedgeData<double> d;
    HalfedgeData<double> alpha;
    HalfedgeData<double> beta;

    EdgeData<char> badEdges;
    EdgeData<double> netEdgeCurvature;
    EdgeData<char> negEdges;

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
        rhs[i] = 2 * angleDefects[mesh->vertex(i - nHalfedges)];
      }
      for (EdgePtr e : mesh->edges()) {
        HalfedgePtr h1 = e.halfedge();
        HalfedgePtr h2 = h1.twin();
        size_t v1 = vInd[h1.vertex()];
        size_t v2 = vInd[h2.vertex()];
        tripletList.emplace_back(nHalfedges + v1, hInd[h1], lengths[e]);
        tripletList.emplace_back(nHalfedges + v1, hInd[h2], lengths[e]);
        tripletList.emplace_back(nHalfedges + v2, hInd[h1], lengths[e]);
        tripletList.emplace_back(nHalfedges + v2, hInd[h2], lengths[e]);

        tripletList.emplace_back(hInd[h1], nHalfedges + v1, lengths[e]);
        tripletList.emplace_back(hInd[h2], nHalfedges + v1, lengths[e]);
        tripletList.emplace_back(hInd[h1], nHalfedges + v2, lengths[e]);
        tripletList.emplace_back(hInd[h2], nHalfedges + v2, lengths[e]);
      }
      d0.setFromTriplets(tripletList.begin(), tripletList.end());
      // cout << "Matrix built" << endl;

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
      }
      for (size_t i = nHalfedges; i < dim; i++) {
        multiplier[i - nHalfedges] = solution[i];
      }
      polyscope::getSurfaceMesh(niceName)->addQuantity("Original Curvature", angleDefects);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Curvature change", finalCurvature);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Lagrange Multiplier", multiplier);
      return;
    }
    // Updates the straight distances between edges, then checks for bad distances
    void updateDistances() {
      size_t bad_halfedges = 0;
      // Initialize all distances first
      for (HalfedgePtr h : mesh->allHalfedges()) {
        theta[h] = lengths[h.edge()] * finalCurvature[h];
        // Basic constraint
        if (theta[h] < 2. * M_PI && theta[h] > -2. * M_PI) {
          d[h] = (theta[h] == 0 ? lengths[h.edge()] : 2 * sin(theta[h] / 2) / finalCurvature[h]);
        } else {
          bad_halfedges++;
          badEdges[h.edge()] = true;
        }
      }
      cout << "Bad angles:" << bad_halfedges << endl;
    }
    // Updates angles based on distances, checks for self intersection
    void updateAngles() {
      size_t bad_halfedges = 0;
      for (HalfedgePtr h : mesh->allHalfedges()) {
        double ij = d[h];
        double jk = d[h.next()];
        double ki = d[h.next().next()];
        double cosAngle = (pow(ij, 2) + pow(ki, 2) - pow(jk, 2)) / (2 * ij * ki);

        if (cosAngle > 1.) {
          alpha[h] = 0;
          bad_halfedges++;
          badEdges[h.edge()] = true;
        } else if (cosAngle < -1.) {
          alpha[h] = M_PI;
          // alpha[h] = 10000;
          bad_halfedges++;
          badEdges[h.edge()] = true;
        } else {
          alpha[h] = acos(cosAngle);
        }
      }
      cout << "Bad distances:" << bad_halfedges << endl;
      for (VertexPtr v : mesh->vertices()) {
        angleDefects[v] = 2 * M_PI;
      }

      bad_halfedges = 0;
      for (HalfedgePtr h : mesh->allHalfedges()) {
        double betaA = alpha[h] + (theta[h] + theta[h.next().next()]) / 2.;
        angleDefects[h.vertex()] -= betaA;
        beta[h] = betaA;
        if (betaA > 2 * M_PI || betaA < 0.) {
          // cout << betaA << endl;
          bad_halfedges++;
          badEdges[h.edge()] = true;
        }
      }
      cout << "Self intersection:" << bad_halfedges << endl;
    }

    double averageAngleDefect() {
      double accum = 0;
      for (size_t i = 0; i < nVerts; i++) {
        accum += abs(angleDefects[i]);
      }
      return accum / nVerts;
    }

    void checkNegedges() {
      size_t neg_edges = 0;
      for (EdgePtr E : mesh->edges()) {
        netEdgeCurvature[E] = finalCurvature[E.halfedge()] + finalCurvature[E.halfedge().twin()];
        if (netEdgeCurvature[E] < 0) {
          negEdges[E] = true;
          neg_edges++;
        }
      }
      cout << "Neg edges:" << neg_edges << endl;
    }
    CatDataOld(std::string filename) {
      niceName = polyscope::utilities::guessNiceNameFromPath(filename);
      mesh = new HalfedgeMesh(PolygonSoupMesh(filename), geom);
      polyscope::registerSurfaceMesh(niceName, geom);

      vInd = mesh->getVertexIndices();
      nVerts = mesh->nVertices();
      hInd = mesh->getHalfedgeIndices();
      nHalfedges = mesh->nHalfedges();
      dim = nVerts + nHalfedges;

      finalCurvature = HalfedgeData<double>(mesh);
      multiplier = VertexData<double>(mesh);
      theta = HalfedgeData<double>(mesh);
      d = HalfedgeData<double>(mesh);
      alpha = HalfedgeData<double>(mesh);
      beta = HalfedgeData<double>(mesh);
      badEdges = EdgeData<char>(mesh, false);
      netEdgeCurvature = EdgeData<double>(mesh);
      negEdges = EdgeData<char>(mesh, false);

      geom->getVertexAngleDefects(angleDefects);
      geom->getEdgeLengths(lengths);
      for (size_t i = 0; i < 1000; i++) {
        cout << "Starting Iteration " << i << endl;
        solveOptMatrix();
        updateDistances();
        updateAngles();
        cout << "Average angle defect: " << averageAngleDefect() << endl;
        cout << "Done" << endl;
      }
      checkNegedges();
      polyscope::getSurfaceMesh(niceName)->addQuantity("Central angles", theta);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Interior angles", alpha);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Straight distances", d);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Exterior angles", beta);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Final Angle Defect", angleDefects);
      polyscope::getSurfaceMesh(niceName)->addQuantity("Net Edge Curvature", netEdgeCurvature);
      polyscope::getSurfaceMesh(niceName)->addSubsetQuantity("Bad edges", badEdges);
      polyscope::getSurfaceMesh(niceName)->addSubsetQuantity("Neg edges", negEdges);
      delete geom;
      delete mesh;
    }
}; */