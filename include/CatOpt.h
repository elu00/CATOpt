#pragma once
#include "geometrycentral/numerical/linear_solvers.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "imgui.h"

#include "glm/vec3.hpp"

using std::max;
using std::shared_ptr;

using Eigen::VectorXd;
using std::make_tuple;
using std::string;
using std::tuple;
using std::vector;

class CatOpt {
  public:
    void conformalFlatten();
    void generateConstraints();
    CatOpt(string filename);
    void polyscopeInit();
    // == Geometry-central data
    string inputMeshPath;

    shared_ptr<ManifoldSurfaceMesh> CATmesh;
    shared_ptr<EdgeLengthGeometry> intrinsicGeometry;

    shared_ptr<ManifoldSurfaceMesh> flatmesh;
    shared_ptr<VertexPositionGeometry> flatGeometry;

    // unique_ptr<VertexPositionGeometry> CATgeometry;

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
    // Vector<double> x_init;
    vector<double> sol;
    vector<double> rhs;
    vector<double> ineqRHS0;
    vector<double> ineqRHS1;

    size_t iter = 0;
    // specifies number of
    size_t subdiv_level = 5;
    vector<Vector3> subdiv_points;
    // vector<Vector3> grad;
    vector<tuple<size_t, size_t, double>> correct_dist;
    Eigen::SparseMatrix<double> bendingMatrix;
    // tuning parameters for gradient descent
    // maybe change this back to 0.1
    double lineSearchAlpha = 0.1;
    double beta = 0.5;
    double ep = 1e-3;
    // weight for the bending energy
    double bendingWeight = 1e-8;

    // stuff for bff
    vector<Vector2> flattened;
    vector<double> alphas;
    CornerData<double> targetAngles;

    // Polyscope visualization handle, to quickly add data to the surface
    polyscope::SurfaceMesh *psMesh;
    polyscope::SurfaceMesh *CATpsMesh;
    void circlePatterns();

    void setOffsets();

  private:
    double sqr(double x) { return x * x; }
    // convenience function to return square of norm of gradient
    double grad_norm_sq(const VectorXd &grad1, const VectorXd &grad2,
                        const VectorXd &grad3) {
        return grad1.squaredNorm() + grad2.squaredNorm() + grad3.squaredNorm();
    }

    Vector3 bary(Face f, double a, double b, double c);

    // Intrinsic angle optimization stuff
    // Flattening optimization stuff

    void dbgSVG(string filename);
    double confObjective(vector<Vector2> &flattened, vector<double> &alphas);
    void confGradient(vector<Vector2> &x, vector<double> &alphas,
                      vector<Vector2> &flattened_grad,
                      vector<double> &alphas_grad);
    void confStep(int n);
    void loadModel(const std::string &inputPath, bff::Model &model,
                   std::vector<bool> &surfaceIsClosed);

    void flatten(bff::Model &model, const std::vector<bool> &surfaceIsClosed,
                 int nCones, bool flattenToDisk, bool mapToSphere);
    void buildNewGeometry();
    void testFlatteningDerivatives();
    // just for validating the SVG formula I'm using
    void testSVG();

    // Circle pattern stuff?
};
