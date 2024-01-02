#pragma once

#include <iostream>
#include <memory>
#include <tuple>
#include <vector>
using std::cout;
using std::endl;
using std::pair;
using std::shared_ptr;
using std::tuple;
using std::vector;

#include "Common.h"

typedef tuple<int, int, double> T;
// a member of type [A,b] represents the linear system Ax = b
typedef pair<vector<T>, vector<double>> SparseSystem;
class IntrinsicFlattening {
  public:
    // constructor from geometry-central data
    // l[e] is the target edge length
    IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,
                        EdgeData<double> l);

    pair<CornerData<double>, CornerData<double>>
    CoherentAngleSystem(VertexData<double> targetCurvatures,
                        CornerData<double> targetBetas);
    CornerData<double> SolveIntrinsicOnly();
    void CheckConstraintsIntrinsicOnly(CornerData<double> beta);
    // TODO: think about this
    /* pair<EdgeData<double>, CornerData<double>>
    solveFromPlane(double flatWeight);
    */

  private:
    static constexpr double EPS = 1e-6;
    // convenience functions
    void nasoqTest();
    Eigen::VectorXd
    QPSolve(Eigen::SparseMatrix<double, Eigen::ColMajor, int> &A,
            Eigen::Matrix<double, Eigen::Dynamic, 1> &b,
            Eigen::SparseMatrix<double, Eigen::ColMajor, int> &C,
            Eigen::Matrix<double, Eigen::Dynamic, 1> &d,
            Eigen::SparseMatrix<double, Eigen::ColMajor, int> &E,
            Eigen::Matrix<double, Eigen::Dynamic, 1> &f);
    void addTriples(vector<Eigen::Triplet<double>> &triples, vector<T> &tuples,
                    int i = 0, int j = 0);
    // convenience function that concatenates a variable number of
    // std::vector<double> into an Eigen Vector
    template <typename... Args>
    Eigen::VectorXd concat(size_t size, Args &...args);

    // geometric helper functions
    double CornerAngle(double l_ij, double l_jk, double l_ki);
    VertexData<double> ComputeTargetCurvatures();

    // constraints used for the optimization
    Eigen::SparseMatrix<double, Eigen::ColMajor, int>
    constructMatrix(vector<Eigen::Triplet<double>> &triples, int m, int n);
    // constraints used for intrinsic only optimization
    SparseSystem CATValidityConstraint();
    SparseSystem VertexAngleSumConstraint(VertexData<double> curvatures);
    SparseSystem OffsetConstraints();

    bool CheckCATValidityConstraint(CornerData<double> beta);
    bool CheckVertexAngleSumConstraint(CornerData<double> beta,
                                       VertexData<double> curvatures);

    // constraints used for CAS optimization
    SparseSystem PositiveAngleConstraint();
    SparseSystem FaceAngleSumConstraint();
    SparseSystem EdgeDelaunayConstraint();
    SparseSystem EdgeIntersectionAngleConstraint();

    // TODO: think about what the signature for these should be
    bool CheckPositiveAngleConstraint(CornerData<double> beta);
    bool CheckFaceAngleSumConstraint(CornerData<double> beta);
    bool CheckEdgeDelaunayConstraint(CornerData<double> beta);
    bool CheckEdgeIntersectionAngleConstraint(CornerData<double> beta);

    // Matrices defining the energy
    SparseSystem AngleDeviationPenalty(CornerData<double> beta);

    // pointers to geometric data
    shared_ptr<ManifoldSurfaceMesh> mesh;
    // quantities to initialize from the mesh
    size_t nVertices;
    size_t nEdges;
    size_t nCorners;
    size_t nFaces;
    CornerData<size_t> c_;
    CornerData<double> cornerAngles;
    FaceData<size_t> f_;
    EdgeData<size_t> e_;
    VertexData<size_t> v_;
};
