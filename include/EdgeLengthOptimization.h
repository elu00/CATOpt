#pragma once

#include <algorithm>
#include <tuple>
#include <map>
#include <iostream>
#include <vector>
#include <memory>
#include <tuple>
using std::cout;
using std::vector;
using std::shared_ptr;
using std::endl;
using std::tuple;


#include "Common.h"

class EdgeLengthOptimization {
    public:
        // initializes mesh with connectivity and desired intersection angles (etas) per edge
        EdgeLengthOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> eta);
        void solve();
        void planarLayout();
    private:
        int n;
        shared_ptr<ManifoldSurfaceMesh> mesh;
        shared_ptr<VertexPositionGeometry> geometry;
        EdgeData<double> eta;
        size_t nVertices;
        size_t nEdges;
        size_t nCorners;
        size_t nFaces;
        CornerData<size_t> c_;
        FaceData<size_t> f_;
        EdgeData<size_t> e_;
        VertexData<size_t> v_;

        double CornerAngle(double lij, double ljk, double lki);
        double ExpCornerAngle(double lambdaij, double lambdajk, double lambdaki);
        CornerData<double> InitializeCornerAngles(Eigen::VectorXd& lambda);
        double LengthEnergy(Eigen::VectorXd& lambda);
        Eigen::VectorXd LengthEnergyGradient(Eigen::VectorXd& lambda);
        Eigen::SparseMatrix<double> LengthEnergyHessian(Eigen::VectorXd& lambda);
        void MinimizeLengthEnergy();
        void PrescribeCurvature();
        void LayoutMesh();
};
