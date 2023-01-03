#pragma once

#include <algorithm>
#include <tuple>
#include <map>
#include <iostream>
#include <vector>
#include <memory>
#include <tuple>
#include "geometrycentral/numerical/linear_solvers.h"
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
        // additional constructor for specifying initial guess at optimal edge lengths
        EdgeLengthOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l, EdgeData<double> eta);
        Eigen::VectorXd MinimizeLengthEnergy();
    private:
        // number of edges
        int n;
        shared_ptr<ManifoldSurfaceMesh> mesh;
        // initial guess for edge lengths
        Eigen::VectorXd lambda0;
        EdgeData<double> eta;
        size_t nEdges;
        EdgeData<size_t> e_;

        inline double Cl2(double x);
        double CornerAngle(double lij, double ljk, double lki);
        double ExpCornerAngle(double lambdaij, double lambdajk, double lambdaki);
        CornerData<double> InitializeCornerAngles(Eigen::VectorXd& lambda);
        double LengthEnergy(Eigen::VectorXd& lambda);
        Eigen::VectorXd LengthEnergyGradient(Eigen::VectorXd& lambda);
        Eigen::SparseMatrix<double> LengthEnergyHessian(Eigen::VectorXd& lambda);
};
