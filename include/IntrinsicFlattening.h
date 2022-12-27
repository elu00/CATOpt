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
using std::pair;
using std::tuple;


#include "Common.h"


typedef tuple<int, int, double> T;

class IntrinsicFlattening {
    public:
        // constructor from geometry-central data
        IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,EdgeData<double> l);

        pair<CornerData<double>, CornerData<double>> CoherentAngleSystem(VertexData<double> targetCurvatures, CornerData<double> targetBetas);
        CornerData<double> solveIntrinsicOnly();
        pair<EdgeData<double>, CornerData<double>> solveFromPlane(double flatWeight);
    private:
        // convenience function

        double CornerAngle(double l_ij, double l_jk, double l_ki);
        void nasoqTest();
        Eigen::VectorXd QPSolve( Eigen::SparseMatrix<double,Eigen::ColMajor,int>& A,
                Eigen::Matrix<double,Eigen::Dynamic,1>& b,
                Eigen::SparseMatrix<double,Eigen::ColMajor,int>& C,
                Eigen::Matrix<double,Eigen::Dynamic,1>& d,
                Eigen::SparseMatrix<double,Eigen::ColMajor,int>& E,
                Eigen::Matrix<double,Eigen::Dynamic,1>& f);
        void addTriples(vector<Eigen::Triplet<double>>& triples, vector<T>& tuples, int i = 0, int j = 0); 

        Eigen::SparseMatrix<double,Eigen::ColMajor,int> constructMatrix(vector<Eigen::Triplet<double>>& triples, int m, int n);
        pair<vector<T>, vector<double>> PositiveAngleConstraint();
        pair<vector<T>, vector<double>> FaceAngleSumConstraint();
        pair<vector<T>, vector<double>> VertexAngleSumConstraint(VertexData<double> curvatures);
        pair<vector<T>, vector<double>> EdgeDelaunayConstraint();
        pair<vector<T>, vector<double>> EdgeIntersectionAngleConstraint();
        pair<vector<T>, vector<double>> CATValidityConstraint();
        pair<vector<T>, vector<double>> AngleDeviationPenalty(CornerData<double> beta);
        pair<vector<T>, vector<double>> OffsetConstraints();


        // pointers to geometric data
        shared_ptr<ManifoldSurfaceMesh> mesh;
        // quantities to read off from the mesh
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
