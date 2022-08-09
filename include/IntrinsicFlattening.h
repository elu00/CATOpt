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
        IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,shared_ptr<VertexPositionGeometry> geometry);
        pair<CornerData<double>, CornerData<double>> CoherentAngleSystem(VertexData<double> targetCurvatures, CornerData<double> targetBetas);
        CornerData<double> solveIntrinsicOnly();
        pair<CornerData<double>, CornerData<double>> solveFromPlane(double flatWeight);
    private:
        // convenience function
        void shiftTriples(vector<T>& tripletList, int i, int j); 
        pair<vector<T>, vector<double>> PositiveAngleConstraint();
        pair<vector<T>, vector<double>> FaceAngleSumConstraint();
        pair<vector<T>, vector<double>> VertexAngleSumConstraint(VertexData<double> curvatures);
        pair<vector<T>, vector<double>> EdgeDelaunayConstraint();
        pair<vector<T>, vector<double>> EdgeIntersectionAngleConstraint();
        pair<vector<T>, vector<double>> CATValidityConstraint();
        pair<vector<T>, vector<double>> AngleDeviationPenalty(CornerData<double> beta);
        /*
        void buildOffsetConstraints(Model::t& M, Variable::t& alpha, Variable::t& beta);
        void buildBoundaryObjective(Model::t& M, Variable::t& a, Variable::t& t, size_t excl, double interpolationWeight);
        monty::rc_ptr<mosek::fusion::Matrix> sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values);
        */


        // pointers to geometric data
        shared_ptr<ManifoldSurfaceMesh> mesh;
        shared_ptr<VertexPositionGeometry> geometry;
        // quantities to read off from the mesh
        size_t nVertices;
        size_t nEdges;
        size_t nCorners;
        size_t nFaces;
        CornerData<size_t> c_;
        FaceData<size_t> f_;
        EdgeData<size_t> e_;
        VertexData<size_t> v_;
};
