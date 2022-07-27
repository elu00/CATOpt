#pragma once

#include <algorithm>
#include <tuple>
#include <map>
#include <vector>
#include <memory>
#include <tuple>
using std::cout;
using std::vector;
using std::shared_ptr;
using std::endl;
using std::tuple;


#include "Common.h"

#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

class IntrinsicFlattening {
    public:
        // constructor from geometry-central data
        IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,shared_ptr<VertexPositionGeometry> geometry);
        // 
        SolutionData solve();
        EdgeData<double> solveKSS();
        CornerData<double> solveIntrinsicOnly();
        SolutionData solveFromPlane(double flatWeight);
    private:
        // convenience function
        void buildOffsetConstraints(Model::t& M, Variable::t& alpha, Variable::t& beta);
        void buildIntersectionAngleConstraints(Model::t& M, CornerData<double>& beta, Variable::t& a);
        void buildFaceConstraints(Model::t& M, Variable::t& a);
        void buildVertexConstraints(Model::t& M, Variable::t& a);
        void buildDelaunayConstraints(Model::t& M, Variable::t& a);
        void buildBoundaryObjective(Model::t& M, Variable::t& a, Variable::t& t, size_t excl, double interpolationWeight);
        monty::rc_ptr<mosek::fusion::Matrix> sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values);


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
