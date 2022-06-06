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
using namespace monty;

class IntrinsicFlattening {
    public:
        IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,shared_ptr<VertexPositionGeometry> geometry);
        SolutionData solve();
        EdgeData<double> solveKSS();
        SolutionData solveFromPlane();
    private:
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
        // convenience function
        monty::rc_ptr<mosek::fusion::Matrix> sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, 
                    vector<double>& values);

};