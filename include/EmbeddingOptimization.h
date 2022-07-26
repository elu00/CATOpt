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

class EmbeddingOptimization {
    public:
        EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry, EdgeData<double> beta);
        std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry>> solve(int N);
    private:
        // pointers to geometric data
        int n;
        shared_ptr<ManifoldSurfaceMesh> mesh;
        shared_ptr<ManifoldSurfaceMesh> submesh;
        shared_ptr<VertexPositionGeometry> geometry;
        shared_ptr<VertexPositionGeometry> subgeometry;
        // quantities to read off from the mesh
        size_t nVertices;
        size_t nEdges;
        size_t nCorners;
        size_t nFaces;
        CornerData<size_t> c_;
        FaceData<size_t> f_;
        EdgeData<size_t> e_;
        VertexData<size_t> v_;


        // Union find functions and data structures
        vector<int> top;
        vector<int> next;
        void merge(int a, int b);
        int find(int a);
        void buildEquivalenceClasses();

        vector<int> finalIndices;
        int buildFinalIndices();

        void buildSubdivision();

        // energy stuff
        void evaluateEnergy(double& energy, const Eigen::VectorXd& v);
        void evaluateGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& v);
        Eigen::VectorXd gradientDescent();

        // barycentric coordinates for each quad
        Vector3 bary(Corner c, int x, int y);

        void buildIntrinsicCheckboard();
        monty::rc_ptr<mosek::fusion::Matrix> sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values);
};
