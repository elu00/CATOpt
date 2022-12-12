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

#include "unsupported/Eigen/LevenbergMarquardt"

class EmbeddingOptimization {
    public:
        EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry, CornerData<double> beta);
        std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry>> solve(int N);
        void optimize();
        // Optimization procedures
        void evaluateEnergy(const Eigen::VectorXd& v, Eigen::VectorXd& energy);
        void evaluateJacobian(const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& J);
        Eigen::VectorXd gradientDescent();
        Eigen::VectorXd x;


    private:
        // pointers to geometric data associated to the original mesh
        int n;
        shared_ptr<ManifoldSurfaceMesh> mesh;
        shared_ptr<VertexPositionGeometry> geometry;
        CornerData<double> beta;
        // quantities to read off from the mesh
        size_t nVertices;
        size_t nEdges;
        size_t nCorners;
        size_t nFaces;
        CornerData<size_t> c_;
        FaceData<size_t> f_;
        EdgeData<size_t> e_;
        VertexData<size_t> v_;


        // pointers to geometric data associated to the subdivided mesh
        //
        shared_ptr<ManifoldSurfaceMesh> submesh;
        shared_ptr<VertexPositionGeometry> subgeometry;
        vector<double> c_iso_0;
        vector<double> c_iso_1;
        vector<double> c_iso_2;
        size_t nSubdividedVertices;


        // Union find functions and data structures
        vector<int> top;
        vector<int> next;
        void merge(int a, int b);
        int find(int a);
        void buildEquivalenceClasses();

        // final indices built by reindexing union find stuff
        vector<int> finalIndices;
        int buildFinalIndices();

        // Subdivision routines
        void buildSubdivision();
        void buildIntrinsicCheckerboard();

        
        void basisFunctionDebugging();

        // barycentric coordinates for each quad
        Vector3 bary(Corner c, int x, int y);
        tuple<double, double, double> baryCoords(int x, int y);
};
struct EnergyFunctor : Eigen::SparseFunctor<double,int>
{
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorType;
    typedef Eigen::SparseFunctor<double,int> Base;
    typedef typename Base::JacobianType JacobianType;
    EnergyFunctor(int inputs, int values, EmbeddingOptimization* _E) : SparseFunctor<double,int>(inputs,values) { 
        E = _E;
    }

    int operator()(const VectorType& v, VectorType& energy) {
        E->evaluateEnergy(v, energy);
        return 0;
    };
    int df(const VectorType& uv, JacobianType& fjac){ 
        E->evaluateJacobian(uv, fjac);
        cout << "JACOBIAN NORM:" << fjac.norm() << endl;
        return 0;
    };
    EmbeddingOptimization* E;
};



