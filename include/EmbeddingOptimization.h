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
#include "polyscope/curve_network.h"

class EmbeddingOptimization {
    public:
        EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry, CornerData<double> beta);
        std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry>> initializeSubdivision(int N);
        void initializeLM();
        void optimizeOneStep(int MAX_ITERS = 1000);
        double fairnessNormalization;
    private:
        // main methods
        void evaluateEnergy(const Eigen::VectorXd& v, Eigen::VectorXd& energy);
        void evaluateJacobian(const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& J);


        void LMOneStep(int MAX_ITERS);
        // current status of the solution
        Eigen::VectorXd initialSolution;
        Eigen::VectorXd currentSolution;


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
        vector<tuple<size_t, size_t, size_t, size_t>> quads;
        vector<vector<size_t>> fairVertices;
        size_t nSubdividedVertices;


        // internal convenience functions
        double Angle(Vector2 u, Vector2 v); 
        std::tuple<double, double, double> bendAngles(double t1, double t2, double t3, double b1, double b2, double b3);
        std::tuple<Vector2,Vector2,Vector2> projectToPlane(Vector3 i, Vector3 j, Vector3 k);

        BezierTriangle Coefficients (Vector3 I, Vector3 J, Vector3 K, double Bi, double Bj, double Bk);
        BezierTriangle rotIndices(BezierTriangle T);

        Vector2 RationalBezierTriangle(BezierTriangle T, std::tuple<double,double,double> coords);


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
        bool intrinsicQuantitiesInitialized;

        // Optimization routines
        double sqr(double x);
        inline void addLengthTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, 
                size_t energyIndex, size_t iIndex, size_t jIndex, double target);
        inline void addLengthGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, 
                size_t energyIndex, size_t iIndex, size_t jIndex, double target);
        inline void addAngleTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, 
                size_t energyIndex, size_t iIndex, size_t jIndex, size_t kIndex, size_t lIndex, double target);
        inline void addAngleGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, 
                size_t energyIndex, size_t iIndex, size_t jIndex, size_t kIndex, size_t lIndex, double target);
        inline void addCenterTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, 
                size_t energyIndex, vector<size_t>& indices);
        inline void addCenterGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, 
                size_t energyIndex, vector<size_t>& indices);


        // optimization values
        bool LMInitialized;
        size_t LMInputs;
        size_t LMValues;
        
        void basisFunctionDebugging();

        // barycentric coordinates for each quad
        Vector3 bary(Corner c, int x, int y);
        tuple<double, double, double> baryCoords(int x, int y);
};
