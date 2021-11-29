#ifndef CIRCLE_PATTERNS_H
#define CIRCLE_PATTERNS_H

#include "Solver.h"
#include "Common.h"
#include <stack>
#include <string>

using std::vector;
using std::shared_ptr;
using std::endl;

class CirclePatterns{
public:
    shared_ptr<ManifoldSurfaceMesh> mesh;
    // constructor
    CirclePatterns(shared_ptr<ManifoldSurfaceMesh> mesh0, Vertex infVertex,
    EdgeData<bool> eMask, EdgeData<bool> eBdry, FaceData<bool> fMask,
    int optScheme0, Eigen::VectorXd thetas);
    
    // parameterize
    VertexData<Eigen::Vector2d> parameterize();


    
protected:

        
    // sets thetas
    //void setThetas();
    
    // compute angles
    //bool computeAngles();
    
    // computes energy, gradient and hessian
    void computeEnergy(double& energy, const Eigen::VectorXd& rho);
    void computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho);
    void computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho);

    // sets radii
    void setRadii();
    
    // compute radii
    bool computeRadii();
    
    // computes angles and edge lengths
    void computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths);
    
    // determines position of unfixed face vertex
    void performFaceLayout(Halfedge he, const Eigen::Vector2d& dir, Eigen::VectorXd& lengths,
                           std::unordered_map<int, bool>& visited, std::stack<Edge>& stack);

    
    // sets uvs
    void setUVs();

    void normalize();
    double uvArea(Face f);


    Eigen::Vector2d uvBarycenter(Face f);

    
    // member variables
    //Vector<double> sol;
    Eigen::VectorXd angles;
    Eigen::VectorXd thetas;
    Eigen::VectorXd radii;
    Eigen::VectorXi eIntIndices;
    int imaginaryHe;
    SolverF solver;
    int OptScheme;

    Vertex infVertex;
    vector<double> sol;
    EdgeData<size_t> eInd;
    VertexData<size_t> vInd;
    FaceData<size_t> fInd;
    HalfedgeData<size_t> hInd;
    // whether or not each edge is considered in the solve
    EdgeData<bool> eMask;
    EdgeData<bool> eBdry;
    FaceData<bool> fMask;
    VertexData<Eigen::Vector2d> uv;

};

#endif 
