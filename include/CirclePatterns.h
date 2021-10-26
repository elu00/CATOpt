#ifndef CIRCLE_PATTERNS_H
#define CIRCLE_PATTERNS_H

#include "Solver.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include <stack>
#include <string>

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::vector;
using std::shared_ptr;
using std::endl;

class CirclePatterns{
public:
    shared_ptr<ManifoldSurfaceMesh> mesh;
    // constructor
    CirclePatterns(shared_ptr<ManifoldSurfaceMesh> mesh0, int optScheme0, vector<double>& solve,
    EdgeData<size_t> eInd, VertexData<size_t> vInd, FaceData<size_t> fInd, Eigen::VectorXd thetas);
    
    // parameterize
    void parameterize();

    void dbgSVG(std::string filename);
    
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

    void setOffsets();

    Eigen::Vector2d uvBarycenter(Face f);

    
    // member variables
    Eigen::VectorXd angles;
    Eigen::VectorXd thetas;
    Eigen::VectorXd radii;
    Eigen::VectorXi eIntIndices;
    int imaginaryHe;
    SolverF solver;
    int OptScheme;

    vector<double> sol;
    EdgeData<size_t> eInd;
    VertexData<size_t> vInd;
    FaceData<size_t> fInd;
    VertexData<Eigen::Vector2d> uv;
    
};

#endif 
