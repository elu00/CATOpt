#ifndef SOLVER_H
#define SOLVER_H
#include <stdlib.h>
#include <functional>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include "math.h"
#include "geometrycentral/numerical/linear_solvers.h"
#define SCP 0
#define LSCM 1
#define CIRCLE_PATTERNS 2
#define CETM 3
#define GRAD_DESCENT 0
#define NEWTON 1
#define TRUST_REGION 2
#define LBFGS 3
using namespace std::placeholders;

struct MeshHandle {
    // typedefs
    typedef std::function<void(double&, const Eigen::VectorXd&)> ComputeEnergy;
    typedef std::function<void(Eigen::VectorXd&, const Eigen::VectorXd&)> ComputeGradient;
    typedef std::function<void(Eigen::SparseMatrix<double>&, const Eigen::VectorXd&)> ComputeHessian;
    
    // constructor
    MeshHandle() {}
    
    // member variables
    ComputeEnergy computeEnergy;
    ComputeGradient computeGradient;
    ComputeHessian computeHessian;
};

class SolverF {
public:
    // constructor
    SolverF();
    
    // gradient descent
    void gradientDescent();
    
    // newton
    void newton();
    
    // trust region
    void trustRegion();
    
    // lbfgs
    void lbfgs(int m = 10);
    
    // member variables
    MeshHandle *handle;
    Eigen::VectorXd x;
    std::vector<double> obj;
    int n;
};

#endif
