#include "CatOpt.h"
#include "fusion.h"


using namespace mosek::fusion;
using namespace monty;
void CatOpt::initializeQuantities() {
    // Initialization
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    nFaces = mesh->nFaces();
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    cInd = mesh->getCornerIndices();
    fInd = mesh->getFaceIndices();
    geometry->requireEdgeLengths();
    geometry->requireVertexGaussianCurvatures();
    geometry->requireVertexAngleSums();
    angleDefects = geometry->vertexGaussianCurvatures;
    sol = vector<double>(nEdges);
}

// calls mosek to generate optimal alphas, and places them into sol
void CatOpt::generateConstraints() {
    // Model initialization
    // TODO: check for optimality
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", nEdges, Domain::inRange(-2 * PI, 2 * PI));

    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    // Equality constraint initialization
    rhs = vector<double>(nVertices);
    for (size_t i = 0; i < nVertices; i++) {
        if (mesh-> vertex(i).isBoundary()) {
            rhs[i] = PI - geometry->vertexAngleSums[mesh->vertex(i)];
        } else {
            rhs[i] =  angleDefects[mesh->vertex(i)];                
        }
    }
    for (Edge e : mesh->edges()) {
        double weight = e.isBoundary() ? 1. : 2. ;
        rows.emplace_back(vInd[e.halfedge().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(weight);
        rows.emplace_back(vInd[e.halfedge().twin().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(weight);
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto Meq = Matrix::sparse(nVertices, nEdges, r, c, v);
    auto eqRHS = new_array_ptr(rhs);
    M->constraint("eq constraints", Expr::mul(Meq, x), Domain::equalsTo(eqRHS));
    // inequality constraints
    ineqRHS0 = vector<double>(nCorners);
    ineqRHS1 = vector<double>(nCorners);
    rows.clear();
    cols.clear();
    values.clear();
    for (Corner c : mesh->corners()) {
        ineqRHS0[cInd[c]] = -geometry->cornerAngle(c);
        ineqRHS1[cInd[c]] = 2 * PI - geometry->cornerAngle(c);
        Halfedge h = c.halfedge();
        rows.emplace_back(cInd[c]);
        cols.emplace_back(eInd[h.edge()]);
        values.emplace_back(1.);
        rows.emplace_back(cInd[h.next().corner()]);
        cols.emplace_back(eInd[h.edge()]);
        values.emplace_back(1.);
    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto Mineq = Matrix::sparse(nCorners, nEdges, r, c, v);
    auto inRHS0 = new_array_ptr(ineqRHS0);
    auto inRHS1 = new_array_ptr(ineqRHS1);
    // sum is greater than 0
    M->constraint("ineq0 constraints", Expr::mul(Mineq, x), Domain::greaterThan(inRHS0));
    // sum is less than 2pi
    M->constraint("ineq1 constraints", Expr::mul(Mineq, x), Domain::lessThan(inRHS1));

    auto ones = std::make_shared<ndarray<double, 1>>(shape(nEdges), 1.);
    //cout << "ineq generated" << endl;
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::vstack(t, x), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << M->getProblemStatus() << endl;
    //cout << x->level() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = x->getSize();
    auto xVal = x->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    for (int i = 0; i < xsize; ++i) {
        sol[i] = (*xVal)[i];
    }
}
