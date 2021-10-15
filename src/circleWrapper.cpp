#include "CatOpt.h"
#include "Solver.h"
#include "CirclePatterns.h"
#include "fusion.h"
using namespace mosek::fusion;
using namespace monty;

void CatOpt::circlePatterns() {
    targetAngles = CornerData<double>(*mesh);
    Eigen::VectorXd thetas(mesh->nEdges());

    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    Variable::t Theta = M->variable("Theta", mesh->nExteriorHalfedges(), Domain::inRange(0, PI));
    Variable::t alpha = M->variable("alpha", mesh->nCorners(), Domain::inRange(0, PI));
    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    vector<double> rhs(nEdges);
    // Interior agreement constraint
    for (Edge e : mesh->edges()) {
        if (!e.isBoundary()) {
            rows.emplace_back(eInd[e]);
            cols.emplace_back(eInd[e]);
            values.emplace_back(1);
            // TODO: calculate what this should be
            //double angle = geometry->cornerAngle(C);
            
            //cout << "Inital:" << angle << " Next:";
            //angle += sol[eInd[h.edge()]];
            //angle += sol[eInd[h.next().next().edge()]];
            //rhs[eInd[e]] = ;
        }
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto boundaryConst = Matrix::sparse(nEdges, nEdges, r, c, v);
    auto intThetas = new_array_ptr(rhs);
    M->constraint("Interior Agreement", Expr::mul(boundaryConst, alpha), Domain::equalsTo(intThetas));
    rows.clear();
    cols.clear();
    values.clear();
    // Intersection angle constraint
    rhs = vector<double>(nEdges, PI);
    for (Edge e : mesh->edges()) {
        Halfedge he = e.halfedge();
        if (he.isInterior()) {
            rows.emplace_back(eInd[e]);
            cols.emplace_back(cInd[he.next().next().corner()]);
            values.emplace_back(1);
        } 
        if (he.twin().isInterior()) {
            rows.emplace_back(eInd[e]);
            cols.emplace_back(cInd[he.twin().next().next().corner()]);
            values.emplace_back(1);
        }
    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto intersectionConst = Matrix::sparse(nEdges, mesh->nCorners(), r, c, v);
    auto intersectionPI = new_array_ptr(rhs);
    M->constraint("Intersection Agreement", Expr::add(Theta, Expr::mul(intersectionConst, alpha)), Domain::equalsTo(intersectionPI));
    rows.clear();
    cols.clear();
    values.clear();
    // Flat Boundary Constraint
    rows.clear();
    cols.clear();
    values.clear();
    // Sum to pi in each triangle constraint


    auto ones = std::make_shared<ndarray<double, 1>>(shape(nEdges), 1.);
    //Variable::t t = M->variable("t", 1, Domain::unbounded());
    //M->constraint(Expr::vstack(t, x), Domain::inQCone());
    //M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << "LP Done" << endl;
    cout << M->getProblemStatus() << endl;
    auto xsize = Theta->getSize();
    auto xVal = Theta->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    for (int i = 0; i < xsize; ++i) {
        thetas[i] = (*xVal)[i];
    }


    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        if (h.isInterior()) {
            

            //targetAngles[C] = angle;
        }
    }
    /*
    for (Vertex v: mesh->vertices()) {
        double accum = 0;
        for (Corner C : v.adjacentCorners()) {
            accum += targetAngles[C];
        }
        if (v.isBoundary()) {
            cout << "boundary" << accum << endl;
        } else {
            cout << "not boundary" << accum << endl;
        }
    }
    */
    CirclePatterns prob(mesh, 0, sol, eInd, vInd, fInd, thetas);
    cout << "starting parameterization" << endl;
    prob.parameterize();
    cout << "parameterization done" << endl;
    prob.dbgSVG("wog.svg");
}