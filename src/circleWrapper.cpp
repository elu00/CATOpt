#include "CatOpt.h"
#include "Solver.h"
#include "CirclePatterns.h"
#include "fusion.h"
using namespace mosek::fusion;
using namespace monty;

void CatOpt::circlePatterns() {
    //targetAngles = CornerData<double>(*mesh);
    EdgeData<size_t> eInd = flatmesh->getEdgeIndices();
    VertexData<size_t> vInd = flatmesh->getVertexIndices();
    CornerData<size_t> cInd = flatmesh->getCornerIndices();
    FaceData<size_t> fInd = flatmesh->getFaceIndices();
    size_t vi, vj, vk;
    for (Vertex v: flatmesh->vertices()) {
        if (v.isBoundary()) {
            int count = 0;
            for (Vertex u: v.adjacentVertices()) {
                if (u.isBoundary()) {
                    if (count == 0) vi = vInd[u];
                    if (count == 1) {
                        vk = vInd[u];
                        break;
                    }
                    count++;
                }
            }
            vj = vInd[v];
            //if (flatmesh->removeVertex(v) == Face()) cout << "BAD BOUNDARY" << endl;
            break;
        }
    }



    Eigen::VectorXd thetas(flatmesh->nEdges());

    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    Variable::t Theta = M->variable("Theta", flatmesh->nEdges(), Domain::inRange(-PI, PI));
    Variable::t alpha = M->variable("alpha", flatmesh->nCorners(), Domain::inRange(0, PI));
    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    vector<double> rhs(nEdges, 0);


    // Interior agreement constraint
    for (Edge e : flatmesh->edges()) {
        if (!e.isBoundary()) {
            rows.emplace_back(eInd[e]);
            cols.emplace_back(eInd[e]);
            values.emplace_back(1);
            rhs[eInd[e]] = PI - flatGeometry->cornerAngle(e.halfedge().next().next().corner()) 
                            - flatGeometry->cornerAngle(e.halfedge().twin().next().next().corner());
            //cout << rhs[eInd[e]] << endl;
        }
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto boundaryConst = Matrix::sparse(nEdges, nEdges, r, c, v);
    auto intThetas = new_array_ptr(rhs);
    M->constraint("Interior Agreement", Expr::mul(boundaryConst, Theta), Domain::equalsTo(intThetas));
    rows.clear();
    cols.clear();
    values.clear();



    // Intersection angle constraint
    rhs = vector<double>(nEdges, PI);
    for (Edge e : flatmesh->edges()) {
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
    auto intersectionConst = Matrix::sparse(nEdges, flatmesh->nCorners(), r, c, v);
    auto intersectionPI = new_array_ptr(rhs);
    M->constraint("Intersection Agreement", Expr::add(Theta, Expr::mul(intersectionConst, alpha)), Domain::equalsTo(intersectionPI));
    rows.clear();
    cols.clear();
    values.clear();




    // Flat Boundary Constraint
    size_t count = 0;
    for (Vertex v: flatmesh->boundaryLoop(0).adjacentVertices()) {
        if (!(vInd[v] != vi && vInd[v] != vj && vInd[v] != vk)) {
            for (Corner c: v.adjacentCorners()) {
                rows.emplace_back(count);
                cols.emplace_back(cInd[c]);
                values.emplace_back(1.);
            }
            count++;
        }
    }
    rhs = vector<double>(count, PI);
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto bdryFlat = Matrix::sparse(count, flatmesh->nCorners(), r, c, v);
    auto bdryPI = new_array_ptr(rhs);
    //M->constraint("Flat boundary", Expr::mul(bdryFlat, alpha), Domain::equalsTo(bdryPI));
    rows.clear();
    cols.clear();
    values.clear();



    // Sum to pi in each triangle constraint
    rhs = vector<double>(nFaces, PI);
    for (Face f : flatmesh->faces()) {
        for (Vertex v: f.adjacentVertices()) {
            rows.emplace_back(fInd[f]);
            cols.emplace_back(vInd[v]);
            values.emplace_back(1);
        }
    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto sumConst = Matrix::sparse(nFaces, flatmesh->nCorners(), r, c, v);
    auto sumPI = new_array_ptr(rhs);
    M->constraint("sum constraint", Expr::mul(sumConst, alpha), Domain::equalsTo(sumPI));
    rows.clear();
    cols.clear();
    values.clear();


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
    CirclePatterns prob(flatmesh, 0, sol, eInd, vInd, fInd, thetas);
    cout << "starting parameterization" << endl;
    prob.parameterize();
    cout << "parameterization done" << endl;
    prob.dbgSVG("flat.svg");
}