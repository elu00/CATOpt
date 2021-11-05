#include "CatOpt.h"
#include "Solver.h"
#include "CirclePatterns.h"
#include "fusion.h"
using namespace mosek::fusion;
using namespace monty;

void CatOpt::circlePatterns() {
    //targetAngles = CornerData<double>(*mesh);
    cout << "Circle pattern LP" << endl;
    EdgeData<size_t> eInd = flatmesh->getEdgeIndices();
    VertexData<size_t> vInd = flatmesh->getVertexIndices();
    CornerData<size_t> cInd = flatmesh->getCornerIndices();
    FaceData<size_t> fInd = flatmesh->getFaceIndices();
    VertexData<bool> mark(*flatmesh, false);
    size_t vi, vj, vk;
    /*
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
    */



    Eigen::VectorXd thetas(flatmesh->nEdges());

    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    //Variable::t Theta = M->variable("Theta", flatmesh->nEdges(), Domain::inRange(-PI, PI));
    Variable::t alpha = M->variable("alpha", flatmesh->nCorners(), Domain::inRange(0, PI));
    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    vector<double> rhs(nEdges, 0);

    /*
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

    */

    // Intersection angle constraint
    rhs = vector<double>(nEdges, 0);
    for (Edge e : flatmesh->edges()) {
        Halfedge he = e.halfedge();
        if (!e.isBoundary()) {
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
            rhs[eInd[e]] = flatGeometry->cornerAngle(e.halfedge().next().next().corner()) 
                            + flatGeometry->cornerAngle(e.halfedge().twin().next().next().corner());
        }
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto intersectionConst = Matrix::sparse(nEdges, flatmesh->nCorners(), r, c, v);
    auto intersectionPI = new_array_ptr(rhs);
    M->constraint("Intersection Agreement", Expr::mul(intersectionConst, alpha), Domain::equalsTo(intersectionPI));
    rhs.clear();
    rows.clear();
    cols.clear();
    values.clear();



    // Flat Boundary Constraint
    size_t excl = 3;
    size_t count = 0;
    for (Vertex v: flatmesh->boundaryLoop(0).adjacentVertices()) {
        if (count >= excl) {
        //if (!(vInd[v] != vi && vInd[v] != vj && vInd[v] != vk)) {
            for (Corner c: v.adjacentCorners()) {
                rows.emplace_back(count);
                cols.emplace_back(cInd[c]);
                values.emplace_back(1.);
            }
            rhs.push_back(PI);
        } else {
            rhs.push_back(0);
        }
        count++;
    }
    for (size_t i = 0; i < excl; i++) rhs[i] = 0;
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto bdryFlat = Matrix::sparse(count, flatmesh->nCorners(), r, c, v);
    auto bdryPI = new_array_ptr(rhs);
    M->constraint("Flat boundary", Expr::mul(bdryFlat, alpha), Domain::equalsTo(bdryPI));
    rows.clear();
    cols.clear();
    values.clear();



    // Sum to pi in each triangle constraint
    rhs = vector<double>(nFaces, PI);
    for (Face f : flatmesh->faces()) {
        for (Corner c: f.adjacentCorners()) {
            rows.emplace_back(fInd[f]);
            cols.emplace_back(cInd[c]);
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

    auto avg = std::make_shared<ndarray<double, 1>>(shape(nCorners), PI/3);
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::vstack(t, Expr::sub(alpha,avg)), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << "LP Done" << endl;
    cout << M->getProblemStatus() << endl;
    //auto xsize = Theta->getSize();
    //auto xVal = Theta->level();
    auto asize = alpha->getSize();
    auto asol = alpha->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    /*
    for (int i = 0; i < xsize; ++i) {
        thetas[i] = (*xVal)[i];
    }
    */
   for (Edge e : flatmesh->edges()) {
        double a1 = e.halfedge().isInterior() ? (*asol)[cInd[e.halfedge().next().next().corner()]] : 0;
        double a2 = e.halfedge().twin().isInterior() ? (*asol)[cInd[e.halfedge().twin().next().next().corner()]] : 0;
        thetas[eInd[e]] = PI - a1 - a2;
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
    uv = prob.parameterize();
    cout << "parameterization done" << endl;
    // calculate circle center to invert around
    count = 0;
    invRadius = 0.5;
    Eigen::Vector2d b1;
    Eigen::Vector2d b2;
    for (Vertex v: flatmesh->boundaryLoop(0).adjacentVertices()) {
        if (count == 0) {
            b1 = uv[v];
        } else if (count == excl - 1) {
            b2 = uv[v];
            break;
        }
        count++;
    }
    auto offset = Eigen::Vector2d((b2-b1).y(), -(b2-b1).x());
    center = (b1 + b2)/2. + invRadius * offset/offset.norm();
    circleSol = Vector<double>(flatmesh->nEdges());
    for (int i = 0; i < flatmesh->nEdges(); i++) circleSol[i] = 0;
    uvSVG("flat.svg");
    circleInversion();
    uvSVG("inverted.svg");
}

void CatOpt::uvSVG(std::string filename) {

    std::ofstream ss(filename, std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
       << "<svg width=\"2000\" height=\"2000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
    for (Vertex v : flatmesh->vertices()) {
        Eigen::Vector2d thing = uv[v];
        ss << "<circle cx=\"" << shift(thing.x()) << "\" cy=\"" << shift(thing.y()) << "\" r=\"1\"/>" << endl;
    }
    ss << "<circle cx=\"" << shift(center.x()) << "\" cy=\"" << shift(center.y()) << "\" r=\"2\"/>" << endl;
    ss << "<circle cx=\"" << shift(center.x()) << "\" cy=\"" << shift(center.y()) << "\" r=\"" << 500 * invRadius <<"\" stroke=\"black\" fill-opacity=\"0\"/>" << endl;
    for (Edge e : flatmesh->edges()) {
        
        Eigen::Vector2d i = uv[e.halfedge().vertex()];
        Eigen::Vector2d j = uv[e.halfedge().twin().vertex()];
        double angle = -circleSol[eInd[e]];
        // FROM NORMALIZATION
        double radius = 500 * (i - j).norm() / abs(2 * sin(angle));
        if (abs(angle) > 1e-7) {
            std::string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
            // sweep flag is 1 if going outward, 0 if going inward
            std::string sweepFlag = angle < 0 ? "0" : "1";
            ss << "<path d=\"M" << shift(i.x()) << "," << shift(i.y()) << " A" << radius << ","
            << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
            << shift(j.x()) << "," << shift(j.y()) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
        } else {
            ss << "<line x1=\"" << shift(i.x()) << "\" x2=\"" << shift(j.x()) 
            << "\" y1=\"" << shift(i.y()) <<"\" y2=\"" << shift(j.y()) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
        }
        
    }
  
    // footer
    ss << "</svg>";
    std::cout << "done";
}
void CatOpt::circleInversion() {
    for (Vertex v : flatmesh->vertices()) {
        uv[v] = invRadius * invRadius * (uv[v] -center)/((uv[v] - center).squaredNorm()) + center;
    }
}
void CatOpt::setOffsets() {
    EdgeData<size_t> eInd = flatmesh->getEdgeIndices();
    VertexData<size_t> vInd = flatmesh->getVertexIndices();
    CornerData<size_t> cInd = flatmesh->getCornerIndices();
    FaceData<size_t> fInd = flatmesh->getFaceIndices();
    geometrycentral::SparseMatrix<double> A(flatmesh->nCorners(), flatmesh->nEdges());
    Vector<double> rhs(flatmesh->nCorners());
    vector<Eigen::Triplet<double>> triplets;
    for (Corner C: flatmesh->corners()) {
        Halfedge h = C.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        if (h.isInterior()) {
            const Vector2 &a = flattened[vInd[ab.vertex()]];
            const Vector2 &b = flattened[vInd[bc.vertex()]];
            const Vector2 &c = flattened[vInd[ca.vertex()]];
            auto u = b - a;
            auto v = c - a;
            double tAngle = orientedAngle(u,v);
            // alphas are offsets from canonical halfedge orientation
            if (h.edge().halfedge() == h) {
                tAngle += alphas[eInd[h.edge()]];
                triplets.push_back({(int)cInd[C], (int)eInd[h.edge()],1});
            } else {
                tAngle -= alphas[eInd[h.edge()]];
                triplets.push_back({(int)cInd[C], (int)eInd[h.edge()],-1});
            }
            if (ca.edge().halfedge() == ca) {
                tAngle += alphas[eInd[ca.edge()]];
                triplets.push_back({(int)cInd[C], (int)eInd[ca.edge()],1});
            } else {
                tAngle -= alphas[eInd[ca.edge()]];
                triplets.push_back({(int)cInd[C], (int)eInd[ca.edge()],-1});
            }   
            auto a0 = uv[ab.vertex()];
            auto b0 = uv[bc.vertex()];
            auto c0 = uv[ca.vertex()];
            auto u0 = b0 - a0;
            auto v0 = c0 - a0;
            double actAngle = -orientedAngle(Vector2({u0.x(), u0.y()}),Vector2({v0.x(), v0.y()}));
            cout << "act" << actAngle << endl;
            rhs[cInd[C]] = tAngle - actAngle;
        } else {cout << "WAAAAH";}
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    // Build the solver
    Solver<double> solver(A);

    // Solve a problem


    circleSol = solver.solve(rhs);
    // Some solvers have extra powers.
    // Solver<> can compute matrix rank, since it uses QR under the hood.
    std::cout << "matrix rank is " << solver.rank() << std::endl;
    std::cout << "nEdges is " << flatmesh->nEdges() << std::endl;

    uvSVG("offset.svg");
    for (Edge e: flatmesh->boundaryLoop(0).adjacentEdges()) {
        cout << circleSol[eInd[e]];
    }
}