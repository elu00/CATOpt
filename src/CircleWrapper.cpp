#include "CircleWrapper.h"
#include "CirclePatterns.h"
#include "Solver.h"
#include <cmath>



inline double shift(double c) {
    return (c + 2.) * 500;
}
CircleWrapper::CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, CornerData<double> betas, polyscope::SurfaceMesh *psMesh): 
    mesh(mesh), psMesh(psMesh) {
        beta = betas;
        thetas = Eigen::VectorXd::Zero(mesh->nEdges());
        EdgeData<size_t> e_ = mesh->getEdgeIndices();
        for (Edge e: mesh->edges()) {
            size_t ind = e_[e];
            Halfedge h23 = e.halfedge();
            Halfedge h30 = h23.next();
            Halfedge h02 = h30.next();
            Halfedge h45 = h23.twin();
            Halfedge h51 = h45.next();
            Halfedge h14 = h51.next();
            Corner c0 = h02.corner();
            Corner c1 = h14.corner();
            Corner c2 = h23.corner();
            Corner c3 = h30.corner();
            Corner c4 = h45.corner();
            Corner c5 = h51.corner();
            thetas[ind] += beta[c0]/2.;
            thetas[ind] += beta[c1]/2.;
            thetas[ind] -= beta[c2]/2.;
            thetas[ind] -= beta[c3]/2.;
            thetas[ind] -= beta[c4]/2.;
            thetas[ind] -= beta[c5]/2.;
        }
    }
CircleWrapper::CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> intersectionAngles, polyscope::SurfaceMesh *psMesh): 
    mesh(mesh), psMesh(psMesh) {
        theta = intersectionAngles;
        thetas = Eigen::VectorXd::Zero(mesh->nEdges());
        EdgeData<size_t> e_ = mesh->getEdgeIndices();
        for (Edge e: mesh->edges()) {
            thetas[e_[e]] = theta[e];
        }
    }

void CircleWrapper::solve(std::string name) {
    //TODO: change this
    std::cout << "Circle pattern stuff" << endl;
    CirclePatterns prob(mesh, 0, thetas);
    std::cout << "starting parameterization" << std::endl;
    uv = prob.parameterize();
    std::cout << "parameterization done" << std::endl;
    circleSol = Vector<double>(mesh->nEdges());
    for (int i = 0; i < mesh->nEdges(); i++) circleSol[i] = 0;
    uvSVG("flat" + name + ".svg" );
    setOffsets();
    uvSVG(name+".svg");
}

void CircleWrapper::uvSVG(std::string filename) {
    EdgeData<size_t> e_ = mesh->getEdgeIndices();
    std::ofstream ss(filename, std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
        << endl
        << "<svg width=\"2000\" height=\"2000\" "
        "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >"
        << "<rect width=\"100%\" height =\"100%\" fill=\"#ffffff\"/>"
        << endl;
    // point labels
    for (Vertex v : mesh->vertices()) {
        Eigen::Vector2d thing = uv[v];
        ss << "<circle cx=\"" << shift(thing.x()) << "\" cy=\""
            << shift(thing.y()) << "\" r=\"1\"/>" << endl;
    }

    // circle center thing
    /*
       ss << "<circle cx=\"" << shift(center.x()) << "\" cy=\""
       << shift(center.y()) << "\" r=\"2\"/>" << endl;
       ss << "<circle cx=\"" << shift(center.x()) << "\" cy=\""
       << shift(center.y()) << "\" r=\"" << 500 * invRadius
       << "\" stroke=\"black\" fill-opacity=\"0\"/>" << endl;
       */
    for (Edge e : mesh->edges()) {
        Eigen::Vector2d i = uv[e.halfedge().vertex()];
        Eigen::Vector2d j = uv[e.halfedge().twin().vertex()];
        double angle = -circleSol[e_[e]];
        // FROM NORMALIZATION
        double radius = 500 * (i - j).norm() / abs(2 * sin(angle));
        if (abs(angle) > 1e-7) {
            std::string largeArcFlag =
                std::abs(2 * angle) <= 3.14159265358979323846264 ? "0"
                : "1";
            // sweep flag is 1 if going outward, 0 if going inward
            // debug: maybe circle inversion stuff changes this?
            std::string sweepFlag = angle >= 0 ? "0" : "1";
            //std::string sweepFlag = angle < 0 ? "0" : "1";
            ss << "<path d=\"M" << shift(i.x()) << "," << shift(i.y())
                << " A" << radius << "," << radius << " 0 " << largeArcFlag
                << " " << sweepFlag << " " << shift(j.x()) << ","
                << shift(j.y())
                << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />"
                << endl;
        } else {
            ss << "<line x1=\"" << shift(i.x()) << "\" x2=\""
                << shift(j.x()) << "\" y1=\"" << shift(i.y()) << "\" y2=\""
                << shift(j.y()) << "\" stroke=\"blue\" stroke-width=\"2\"/>"
                << endl;
        }
    }
    size_t count = 0;
    size_t excl = 1;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        Eigen::Vector2d thing = uv[v];
        if (count < excl) {
            ss << "<circle fill=\"red\" cx=\"" << shift(thing.x()) << "\" cy=\""
                << shift(thing.y()) << "\" r=\"5\"/>" << endl;
        } 
        count++;
    }
    // footer
    ss << "</svg>";
    std::cout << "done" << endl;
}
// Builds solver for alphas, storing them in circleSol
void CircleWrapper::setOffsets() {
    EdgeData<size_t> e_ = mesh->getEdgeIndices();
    VertexData<size_t> v_ = mesh->getVertexIndices();
    CornerData<size_t> c_ = mesh->getCornerIndices();
    FaceData<size_t> f_ = mesh->getFaceIndices();
    geometrycentral::SparseMatrix<double> A(mesh->nCorners(),
            mesh->nEdges());


    // DEBUG: remove this at some point
    VertexData<bool> bad(*mesh);
    Vector<double> rhs(mesh->nCorners());
    vector<Eigen::Triplet<double>> triplets;
    for (Corner c : mesh->corners()) {
        Halfedge h = c.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        // if (h.isInterior()) {
        // alphas are offsets from canonical halfedge orientation
        if (h.edge().halfedge() == h) {
            triplets.push_back({(int)c_[c], (int)e_[h.edge()], 1});
        } else {
            triplets.push_back({(int)c_[c], (int)e_[h.edge()], -1});
        }
        if (ca.edge().halfedge() == ca) {
            triplets.push_back({(int)c_[c], (int)e_[ca.edge()], 1});
        } else {
            triplets.push_back({(int)c_[c], (int)e_[ca.edge()], -1});
        }
        auto a0 = uv[ab.vertex()];
        auto b0 = uv[bc.vertex()];
        auto c0 = uv[ca.vertex()];
        auto u0 = b0 - a0;
        auto v0 = c0 - a0;
        double actAngle = orientedAngle(Vector2({u0.x(), u0.y()}),
                Vector2({v0.x(), v0.y()}));
        // circle inversion changes orientation
        /*
           double actAngle = -orientedAngle(Vector2({u0.x(), u0.y()}),
           Vector2({v0.x(), v0.y()}));
           */
        //cout << "act" << actAngle << endl;
        rhs[c_[c]] = beta[c] - actAngle;
        if(!std::isfinite(rhs[c_[c]])) {
            cout << "a index" << v_[ab.vertex()] << endl;
            cout << "b index" << v_[bc.vertex()] << endl;
            cout << "c index" << v_[ca.vertex()] << endl;
            cout << "a0:  " << a0.x() << "  " << a0.y() << endl;
            cout << "b0:  " << b0.x() << "  " << b0.y() << endl;
            cout << "c0:  " << c0.x() << "  " << c0.y() << endl;
            cout << "u0:  " << u0.x() << "  " << u0.y() << endl;
            cout << "v0:  " << v0.x() << "  " << v0.y() << endl;
            bad[ab.vertex()] = true;
            bad[bc.vertex()] = true;
            bad[ca.vertex()] = true;
        }
    }
    A.setFromTriplets(triplets.begin(), triplets.end());

    //DEBUG
    //psMesh->addVertexScalarQuantity("bad vertices", bad);
    // Build the solver
    Solver<double> solver(A);
    circleSol = solver.solve(rhs);
    std::cout << "matrix rank is " << solver.rank() << std::endl;
    std::cout << "nEdges is " << mesh->nEdges() << std::endl;
    cout << "residual is" << residual(A,circleSol,rhs) << endl;

    //uvSVG("offset.svg", EdgeData<bool>(*mesh, true));
    /*
       for (Edge e : mesh->boundaryLoop(0).adjacentEdges()) {
       std::cout << circleSol[e_[e]];
       }
       */
    // Double check intersection angle condition
    for (Edge e: mesh->edges()){
        if (e.isBoundary()) continue;
        vector<Halfedge> temp = {e.halfedge(), e.halfedge().twin()};
        double actual = PI;
        double expected = 0;
        for (auto h: temp) {
            Halfedge ab = h;
            Halfedge bc = h.next();
            Halfedge ca = h.next().next();
            auto a0 = uv[ab.vertex()];
            auto b0 = uv[bc.vertex()];
            auto c0 = uv[ca.vertex()];
            auto u0 = a0 - c0;
            auto v0 = b0 - c0;
            actual -= orientedAngle(Vector2({u0.x(), u0.y()}),
                    Vector2({v0.x(), v0.y()}));
            expected += beta[h.corner()] + beta[h.next().corner()] - beta[h.next().next().corner()];
        }
        expected /= 2;
        //cout << "actual intersection angle: " << actual << " expected: " << expected << endl;
        if (abs(actual-expected) > 1e-3) cout << "diff: " << actual - expected << endl;
    }

    // check boundary conditions of adding to PI
    VertexData<double> boundarySum(*mesh);
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        double accum = 0;
        double betaSum = 0;
        for (Corner c : v.adjacentCorners()) {
            double cornerBeta = 0;
            betaSum += beta[c];
            Halfedge h = c.halfedge();
            Halfedge ab = h;
            Halfedge bc = h.next();
            Halfedge ca = h.next().next();
            auto a0 = uv[ab.vertex()];
            auto b0 = uv[bc.vertex()];
            auto c0 = uv[ca.vertex()];
            auto u0 = b0 - a0;
            auto v0 = c0 - a0;
            // circle inversion changes orientation
            cornerBeta += orientedAngle(Vector2({u0.x(), u0.y()}),
                    Vector2({v0.x(), v0.y()}));
            /*
               cornerBeta += -orientedAngle(Vector2({u0.x(), u0.y()}),
               Vector2({v0.x(), v0.y()}));
               */
            if (h.edge().halfedge() == h) {
                cornerBeta += circleSol[e_[h.edge()]];
            } else {
                cornerBeta -= circleSol[e_[h.edge()]];
            }
            if (ca.edge().halfedge() == ca) {
                cornerBeta += circleSol[e_[ca.edge()]];
            } else {
                cornerBeta -= circleSol[e_[ca.edge()]];
            }
            //cout << cornerBeta << "expected " << beta[c] << endl;
            accum += cornerBeta;
        }
        /*
           for (Edge e: v.adjacentEdges()){
           if (e.isBoundary()) {
           if (e.halfedge().isInterior()) {
           accum += circleSol[e_[e]];
           } else {
           accum -= circleSol[e_[e]];
           }
           }

           }
           */
        boundarySum[v] = accum;
    }
    //psMesh->addVertexScalarQuantity("boundary sum", boundarySum);
    }