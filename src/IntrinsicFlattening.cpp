#include "fusion.h"
#include "IntrinsicFlattening.h"


using namespace mosek::fusion;
using namespace monty;

IntrinsicFlattening::IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,shared_ptr<VertexPositionGeometry> geometry):
    mesh(mesh), geometry(geometry) {
    geometry->requireEdgeLengths();
    geometry->requireVertexGaussianCurvatures();
    geometry->requireVertexAngleSums();
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    nFaces = mesh->nFaces();
    c_ = mesh->getCornerIndices();
    e_ = mesh->getEdgeIndices();
    v_ = mesh->getVertexIndices();
    
    f_ = mesh->getFaceIndices();
    /* 
    angleDefects = geometry->vertexGaussianCurvatures;
    */

}


// convenience function for making a sparse matrix
monty::rc_ptr<mosek::fusion::Matrix> IntrinsicFlattening::sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values ) {
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto res =  Matrix::sparse(m,n,r,c,v);
    rows.clear();
    cols.clear();
    values.clear();
    return res;
}

// calls mosek to generate and return optimal betas
SolutionData IntrinsicFlattening::solve() {


    // Model initialization
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    // the fixed offset from the Euclidean angles at each edge
    Variable::t alpha = M->variable("alpha", nEdges, Domain::inRange(-PI, PI));
    // the "total angle" at each corner
    Variable::t beta = M->variable("beta", nCorners, Domain::inRange(0, PI));
    // the witness of our "coherent angle system"
    Variable::t a = M->variable("a", nCorners, Domain::inRange(0, 2 * PI));

    // dummy vectors for building matrices
    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    // Equality constraints: beta has to sum to 2 \pi at interior vertices and 
    // \pi at exterior vertices!
    vector<double> rhs(nVertices);
    cout << "nCorners:" << nCorners << endl;
    for (Vertex v: mesh->vertices()) {
        size_t index = v_[v];
        for (Corner c: v.adjacentCorners()) {
            rows.emplace_back(index);
            cols.emplace_back(c_[c]);
            values.emplace_back(1.);
        }
        if (v.isBoundary()) {
            rhs[index] = (PI);
        } else {
            rhs[index] = (2*PI);
        }
    }
    auto sumEquality = sMatrix(nVertices, nCorners, rows, cols, values);
    M->constraint("angle sum equalities", Expr::mul(sumEquality, beta), Domain::equalsTo(new_array_ptr(rhs)));
    rhs.clear();

    // More equality constraints: each beta has to be "explained by" the corresponding alphas
    for (Corner c: mesh->corners()) {
        rhs.push_back(geometry->cornerAngle(c));

        // the two adjacent edges to the corner
        rows.emplace_back(c_[c]);
        cols.emplace_back(e_[c.halfedge().next().next().edge()]);
        values.emplace_back(1.);

        rows.emplace_back(c_[c]);
        cols.emplace_back(e_[c.halfedge().edge()]);
        values.emplace_back(1.);
    }
    auto adjacency = sMatrix(nCorners, nEdges, rows, cols, values);
    // beta = corner angle + alphas
    M->constraint("beta-alpha equalities", Expr::sub(beta, Expr::mul(adjacency,alpha)), Domain::equalsTo(new_array_ptr(rhs)));
    rhs.clear();
    

    // coherent angle system stuff
    //VertexData<bool> mark(*mesh, false);
    Vertex infVertex;

    // Flat Boundary Constraint
    size_t excl = 4;
    size_t count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        if (count == 1) infVertex = v;
        if (count >= excl) {
            for (Corner c : v.adjacentCorners()) {
                rows.emplace_back(count);
                cols.emplace_back(c_[c]);
                values.emplace_back(1.);
            }
            rhs.push_back(PI);
        } else {
            rhs.push_back(0);
        }
        count++;
    }
    auto bdryFlat = sMatrix(count, nCorners, rows, cols, values);
    auto bdryPI = new_array_ptr(rhs);
    M->constraint("Flat boundary", Expr::mul(bdryFlat, a), Domain::equalsTo(bdryPI));
    rhs.clear();

    FaceData<bool> fMask(*mesh, true);
    VertexData<bool> vNewBdry(*mesh, false);
    // halfedges that are now exterior to the mesh after deleting the vertex at infinity
    HalfedgeData<bool> hOutsideBdry(*mesh, false);
    EdgeData<bool> eMask(*mesh, true);
    EdgeData<bool> eBdry(*mesh, false);
    
    // mark all the boundary edges as the new boundary
    for (Edge e : mesh->boundaryLoop(0).adjacentEdges()) {
        eBdry[e] = true;
    }
    // the edges adjacent to the infinite vertex shouldn't be included in
    // future calculations of the boundary
    /*
    for (Halfedge h: infVertex.incomingHalfedges()) {
        vNewBdry[h.vertex()] = true;
    }
    for (Edge e: infVertex.adjacentEdges()) {
        eMask[e] = false;
        eBdry[e] = false;
    }

    for (Face f: infVertex.adjacentFaces()) {
        fMask[f] = false;
        // extra boundary edges
        for (Edge e: f.adjacentEdges()) {
            if (eMask[e]) eBdry[e] = true;
        }
        for (Halfedge h: f.adjacentHalfedges()) {
            hOutsideBdry[h] = true;
        }
    }
    */
    // Intersection angle constraint: 
    // the betas and the as induce the same intersection angle: that is:
    // ++++ -- on the betas = theta =  PI - sum of opposite alphas
    rhs = vector<double>(nEdges,PI);
    vector<int> aRows;
    vector<int> aCols;
    vector<double> aVals;
    for (Edge e : mesh->edges()) {
        size_t ind = e_[e];
        // we don't care about coherence on the boundary
        if (eMask[e] && !e.isBoundary()) {
            vector<Halfedge> temp = {e.halfedge(), e.halfedge().twin()};
            for (auto he: temp) {
                // shared halfedge
                rows.push_back(ind);
                cols.push_back(c_[he.corner()]);
                values.push_back(0.5);
                he = he.next();
                // other half edge
                rows.push_back(ind);
                cols.push_back(c_[he.corner()]);
                values.push_back(0.5);
                he = he.next();
                // corner corresponding to opposite vertex
                rows.push_back(ind);
                cols.push_back(c_[he.corner()]);
                values.push_back(-0.5);

                if (!hOutsideBdry[he]) {
                    aRows.push_back(ind);
                    aCols.push_back(c_[he.corner()]);
                    aVals.push_back(1.);
                }
            }
        } else rhs[ind] = 0;
    }
    
    auto intersectionConst = sMatrix(nEdges, nCorners, rows, cols, values);
    auto aConst = sMatrix(nEdges, nCorners, aRows, aCols, aVals);
    auto intersectionPI = new_array_ptr(rhs);
    M->constraint("Intersection Agreement", 
            Expr::add(Expr::mul(intersectionConst, beta),
                    Expr::mul(aConst, a)),
                  Domain::equalsTo(intersectionPI));
    rhs.clear();


    // Sum to pi in each triangle constraint
    rhs = vector<double>(nFaces, PI);
    for (Face f : mesh->faces()) {
        if(fMask[f]) {
            for (Corner c : f.adjacentCorners()) {
                rows.emplace_back(f_[f]);
                cols.emplace_back(c_[c]);
                values.emplace_back(1.); 
            }
        } else rhs[f_[f]] = 0.;
    }
    auto sumConst = sMatrix(nFaces, nCorners, rows, cols, values);
    auto sumPI = new_array_ptr(rhs);
    M->constraint("sum constraint", Expr::mul(sumConst, a),
                  Domain::equalsTo(sumPI));
    rhs.clear();

    // Sum to 2pi around each vertex
    rhs = vector<double>(nVertices, 0);
    for (Vertex v : mesh->vertices()) {
        if(!v.isBoundary() && v != infVertex && !vNewBdry[v]) {
            rhs[v_[v]] = 2*PI;
            for (Corner c : v.adjacentCorners()) {
                rows.emplace_back(v_[v]);
                cols.emplace_back(c_[c]);
                values.emplace_back(1.); 
            }
        } 
    }
    auto vSum = sMatrix(nVertices, nCorners, rows, cols, values);
    auto vSumRHS = new_array_ptr(rhs);
    M->constraint("vertex sum constraint", Expr::mul(vSum, a), Domain::equalsTo(vSumRHS));
    rhs.clear();

    // local delaunay constraint
    rhs = vector<double>(nEdges, 0);
    for (Edge e : mesh->edges()) {
        if(eMask[e] && !eBdry[e]) {
            size_t eInd = e_[e];
            rhs[eInd] = PI;
            rows.emplace_back(eInd);
            cols.emplace_back(c_[e.halfedge().next().next().corner()]);
            values.emplace_back(1.); 
            rows.emplace_back(eInd);
            cols.emplace_back(c_[e.halfedge().twin().next().next().corner()]);
            values.emplace_back(1.); 
        } 
    }
    auto delaunay = sMatrix(nEdges, nCorners, rows, cols, values);
    auto delaunayRHS = new_array_ptr(rhs);
    M->constraint("Delaunay", Expr::mul(delaunay, a),Domain::lessThan(delaunayRHS));
    rhs.clear();



    // just setup the objective: \sum alpha^2
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::vstack(t, alpha), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = beta->getSize();
    auto xVal = beta->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    CornerData<double> betaSolve(*mesh);
    for (Corner c: mesh->corners()) {
        betaSolve[c] = ((*xVal)[c_[c]]);
    }

    EdgeData<double> thetaSolve(*mesh,0);
    auto asize = a->getSize();
    auto asol = a->level();
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior() && !hOutsideBdry[e.halfedge()]
                        ? (*asol)[c_[e.halfedge().next().next().corner()]]
                        : 0;
        double a2 =
            e.halfedge().twin().isInterior() && !hOutsideBdry[e.halfedge().twin()]
                ? (*asol)[c_[e.halfedge().twin().next().next().corner()]]
                : 0;
        thetaSolve[e] = PI - a1 - a2;
    }
    return {infVertex, eMask, eBdry, fMask, betaSolve, thetaSolve};
}