#include "IntrinsicFlattening.h"

#include "nasoq/nasoq_eigen.h"
void IntrinsicFlattening::nasoqTest() {
    vector<Eigen::Triplet<double>> HList = {Eigen::Triplet<double>(0,0,1),
        Eigen::Triplet<double>(1,1,1),
        Eigen::Triplet<double>(2,2,1)};
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> H(6,6);
    H.setFromTriplets(HList.begin(), HList.end());

    vector<Eigen::Triplet<double>> AList = {Eigen::Triplet<double>(0,0,1),
        Eigen::Triplet<double>(0,3,-1),
        Eigen::Triplet<double>(0,5,-1),
        Eigen::Triplet<double>(1,1,1),
        Eigen::Triplet<double>(1,4,-1),
        Eigen::Triplet<double>(1,3,-1),
        Eigen::Triplet<double>(2,2,1),
        Eigen::Triplet<double>(2,5,-1),
        Eigen::Triplet<double>(2,4,-1)};
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A(3,6); 
    A.setFromTriplets(AList.begin(), AList.end());

    vector<Eigen::Triplet<double>> CList = {Eigen::Triplet<double>(0,0,1),
        Eigen::Triplet<double>(1,1,1),
        Eigen::Triplet<double>(2,2,1),
        Eigen::Triplet<double>(3,0,-1),
        Eigen::Triplet<double>(6,0,1),
        Eigen::Triplet<double>(4,1,-1),
        Eigen::Triplet<double>(7,1,1),
        Eigen::Triplet<double>(5,2,-1),
        Eigen::Triplet<double>(8,2,1),
        Eigen::Triplet<double>(9,0,-1),
        Eigen::Triplet<double>(9,1,1),
        Eigen::Triplet<double>(9,2,1),
        Eigen::Triplet<double>(12,0,1),
        Eigen::Triplet<double>(12,1,-1),
        Eigen::Triplet<double>(12,2,-1),
        Eigen::Triplet<double>(10,1,-1),
        Eigen::Triplet<double>(10,2,1),
        Eigen::Triplet<double>(10,0,1),
        Eigen::Triplet<double>(13,1,1),
        Eigen::Triplet<double>(13,2,-1),
        Eigen::Triplet<double>(13,0,-1),
        Eigen::Triplet<double>(11,2,-1),
        Eigen::Triplet<double>(11,0,1),
        Eigen::Triplet<double>(11,1,1),
        Eigen::Triplet<double>(14,2,1),
        Eigen::Triplet<double>(14,0,-1),
        Eigen::Triplet<double>(14,1,-1)};
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C(15,6); 
    C.setFromTriplets(CList.begin(), CList.end());

    Eigen::VectorXd q(6);
    q[0] = 0; q[1] = 0; q[2] = 0;
    q[3] = 0; q[4] = 0; q[5] = 0;
    Eigen::VectorXd b(3);
    b[0] = 0; b[1] = 1; b[2] = 0;
    Eigen::VectorXd d(15);
    d[1] = -2;
    for (int i = 0; i < 15; i++) d[i] = 0;


    Eigen::VectorXd x = QPSolve(H,q,C,d,A,b);
    for (int i = 0; i < 2; i++) {
        cout << "wah" << x[i] << endl;
    }
    
}
Eigen::VectorXd IntrinsicFlattening::QPSolve(
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A,
    Eigen::Matrix<double,Eigen::Dynamic,1> b,
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C,
    Eigen::Matrix<double,Eigen::Dynamic,1> d,
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> E,
    Eigen::Matrix<double,Eigen::Dynamic,1> f) {

    Eigen::VectorXd x, y, z;
    // DEBUG
    //cout << A.isApprox(A.triangularView<Eigen::Lower>(),0) << "WASSIM" << endl;
    cout << A.nonZeros() << endl;
    cout << C.nonZeros() << endl;
    cout << E.nonZeros() << "E" << endl;

    nasoq::QPSettings *qs = new nasoq::QPSettings;
    qs->diag_perturb=pow(10,-9);
    qs->eps=pow(10,-3);
    qs->max_iter = 0;
    qs->nasoq_variant = "fixed";


    // TODO: do something if this is bad
    cout << "calling solver..." << endl;
    int status = nasoq::quadprog(A,b,E,f,C,d,x,y,z,qs);
    return x;
}
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


void IntrinsicFlattening::shiftTriples(vector<T>& tripletList, int i, int j) {
    for (auto& [a,b,v]: tripletList) {
        a += i;
        b += j;
    }
}
// |C| x |C|
pair<vector<T>, vector<double>> IntrinsicFlattening::PositiveAngleConstraint () {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(nCorners, 0);
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({-i,-i,1.});
    }
    return {tripletList, rhs};
}
// |F| x |C|
pair<vector<T>, vector<double>> IntrinsicFlattening::FaceAngleSumConstraint () {
    // =========== Equation [2] =================================================

    // For each triangle, add a constraint that says that the
    // interior angles from the coherent angle system sum to π:
    //
    //    α0 + α1 + α2 = π       [2]
    //
    //        *
    //       /0\
    //      /   \
    //     /     \
    //    /1     2\
    //   *---------*
    //

   // arrays for building sparse linear system for Equation [2]
   vector<T> tripletList;
   vector<double> rhs = vector<double>(nFaces, PI); // set all values to π

   for (Face f : mesh->faces()) {
      for (Corner c : f.adjacentCorners()) {
         tripletList.push_back({f_[f],c_[c],1.}); 
      }
   }

    return {tripletList, rhs};
}
// |V| x |C|
pair<vector<T>, vector<double>> IntrinsicFlattening::VertexAngleSumConstraint(VertexData<double> curvatures) {
   // =========== Equation [3] =================================================
   // For each interior vertex, add a constraint that says the
   // angles around the vertex sum to 2π:
   //
   //     _______
   //    /\     /\
   //   /  \   /  \
   //  /   2\1/0   \
   //  ------*------
   //  \   3/…\n   /
   //   \  /   \  /
   //    \/_____\/
   // and similarly for boundary vertices

    vector<T> tripletList;
    vector<double> rhs = vector<double>(nVertices, 0);

    for (Vertex v: mesh->vertices()) {
        if (!v.isBoundary()) {
            for (Corner c: v.adjacentCorners()) {
                tripletList.push_back({v_[v],c_[c],1.});
            }
            rhs[v_[v]] = 2 * PI - curvatures[v];
        } else if (curvatures[v] > PI) {
            // only add these constraints if we want to actually constrain the curvature on the boundary
            for (Corner c: v.adjacentCorners()) {
                tripletList.push_back({v_[v],c_[c],1.});
            }
            rhs[v_[v]] = PI - curvatures[v];
        }
    }
    return {tripletList, rhs};
}
// |E| x |C|
pair<vector<T>, vector<double>> IntrinsicFlattening::EdgeDelaunayConstraint(){
   // =========== Equation [4] =================================================
   // For each interior edge, add a constraint which ensures that
   // the two opposite angles in the CAS satisfy the Delaunay condition , namely
   //
   //   α0 + α1 < π       [4]
   //
   // where α are the angles in the coherent angle system and corners are indexed as below:
   //
   //        *
   //       /0\
   //      /   \
   //     /2   3\
   //    *-------*
   //     \5   4/
   //      \   /
   //       \1/
   //        *
   //

    vector<T> tripletList;
    vector<double> rhs = vector<double>(nEdges, 0);

    // local delaunay constraint
    for (Edge e : mesh->edges()) {
        if(!e.isBoundary()) {
            size_t eInd = e_[e];
            rhs[eInd] = PI;
            tripletList.push_back({eInd,c_[e.halfedge().next().next().corner()],1.}); 
            tripletList.push_back({eInd,c_[e.halfedge().twin().next().next().corner()],1.}); 
        } 
    }
    return {tripletList, rhs};
}

// |E| x 2 |C|
pair<vector<T>, vector<double>>IntrinsicFlattening::EdgeIntersectionAngleConstraint() {
   // =========== Equation [1] =================================================
   // For each interior edge, add a constraint which ensures that
   // the coherent angle system exhibits the same circumcircle
   // intersection angles as the input CAT, namely
   //
   //   α0 + α1 = (β0 + β1 - β2 - β3 - β4 - β5)/2 + π       [1]
   //
   // where α are the angles in the coherent angle system, and β are the
   // CAT corner angles, and corners are indexed as below:
   //
   //        *
   //       /0\
   //      /   \
   //     /2   3\
   //    *-------*
   //     \5   4/
   //      \   /
   //       \1/
   //        *
   //

    vector<T> tripletList;
    vector<double> rhs = vector<double>(nEdges, 0);



   for (Edge e : mesh->edges()) {
      size_t ind = e_[e];

      // only preserve intersection angles on interior edges
      if (!e.isBoundary()) {

         // right-hand side ((β0 + β1 - β2 - β3 - β4 - β5)/2 + π)
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


         tripletList.push_back({ind, nCorners + c_[c0] , -0.5});
         tripletList.push_back({ind, nCorners + c_[c1] , -0.5});

         tripletList.push_back({ind, nCorners + c_[c2] , 0.5});
         tripletList.push_back({ind, nCorners + c_[c3] , 0.5});
         tripletList.push_back({ind, nCorners + c_[c4] , 0.5});
         tripletList.push_back({ind, nCorners + c_[c5] , 0.5});

         // (α0 + α1)
         tripletList.push_back({ind, c_[c0] , 1.});
         tripletList.push_back({ind, c_[c1] , 1.});
         rhs[ind] = PI;
      }
   }

   return {tripletList, rhs};
}
// 4 |C| x |C|
pair<vector<T>, vector<double>>IntrinsicFlattening::CATValidityConstraint() {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(4*nCorners, 0);
    
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i, i, -1});
        tripletList.push_back({nCorners + i, i, 1});
    }
    for (Corner c: mesh->corners()) {
        int i = c_[c]; 
        int j = c_[c.halfedge().next().corner()]; 
        int k = c_[c.halfedge().next().next().corner()]; 
        tripletList.push_back({2*nCorners + i,i,-1});
        tripletList.push_back({2*nCorners + i,j,1});
        tripletList.push_back({2*nCorners + i,k,1});

        tripletList.push_back({3*nCorners + i,i,1});
        tripletList.push_back({3*nCorners + i,j,-1});
        tripletList.push_back({3*nCorners + i,k,-1});
    }
    for (int i = 0; i < nCorners; i++) {
        rhs[i] = 0;
        rhs[nCorners + i] = 2 * PI;
        rhs[2*nCorners + i] = 3 * PI;
        rhs[3*nCorners + i] = PI;
    }
    return {tripletList, rhs};
}
//|C| x |C|
pair<vector<T>, vector<double>>IntrinsicFlattening::AngleDeviationPenalty(CornerData<double> beta) {
    vector<T> tripletList;
    vector<double> linearTerm = vector<double>(nCorners, 0);
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i,i,1.});
    }
    for (Corner c: mesh->corners()) {
        linearTerm[c_[c]] = -beta[c];
    }
    return {tripletList, linearTerm};
}

pair<CornerData<double>, CornerData<double>> IntrinsicFlattening::CoherentAngleSystem(VertexData<double> targetCurvatures, CornerData<double> targetBetas) {
    // Initialize A, b
    auto [A0t, bt] = AngleDeviationPenalty(targetBetas);
    vector<Eigen::Triplet<double>> At;
    for (auto [i,j,v] : A0t) {
        At.push_back(Eigen::Triplet<double>(i,j,v));
    }
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A(2*nCorners,2*nCorners);
    A.setFromTriplets(At.begin(), At.end());
    Eigen::VectorXd b(2*nCorners);
    for (int i = 0; i < nCorners; i++) {
        b[i] = bt[i];
        b[nCorners + i] = 0;
    }

    // |C| x |C|
    auto [C0t, d0] = PositiveAngleConstraint();
    // |E| x |C|
    auto [C1t, d1] = EdgeDelaunayConstraint();
    // 4 |C| x |C|
    auto [C2t, d2] = CATValidityConstraint();
    // Initialize C and d
    // C = C0 0
    //     C1 0
    //     0  C2
    shiftTriples(C1t, nCorners, 0);
    shiftTriples(C2t, nCorners + nEdges, nCorners);
    vector<Eigen::Triplet<double>> Ct;
    for (auto [i,j,v]: C0t) Ct.push_back(Eigen::Triplet<double>(i,j,v));
    for (auto [i,j,v]: C1t) Ct.push_back(Eigen::Triplet<double>(i,j,v));
    for (auto [i,j,v]: C2t) Ct.push_back(Eigen::Triplet<double>(i,j,v));
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C(nCorners + nEdges + 4*nCorners,2*nCorners);
    C.setFromTriplets(Ct.begin(), Ct.end());
    Eigen::VectorXd d(nCorners + nEdges + 4*nCorners);
    for (int i = 0; i < nCorners; i++) {
        d[i] = d0[i];
    }
    for (int i = 0; i < nEdges; i++) {
        d[nCorners + i] = d1[i];
    }
    for (int i = 0; i < 4*nCorners; i++) {
        d[nCorners + nEdges + i] = d2[i];
    }

    // |F| x |C|
    auto [E0t, f0] = FaceAngleSumConstraint();
    // |V| x |C|
    auto [E1t, f1] = VertexAngleSumConstraint(targetCurvatures);
    // |E| x  2 |C|
    auto [E2t, f2] = EdgeIntersectionAngleConstraint();
    // Initialize E and f
    // E = E0 0
    //     E1 0
    //     E2
    shiftTriples(E1t, nFaces, 0);
    shiftTriples(E2t, nFaces + nVertices, 0);
    vector<Eigen::Triplet<double>> Et;
    for (auto [i,j,v]: E0t) Et.push_back(Eigen::Triplet<double>(i,j,v));
    for (auto [i,j,v]: E1t) Et.push_back(Eigen::Triplet<double>(i,j,v));
    for (auto [i,j,v]: E2t) Et.push_back(Eigen::Triplet<double>(i,j,v));
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> E(nFaces + nVertices + nEdges, 2*nCorners);
    E.setFromTriplets(Et.begin(), Et.end());
    Eigen::VectorXd f(nFaces + nVertices + nEdges);
    for (int i = 0; i < nFaces; i++) {
        f[i] = f0[i];
    }
    for (int i = 0; i < nVertices; i++) {
        f[nFaces + i] = f1[i];
    }
    for (int i = 0; i < nEdges; i++) {
        f[nFaces + nVertices + i] = f2[i];
    }


    // output vectors
    Eigen::VectorXd x = QPSolve(A,b,C,d,E,f);
    CornerData<double> thetas(*mesh);
    CornerData<double> betas(*mesh);
    for (Corner c: mesh->corners()) {
        size_t index = c_[c];
        thetas[c] = x[index];
        betas[c] = x[nCorners + index];
    }
    return {thetas, betas};
}
// |C| x |C| + |E|
pair<vector<T>, vector<double>>IntrinsicFlattening::OffsetConstraints() {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(nCorners, 0);

    for (Corner c: mesh->corners()) {
        size_t index = c_[c];
        tripletList.push_back({index, index,1.});

        tripletList.push_back({index, nCorners + e_[c.halfedge().edge()], -1.});
        tripletList.push_back({index, nCorners + e_[c.halfedge().twin().next().edge()], -1.});

        rhs[index]=geometry->cornerAngle(c);
    }
    return {tripletList, rhs};
}


CornerData<double> IntrinsicFlattening::solveIntrinsicOnly() {
    // x = beta
    //     alpha
    // 
    CornerData<double> beta(*mesh);
    // DEBUG
    nasoqTest();
    return beta;
    CornerData<double> originalAngles(*mesh);
    for (Corner c: mesh->corners()) {
        originalAngles[c]=geometry->cornerAngle(c);
    }
 
    // Minimizing ||beta - origAngles||
    // Initialize A, b
    auto [A0t, bt] = AngleDeviationPenalty(originalAngles);
    vector<Eigen::Triplet<double>> At;
    cout << "A" << endl;
    for (auto [i,j,v] : A0t) {
        At.push_back(Eigen::Triplet<double>(i,j,v));
        cout << "Eigen::Triplet<double>(" << i << "," << j << "," << v << ")," << endl;
    }
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A(nCorners+nEdges,nCorners+nEdges);
    A.setFromTriplets(At.begin(), At.end());
    Eigen::VectorXd b(nCorners+nEdges);
    for (int i = 0; i < nCorners; i++) {
        b[i] = bt[i];
    }
    for (int i = 0; i < nEdges; i++) {
        b[nCorners + i] = 0;
    }



    VertexData<double> targetCurvatures(*mesh);
    for (Vertex v: mesh->vertices()) {
        if (v.isBoundary()) {
            // unconstrained boundary curvature
            targetCurvatures[v] = 2 * PI;
        } else {
            targetCurvatures[v] = 0;
        }
    }
    // |V| x |C|
    auto [C0t, d0t] = VertexAngleSumConstraint(targetCurvatures);
    // 4 * |E| x |C|
    auto [C1t, d1t] = CATValidityConstraint();
    shiftTriples(C1t, nVertices, 0);

    cout << "C (inequality)" << endl;
    vector<Eigen::Triplet<double>> Ct;
    for (auto [i,j,v] : C0t) {
        Ct.push_back(Eigen::Triplet<double>(i,j,v));
        // DEBUG
        cout << "Eigen::Triplet<double>(" << i << "," << j << "," << v << ")," << endl;
    }
    for (auto [i,j,v] : C1t) {
        Ct.push_back(Eigen::Triplet<double>(i,j,v));
        cout << "Eigen::Triplet<double>(" << i << "," << j << "," << v << ")," << endl;
    }
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C(nVertices + 4*nCorners, nCorners+nEdges);
    C.setFromTriplets(Ct.begin(), Ct.end());

    Eigen::VectorXd d(nVertices + 4*nCorners);
    
    for (int i = 0; i < nVertices; i++) {
        d[i] = d0t[i];
    }
    for (int i = 0; i < 4*nCorners; i++) {
        d[nVertices + i] = d1t[i];
    }

    cout << "E (equality)" << endl;
    auto [E0t, ft] = OffsetConstraints();
    vector<Eigen::Triplet<double>> Et;
    // DEBUG
    for (auto [i,j,v] : E0t) {
        Et.push_back(Eigen::Triplet<double>(i,j,v));
        cout << "Eigen::Triplet<double>(" << i << "," << j << "," << v << ")," << endl;
    }
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> E(nCorners,nCorners+nEdges);
    E.setFromTriplets(Et.begin(), Et.end());
    Eigen::VectorXd f(nCorners);
    for (int i = 0; i < nCorners; i++) {
        f[i] = ft[i];
    }



    // call the wrapper.
    Eigen::VectorXd x = QPSolve(A,b,C,d,E,f);
    for (Corner c: mesh->corners()) {
        size_t index = c_[c];
        beta[c] = x[index];
    }
    for (Vertex v: mesh->vertices()) {
        if (!v.isBoundary()) {
            double accum = 0.;
            for (Corner c: v.adjacentCorners()) {
                accum += beta[c];
            }
            cout << "angle sum: " << accum << endl;
        }
    }
    return beta;
}

// Given a CAT in the plane, and an assignment of new boundary curvatures (defined inline),
// returns intersection angles for a conformally equivalent CAT, and new CAT corner angles (which are the same as input)
pair<CornerData<double>, CornerData<double>> IntrinsicFlattening::solveFromPlane(double interpolationWeight) {
    // get corner angles of input CAT mesh
    CornerData<double> beta(*mesh);
    for (Corner c : mesh->corners()) {
        beta[c] = geometry->cornerAngle(c);
    }

    /*
    // Initialize MOSEK solver
    Model::t M = new Model();
    auto _M = finally([&]()
                      { M->dispose(); });

    // The values a at corners are the "angles" of a coherent angle system
    // compatible with the prescribed intersection angles (including new boundary angles)
    Variable::t a = M->variable("a", nCorners, Domain::inRange(0, 2 * PI));

    // intersection angles from a = intersection angles from beta
    buildIntersectionAngleConstraints(M, beta, a);
    // CAS sums to pi within each triangle
    buildFaceConstraints(M, a);
    // Vertex sums to 2 pi for each interior vertex
    buildVertexConstraints(M, a);
    // Delaunay constraint
    // Since the input mesh is assumed to be Delaunay this shouldn't actually be necessary?
    buildDelaunayConstraints(M, a);

    // Dummy variable for the objective value
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    buildBoundaryObjective(M, a, t, 2, interpolationWeight);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << "Optimization Done" << endl;
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;

    EdgeData<double> thetaSolve(*mesh, 0);
    auto asize = a->getSize();
    auto asol = a->level();
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior() && e.halfedge().isInterior()
                        ? (*asol)[c_[e.halfedge().next().next().corner()]]
                        : 0;
        double a2 =
            e.halfedge().twin().isInterior() && e.halfedge().twin().isInterior()
                ? (*asol)[c_[e.halfedge().twin().next().next().corner()]]
                : 0;
        thetaSolve[e] = PI - a1 - a2;
    }

    */
    FaceData<bool> fMask(*mesh, true);
    EdgeData<bool> eMask(*mesh, true);
    EdgeData<bool> eBdry(*mesh, false);
    for (Edge e : mesh->edges())
    {
        eBdry[e] = e.isBoundary();
    }
    // DEBUG: TEMPORARY. REPLACE THIS
    EdgeData<double> thetaSolve(*mesh, 0);
    return {beta, beta};
}
