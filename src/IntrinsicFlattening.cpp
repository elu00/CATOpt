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
    q << -1.5708,-0.785398,-0.785398,0,0,0;

    Eigen::VectorXd b(3);
    b << 1.5708,0.785398,0.785398;
    Eigen::VectorXd d(15);
    d << -3.14159,-3.14159,-3.14159,0,0,0,6.28319,6.28319,6.28319,9.42478,9.42478,9.42478,3.14159,3.14159,3.14159;

    nasoq::QPSettings *qs = NULL;
    Eigen::VectorXd x,y,z;
    int status = nasoq::quadprog(H,q,A,b,C,d,x,y,z,qs);

    
}
// returns the solution x to
//      min x^T A x + x^T b
//  such that
//      Cx ≤ d
//      Ex = f
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
    /*
    cout << C << endl;
    cout << d << endl;
    cout << E << endl;
    cout << f << endl;
    */


    nasoq::QPSettings *qs = new nasoq::QPSettings;
    qs->diag_perturb=pow(10,-9);
    qs->eps=pow(10,-3);
    qs->max_iter = 0;
    qs->nasoq_variant = "fixed";


    // TODO: do something if this is bad
    cout << "calling nasoq" << endl;
    int status = nasoq::quadprog(A,b,E,f,C,d,x,y,z,qs);
    cout << "status: " << status << endl;
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
// convenience function that concatenates a variable number of
// std::vector<double> into an Eigen Vector
template <typename ...Args>
Eigen::VectorXd concat(size_t size, Args & ... args){
    Eigen::VectorXd res = Eigen::VectorXd::Zero(size);
    size_t i = 0;
    for (const auto &V : {args...}) {
        for (auto v: V) {
            res[i] = v;
            i++;
        }
    }
    assert(i <= size);
    return res;
}


// Convenience function that typecasts and shifts a vector of std::tuples
// into a vector of eigen tuples
void IntrinsicFlattening::addTriples(vector<Eigen::Triplet<double>>& triples, vector<T>& tuples, int i, int j){
    for (auto [a,b,v] : tuples) {
        triples.push_back(Eigen::Triplet<double>(a+i,b+j,v));
    }
}

// constructs a sparse m x n matrix from a list of tuples
Eigen::SparseMatrix<double,Eigen::ColMajor,int> IntrinsicFlattening::constructMatrix(vector<Eigen::Triplet<double>>& triples, int m, int n) {
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A(m,n);
    A.setFromTriplets(triples.begin(), triples.end());
    return A;
}
// |C| x |C|
// Returns a matrix A and a vector b such that
// A θ ≤ b represents the constraint that θ is nonnegative
pair<vector<T>, vector<double>> IntrinsicFlattening::PositiveAngleConstraint () {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(nCorners, 0);
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i,i,-1.});
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
   // angles around the vertex sum to 2π - Ω_i:
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
        } else if (curvatures[v] < PI) {
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
// Returns a matrix A and a vector b such that
// Aβ ≤ b represents the constraint that each CAT is valid


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
    addTriples(At, A0t, nCorners, nCorners);
    // Required by NASOQ for some reason
    for (int i = 0; i < 2*nCorners; i++) {
        At.push_back(Eigen::Triplet<double>(i,i,1e-6));
    }
    vector<double> temp(nCorners, 0);
    auto A = constructMatrix(At,2*nCorners,2*nCorners);
    Eigen::VectorXd b = concat(2*nCorners, temp, bt);

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
    vector<Eigen::Triplet<double>> Ct;
    addTriples(Ct, C0t);
    addTriples(Ct, C1t, nCorners, 0);
    addTriples(Ct, C2t, nCorners + nEdges, nCorners);
    auto C = constructMatrix(Ct, nCorners + nEdges + 4*nCorners,2*nCorners);
    Eigen::VectorXd d = concat(nCorners + nEdges + 4*nCorners, d0, d1, d2);

    // |F| x |C|
    auto [E0t, f0] = FaceAngleSumConstraint();
    // |V| x |C|
    auto [E1t, f1] = VertexAngleSumConstraint(targetCurvatures);
    // |E| x  2 |C|
    auto [E2t, f2] = EdgeIntersectionAngleConstraint();
    // 4 |C| x |C|
    // Initialize E and f
    // E = E0 0
    //     E1 0
    //     E2
    vector<Eigen::Triplet<double>> Et;
    addTriples(Et, E0t);
    addTriples(Et, E1t, nFaces, 0);
    addTriples(Et, E2t, nFaces + nVertices, 0);
    auto E = constructMatrix(Et,nFaces + nVertices + nEdges, 2*nCorners);
    Eigen::VectorXd f = concat(nFaces + nVertices + nEdges, f0, f1, f2);


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
// Returns a matrix A and a vector b such that
// A [β; α] = b represents the constraint that β_i = α_ij + α_ki + θ_i
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
    // x = [β; α]
    CornerData<double> beta(*mesh);
    CornerData<double> originalAngles(*mesh);
    for (Corner c: mesh->corners()) {
        originalAngles[c]=geometry->cornerAngle(c);
    }
 
    // Minimizing ||beta - origAngles||
    // Initialize A, b
    auto [A0t, b0t] = AngleDeviationPenalty(originalAngles);
    vector<Eigen::Triplet<double>> At;
    addTriples(At, A0t);
    for (int i = 0; i < nCorners + nEdges; i++) {
        At.push_back(Eigen::Triplet<double>(i,i,1e-6));
    }
    auto A = constructMatrix(At,nCorners+nEdges,nCorners+nEdges);
    Eigen::VectorXd b = concat(nCorners + nEdges, b0t);


    VertexData<double> targetCurvatures(*mesh);
    for (Vertex v: mesh->vertices()) {
        if (v.isBoundary()) {
            // any value > PI corresponds to
            // unconstrained boundary curvature
            targetCurvatures[v] = 2 * PI;
        } else {
            targetCurvatures[v] = 0;
        }
    }
    // 4 * |E| x |C|
    auto [C0t, d0t] = CATValidityConstraint();
    vector<Eigen::Triplet<double>> Ct;
    addTriples(Ct, C0t);
    auto C = constructMatrix(Ct, 4*nCorners, nCorners+nEdges);

    Eigen::VectorXd d = concat(4*nCorners, d0t);

    // |V| x |C|
    auto [E0t, f0t] = VertexAngleSumConstraint(targetCurvatures);
    // |C| x |C| + |E|
    auto [E1t, f1t] = OffsetConstraints();
    vector<Eigen::Triplet<double>> Et;
    addTriples(Et,E0t);
    addTriples(Et,E1t, nVertices, 0);
    auto E = constructMatrix(Et,nVertices + nCorners,nCorners+nEdges);

    Eigen::VectorXd f = concat(nVertices + nCorners, f0t, f1t);


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
pair<EdgeData<double>, CornerData<double>> IntrinsicFlattening::solveFromPlane(double interpolationWeight) {
    // get corner angles of input CAT mesh
    CornerData<double> targetBetas(*mesh);
    for (Corner c : mesh->corners()) {
        targetBetas[c] = geometry->cornerAngle(c);
    }

    // this next block should be replaced at some point
    size_t count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        count++;
    }
    double bdryCount = (double) count; 


    VertexData<double> targetCurvatures(*mesh);
    for (Vertex v: mesh->vertices()) {
        if (v.isBoundary()) {
            // any value > PI corresponds to
            // unconstrained boundary curvature
            //targetCurvatures[v] = 2 * PI;
            double angleSum = 0.;
            for (Corner c: v.adjacentCorners()) {
                angleSum += geometry->cornerAngle(c);
            }
            //targetCurvatures[v] = PI - angleSum;
            targetCurvatures[v] = interpolationWeight * (2*PI/bdryCount) + (1-interpolationWeight) * (PI - angleSum);
        } else {
            targetCurvatures[v] = 0;
        }
    }

    auto [CAS, beta] = CoherentAngleSystem(targetCurvatures, targetBetas);

    EdgeData<double> thetaSolve(*mesh, 0);
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior() ? CAS[e.halfedge().next().next().corner()] : 0;
        double a2 =
            e.halfedge().twin().isInterior() ? CAS[e.halfedge().twin().next().next().corner()] : 0;
        thetaSolve[e] = PI - a1 - a2;
    }
    double avgError = 0.;
    double maxError = 0.;
    for (Corner c: mesh->corners()) {
        double error = abs(beta[c] - geometry->cornerAngle(c));
        maxError = std::max(maxError, error);
        avgError += error;
    }
    avgError /= nCorners;
    cout << "max error: " << maxError << endl;
    cout << "avg error: " << avgError << endl;
    return {thetaSolve, beta};
}
