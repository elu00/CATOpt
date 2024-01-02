#include "IntrinsicFlattening.h"

#include "nasoq/nasoq_eigen.h"
// Code order:
// - Constructor
// - Intrinsic flattening constraints + problem setup and solver
// - CAS solving stuff
// - Utility functions
// ===============================Constructor==================================
// Constructor for IntrinsicFlattening.
// Initializes all relevant geometric data, l is a vector of target edge lengths
IntrinsicFlattening::IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,
                                         EdgeData<double> l)
    : mesh(mesh), cornerAngles(*mesh), nVertices(mesh->nVertices()),
      nEdges(mesh->nEdges()), nCorners(mesh->nCorners()),
      nFaces(mesh->nFaces()), c_(mesh->getCornerIndices()),
      e_(mesh->getEdgeIndices()), v_(mesh->getVertexIndices()),
      f_(mesh->getFaceIndices()) {

    for (Corner c : mesh->corners()) {
        double l_ij = l[c.halfedge().edge()];
        double l_jk = l[c.halfedge().next().edge()];
        double l_ki = l[c.halfedge().next().next().edge()];
        cornerAngles[c] = CornerAngle(l_ij, l_jk, l_ki);
    }
}

//=========================Intrinsic flattening code===========================
// We are interesting in minimizing ||β - origAngles||₂ subject to the
// following constraints:
// - CATValidityConstraint,
// - VertexAngleSumConstraint
// - OffsetConstraints
// To do, we set x = [β; α] and solve a system
// of the form min xᵀ A x + xᵀ b
// Here β represents the angle at any corner, and α is an offset per edge to the
// original angle
// TODO: check
CornerData<double> IntrinsicFlattening::SolveIntrinsicOnly() {

    CornerData<double> originalAngles(*mesh);
    for (Corner c : mesh->corners()) {
        originalAngles[c] = cornerAngles[c];
    }

    // Containers for Eigen triplet type
    vector<Eigen::Triplet<double>> At;
    vector<Eigen::Triplet<double>> Ct;
    vector<Eigen::Triplet<double>> Et;

    // Initialize A, b such that xᵀ A x + xᵀ b = ||β - origAngles||₂
    auto [A0t, b0t] = AngleDeviationPenalty(originalAngles);
    addTriples(At, A0t);
    // diagonal perturbation needed for nasoq
    for (int i = 0; i < nCorners + nEdges; i++) {
        At.push_back(Eigen::Triplet<double>(i, i, 1e-6));
    }
    auto A = constructMatrix(At, nCorners + nEdges, nCorners + nEdges);

    // Initializing matrices C, d
    // 4 * |E| x |C|
    auto [C0t, d0t] = CATValidityConstraint();
    addTriples(Ct, C0t);
    auto C = constructMatrix(Ct, 4 * nCorners, nCorners + nEdges);

    // Initializing matrices E, f
    // |V| x |C|
    auto [E0t, f0t] = VertexAngleSumConstraint(ComputeTargetCurvatures());
    // |C| x |C| + |E|
    auto [E1t, f1t] = OffsetConstraints();
    addTriples(Et, E0t);
    addTriples(Et, E1t, nVertices, 0);
    auto E = constructMatrix(Et, nVertices + nCorners, nCorners + nEdges);

    // construct the eigen types
    Eigen::VectorXd b = concat(nCorners + nEdges, b0t);
    Eigen::VectorXd d = concat(4 * nCorners, d0t);
    Eigen::VectorXd f = concat(nVertices + nCorners, f0t, f1t);

    // call the wrapper.
    Eigen::VectorXd x = QPSolve(A, b, C, d, E, f);

    // read the solution into geometry-central data types
    CornerData<double> beta(*mesh);

    for (Corner c : mesh->corners()) {
        size_t index = c_[c];
        beta[c] = x[index];
        // alpha[c] = x[nCorners + index];
    }
    return beta;
}

//=========================Intrinsic only constraints===========================
SparseSystem
IntrinsicFlattening::AngleDeviationPenalty(CornerData<double> angles) {
    // The energy we are minimizing is
    // ||β - origAngles||² = ||β||² - 2 ⟨β, origAngles⟩ +||origAngles||² ,
    // which corresponds to xᵀ I_|C| x -  xᵀ 2 origAngles + constant
    vector<T> tripletList;
    vector<double> linearTerm = vector<double>(nCorners, 0);
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i, i, 1.});
    }
    for (Corner c : mesh->corners()) {
        linearTerm[c_[c]] = -2 * angles[c];
    }
    return {tripletList, linearTerm};
}
// 4 |C| x |C|
// Returns a matrix A and a vector b such that
// Aβ ≤ b represents the constraint that each CAT is valid
// TODO: implement
bool IntrinsicFlattening::CheckCATValidityConstraint(CornerData<double> beta) {
    return true;
}
// TODO: check
SparseSystem IntrinsicFlattening::CATValidityConstraint() {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(4 * nCorners, 0);

    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i, i, -1});
        tripletList.push_back({nCorners + i, i, 1});
    }
    for (Corner c : mesh->corners()) {
        int i = c_[c];
        int j = c_[c.halfedge().next().corner()];
        int k = c_[c.halfedge().next().next().corner()];
        tripletList.push_back({2 * nCorners + i, i, -1});
        tripletList.push_back({2 * nCorners + i, j, 1});
        tripletList.push_back({2 * nCorners + i, k, 1});

        tripletList.push_back({3 * nCorners + i, i, 1});
        tripletList.push_back({3 * nCorners + i, j, -1});
        tripletList.push_back({3 * nCorners + i, k, -1});
    }
    for (int i = 0; i < nCorners; i++) {
        rhs[i] = 0;
        rhs[1 * nCorners + i] = 2 * PI;
        rhs[2 * nCorners + i] = 3 * PI;
        rhs[3 * nCorners + i] = PI;
    }
    return {tripletList, rhs};
}
// TODO: implement
bool IntrinsicFlattening::CheckVertexAngleSumConstraint(
    // |V| x |C|
    // =========== Equation [3]=====================
    // For each interior vertex, add a constraint that
    // says the angles around the vertex sum to 2π - Ω_i:
    //     _______
    //    /\     /\
    //   /  \   /  \
    //  /   2\1/0   \
    //  ------*------
    //  \   3/…\n   /
    //   \  /   \  /
    //    \/_____\/
    // and similarly for boundary vertices
    CornerData<double> beta, VertexData<double> curvatures) {
    for (Vertex v : mesh->vertices()) {
        if (!v.isBoundary()) {
            double accum = 0.;
            for (Corner c : v.adjacentCorners()) {
                accum += beta[c];
            }
        }
    }
    return true;
}
// TODO: check
SparseSystem
IntrinsicFlattening::VertexAngleSumConstraint(VertexData<double> curvatures) {

    vector<T> tripletList;
    vector<double> rhs = vector<double>(nVertices, 0);

    for (Vertex v : mesh->vertices()) {
        if (!v.isBoundary()) {
            for (Corner c : v.adjacentCorners()) {
                tripletList.push_back({v_[v], c_[c], 1.});
            }
            rhs[v_[v]] = 2 * PI - curvatures[v];
        } else if (curvatures[v] < PI) {
            // only add these constraints if we want to actually constrain the
            // curvature on the boundary
            for (Corner c : v.adjacentCorners()) {
                tripletList.push_back({v_[v], c_[c], 1.});
            }
            rhs[v_[v]] = PI - curvatures[v];
        }
    }
    return {tripletList, rhs};
}
void IntrinsicFlattening::CheckConstraintsIntrinsicOnly(
    CornerData<double> beta) {
    assert(CheckCATValidityConstraint(beta));
    assert(CheckVertexAngleSumConstraint(beta, ComputeTargetCurvatures()));
}
SparseSystem IntrinsicFlattening::OffsetConstraints() {
    // |C| x |C| + |E|
    // Returns a matrix A and a vector b such that
    // A [β; α] = b represents the constraint that β_i = α_ij + α_ki + θ_i
    // TODO: check

    vector<T> tripletList;
    vector<double> rhs = vector<double>(nCorners, 0);

    for (Corner c : mesh->corners()) {
        size_t index = c_[c];
        tripletList.push_back({index, index, 1.});

        tripletList.push_back({index, nCorners + e_[c.halfedge().edge()], -1.});
        tripletList.push_back(
            {index, nCorners + e_[c.halfedge().twin().next().edge()], -1.});

        rhs[index] = cornerAngles[c];
    }
    return {tripletList, rhs};
}
//=========================CAS code===========================
// TODO: check, see where this fits into the pipeline
pair<CornerData<double>, CornerData<double>>
IntrinsicFlattening::CoherentAngleSystem(VertexData<double> targetCurvatures,
                                         CornerData<double> targetBetas) {
    // Initialize A, b
    auto [A0t, bt] = AngleDeviationPenalty(targetBetas);
    vector<Eigen::Triplet<double>> At;
    addTriples(At, A0t, nCorners, nCorners);
    // Required by NASOQ for some reason
    for (int i = 0; i < 2 * nCorners; i++) {
        At.push_back(Eigen::Triplet<double>(i, i, 1e-6));
    }
    vector<double> temp(nCorners, 0);
    auto A = constructMatrix(At, 2 * nCorners, 2 * nCorners);
    Eigen::VectorXd b = concat(2 * nCorners, temp, bt);

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
    auto C =
        constructMatrix(Ct, nCorners + nEdges + 4 * nCorners, 2 * nCorners);
    Eigen::VectorXd d = concat(nCorners + nEdges + 4 * nCorners, d0, d1, d2);

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
    auto E = constructMatrix(Et, nFaces + nVertices + nEdges, 2 * nCorners);
    Eigen::VectorXd f = concat(nFaces + nVertices + nEdges, f0, f1, f2);

    // output vectors
    Eigen::VectorXd x = QPSolve(A, b, C, d, E, f);
    CornerData<double> thetas(*mesh);
    CornerData<double> betas(*mesh);
    for (Corner c : mesh->corners()) {
        size_t index = c_[c];
        thetas[c] = x[index];
        betas[c] = x[nCorners + index];
    }
    return {thetas, betas};
}
//=========================CAS constraints===========================

// |C| x |C|
// Returns a matrix A and a vector b such that
// A θ ≤ b represents the constraint that θ is nonnegative
// TODO: implement
bool IntrinsicFlattening::CheckPositiveAngleConstraint(
    CornerData<double> beta) {
    return true;
}
// TODO: check
SparseSystem IntrinsicFlattening::PositiveAngleConstraint() {
    vector<T> tripletList;
    vector<double> rhs = vector<double>(nCorners, 0);
    for (int i = 0; i < nCorners; i++) {
        tripletList.push_back({i, i, -1.});
    }
    return {tripletList, rhs};
}
// |F| x |C|
// TODO: implement
bool IntrinsicFlattening::CheckFaceAngleSumConstraint(CornerData<double> beta) {
    return true;
}
// TODO: check
SparseSystem IntrinsicFlattening::FaceAngleSumConstraint() {
    // =========== Equation [2]================================
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
            tripletList.push_back({f_[f], c_[c], 1.});
        }
    }

    return {tripletList, rhs};
}

// |E| x |C|
// TODO: implement
bool IntrinsicFlattening::CheckEdgeDelaunayConstraint(CornerData<double> beta) {
    return true;
}
// TODO: check
SparseSystem IntrinsicFlattening::EdgeDelaunayConstraint() {
    // =========== Equation [4]=============================
    // For each interior
    // edge, add a constraint which ensures that the two opposite angles in the
    // CAS satisfy the Delaunay condition , namely
    //
    //   α0 + α1 < π       [4]
    //
    // where α are the angles in the coherent angle system and corners are
    // indexed as below:
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
        if (!e.isBoundary()) {
            size_t eInd = e_[e];
            rhs[eInd] = PI;
            tripletList.push_back(
                {eInd, c_[e.halfedge().next().next().corner()], 1.});
            tripletList.push_back(
                {eInd, c_[e.halfedge().twin().next().next().corner()], 1.});
        }
    }
    return {tripletList, rhs};
}

// TODO: implement
bool IntrinsicFlattening::CheckEdgeIntersectionAngleConstraint(
    CornerData<double> beta) {
    return true;
}
// |E| x 2 |C|
// TODO: check
SparseSystem IntrinsicFlattening::EdgeIntersectionAngleConstraint() {
    // =========== Equation [1]=============================
    // For each interior edge,
    // add a constraint which ensures that the coherent angle system exhibits
    // the same circumcircle intersection angles as the input CAT, namely
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

            tripletList.push_back({ind, nCorners + c_[c0], -0.5});
            tripletList.push_back({ind, nCorners + c_[c1], -0.5});

            tripletList.push_back({ind, nCorners + c_[c2], 0.5});
            tripletList.push_back({ind, nCorners + c_[c3], 0.5});
            tripletList.push_back({ind, nCorners + c_[c4], 0.5});
            tripletList.push_back({ind, nCorners + c_[c5], 0.5});

            // (α0 + α1)
            tripletList.push_back({ind, c_[c0], 1.});
            tripletList.push_back({ind, c_[c1], 1.});
            rhs[ind] = PI;
        }
    }

    return {tripletList, rhs};
}
// TODO: check

VertexData<double> IntrinsicFlattening::ComputeTargetCurvatures() {
    VertexData<double> targetCurvatures(*mesh);
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            // any value > PI corresponds to
            // unconstrained boundary curvature
            targetCurvatures[v] = 2 * PI;
        } else {
            targetCurvatures[v] = 0;
        }
    }
    return targetCurvatures;
}

//==========================Utility functions and constructor ================
// weird edge case that was breaking nasoq earlier; seems to be working now
void IntrinsicFlattening::nasoqTest() {

    vector<Eigen::Triplet<double>> HList = {Eigen::Triplet<double>(0, 0, 1),
                                            Eigen::Triplet<double>(1, 1, 1),
                                            Eigen::Triplet<double>(2, 2, 1)};
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> H(6, 6);
    H.setFromTriplets(HList.begin(), HList.end());

    vector<Eigen::Triplet<double>> AList = {
        Eigen::Triplet<double>(0, 0, 1),  Eigen::Triplet<double>(0, 3, -1),
        Eigen::Triplet<double>(0, 5, -1), Eigen::Triplet<double>(1, 1, 1),
        Eigen::Triplet<double>(1, 4, -1), Eigen::Triplet<double>(1, 3, -1),
        Eigen::Triplet<double>(2, 2, 1),  Eigen::Triplet<double>(2, 5, -1),
        Eigen::Triplet<double>(2, 4, -1)};
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> A(3, 6);
    A.setFromTriplets(AList.begin(), AList.end());

    vector<Eigen::Triplet<double>> CList = {
        Eigen::Triplet<double>(0, 0, 1),   Eigen::Triplet<double>(1, 1, 1),
        Eigen::Triplet<double>(2, 2, 1),   Eigen::Triplet<double>(3, 0, -1),
        Eigen::Triplet<double>(6, 0, 1),   Eigen::Triplet<double>(4, 1, -1),
        Eigen::Triplet<double>(7, 1, 1),   Eigen::Triplet<double>(5, 2, -1),
        Eigen::Triplet<double>(8, 2, 1),   Eigen::Triplet<double>(9, 0, -1),
        Eigen::Triplet<double>(9, 1, 1),   Eigen::Triplet<double>(9, 2, 1),
        Eigen::Triplet<double>(12, 0, 1),  Eigen::Triplet<double>(12, 1, -1),
        Eigen::Triplet<double>(12, 2, -1), Eigen::Triplet<double>(10, 1, -1),
        Eigen::Triplet<double>(10, 2, 1),  Eigen::Triplet<double>(10, 0, 1),
        Eigen::Triplet<double>(13, 1, 1),  Eigen::Triplet<double>(13, 2, -1),
        Eigen::Triplet<double>(13, 0, -1), Eigen::Triplet<double>(11, 2, -1),
        Eigen::Triplet<double>(11, 0, 1),  Eigen::Triplet<double>(11, 1, 1),
        Eigen::Triplet<double>(14, 2, 1),  Eigen::Triplet<double>(14, 0, -1),
        Eigen::Triplet<double>(14, 1, -1)};
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> C(15, 6);
    C.setFromTriplets(CList.begin(), CList.end());

    Eigen::VectorXd q(6);
    q << -1.5708, -0.785398, -0.785398, 0, 0, 0;

    Eigen::VectorXd b(3);
    b << 1.5708, 0.785398, 0.785398;
    Eigen::VectorXd d(15);
    d << -3.14159, -3.14159, -3.14159, 0, 0, 0, 6.28319, 6.28319, 6.28319,
        9.42478, 9.42478, 9.42478, 3.14159, 3.14159, 3.14159;

    nasoq::QPSettings *qs = NULL;
    Eigen::VectorXd x, y, z;
    int status = nasoq::quadprog(H, q, A, b, C, d, x, y, z, qs);
}
// returns x minimizing
//      xᵀ A x + xᵀ b
//  subject to the constraints
//      Cx ≤ d
//      Ex = f
//  x is a n x 1 vector, d is a m x 1 vector, and f is a k x 1 vector
//  A is nxn, C is m x n, and E is k x n
//  TODO: add more diagnostics/checking?
Eigen::VectorXd IntrinsicFlattening::QPSolve(
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> &A,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &b,
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> &C,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &d,
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> &E,
    Eigen::Matrix<double, Eigen::Dynamic, 1> &f) {

#ifdef CATDEBUG
    // check dimensions of input
    assert(A.cols() == A.rows());
    int n = A.cols();
    int m = d.rows();
    int k = f.rows();
    assert(C.rows() == m);
    assert(C.cols() == n);
    assert(E.rows() == k);
    assert(E.cols() == n);

#endif // CATDEBUG
    //
    // Dummy variables for nasoq I/O
    Eigen::VectorXd x, y, z;

    nasoq::QPSettings *qs = new nasoq::QPSettings;
    qs->diag_perturb = pow(10, -9);
    qs->eps = pow(10, -3);
    qs->max_iter = 0;
    qs->nasoq_variant = "fixed";

    cout << "calling nasoq" << endl;
    int status = nasoq::quadprog(A, b, E, f, C, d, x, y, z, qs);
    assert(status == 1 && "nasoq solve failed");
    return x;
}
// returns the angle at corner i^jk given the edge lengths l_ij, l_jk, l_ki
// TODO: check this formula
double IntrinsicFlattening::CornerAngle(double l_ij, double l_jk, double l_ki) {
    if (l_ij > l_jk + l_ki || l_ki > l_ij + l_jk) {
        assert(false && "bad edge lengths");
        return 0.0;
    }
    if (l_jk > l_ki + l_ij) {
        assert(false && "bad edge lengths");
        return M_PI;
    }
    return acos((l_ij * l_ij - l_jk * l_jk + l_ki * l_ki) / (2 * l_ij * l_ki));
}
// Convenience function that typecasts and shifts a vector of std::tuples
// into a vector of eigen tuples
// for each entry (a,b,v) in tuples, the tuple (a+i, b+j, v) is appended to
// triples
void IntrinsicFlattening::addTriples(vector<Eigen::Triplet<double>> &triples,
                                     vector<T> &tuples, int i, int j) {
    for (auto [a, b, v] : tuples) {
        triples.push_back(Eigen::Triplet<double>(a + i, b + j, v));
    }
}

// Convenience function for constructing a sparse m x n matrix from a list of
// tuples
Eigen::SparseMatrix<double, Eigen::ColMajor, int>
IntrinsicFlattening::constructMatrix(vector<Eigen::Triplet<double>> &triples,
                                     int m, int n) {
    Eigen::SparseMatrix<double, Eigen::ColMajor, int> A(m, n);
    A.setFromTriplets(triples.begin(), triples.end());
    return A;
}
// convenience function that concatenates a variable number of
// std::vector<double>s into an Eigen Vector
template <typename... Args>
Eigen::VectorXd IntrinsicFlattening::concat(size_t size, Args &...args) {
    Eigen::VectorXd res = Eigen::VectorXd::Zero(size);
    size_t i = 0;
    for (const auto &V : {args...}) {
        for (auto v : V) {
            res[i] = v;
            i++;
        }
    }
    assert(i <= size);
    return res;
}

/*
// TODO: rewrite/review this code
// Given a CAT in the plane, and an assignment of new boundary curvatures
// (defined inline), returns intersection angles for a conformally equivalent
// CAT, and new CAT corner angles (which are the same as input)
pair<EdgeData<double>, CornerData<double>>

IntrinsicFlattening::solveFromPlane(double interpolationWeight) {
    // get corner angles of input CAT mesh
    CornerData<double> targetBetas(*mesh);
    for (Corner c : mesh->corners()) {
        targetBetas[c] = cornerAngles[c];
    }

    // this next block should be replaced at some point
    size_t count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        count++;
    }
    double bdryCount = (double)count;

    VertexData<double> targetCurvatures(*mesh);
    for (Vertex v : mesh->vertices()) {
        if (v.isBoundary()) {
            // any value > PI corresponds to
            // unconstrained boundary curvature
            // targetCurvatures[v] = 2 * PI;
            double angleSum = 0.;
            for (Corner c : v.adjacentCorners()) {
                angleSum += cornerAngles[c];
            }
            // targetCurvatures[v] = PI - angleSum;
            targetCurvatures[v] = interpolationWeight * (2 * PI / bdryCount) +
                                  (1 - interpolationWeight) * (PI - angleSum);
        } else {
            targetCurvatures[v] = 0;
        }
    }

    auto [CAS, beta] = CoherentAngleSystem(targetCurvatures, targetBetas);

    EdgeData<double> thetaSolve(*mesh, 0);
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior()
                        ? CAS[e.halfedge().next().next().corner()]
                        : 0;
        double a2 = e.halfedge().twin().isInterior()
                        ? CAS[e.halfedge().twin().next().next().corner()]
                        : 0;
        thetaSolve[e] = PI - a1 - a2;
    }
    double avgError = 0.;
    double maxError = 0.;
    for (Corner c : mesh->corners()) {
        double error = abs(beta[c] - cornerAngles[c]);
        maxError = std::max(maxError, error);
        avgError += error;
    }
    avgError /= nCorners;
    cout << "max error: " << maxError << endl;
    cout << "avg error: " << avgError << endl;
    return {thetaSolve, beta};
}
*/
