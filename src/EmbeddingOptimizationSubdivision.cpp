#include "EmbeddingOptimization.h"

#include "polyscope/surface_mesh.h"
// This file contains the implementation of the subdivision code for the
// Embedding optimization. Given a parameter n, this code performs the following
// operations:
// 1. Subdivides each triangular face of the mesh into 3 (n x n) quad grids,
// with each grid corresponding to a corner, which is done in
// buildEquivalenceClasses() and buildIntrinsicCheckerboard()
// 2. Calculates the intrinsic lengths of each quad, using the CAT angles β
// passed in as an argument.
// 3. Returns a new geometry-central mesh that contains the subdivided mesh
//
//
// Constructor for EmbeddingOptimization
EmbeddingOptimization::EmbeddingOptimization(
    shared_ptr<ManifoldSurfaceMesh> mesh,
    shared_ptr<VertexPositionGeometry> geometry, CornerData<double> beta)
    : mesh(mesh), geometry(geometry), LMInitialized(false),
      beta(beta), // the target CAT corner angles
      intrinsicQuantitiesInitialized(
          false), // whether the checkboard has been initialized or not
      nVertices(mesh->nVertices()), // initialized from data
      nEdges(mesh->nEdges()), nCorners(mesh->nCorners()),
      nFaces(mesh->nFaces()),
      v_(mesh->getVertexIndices()), // geometry-central containers for indexing
      e_(mesh->getEdgeIndices()), c_(mesh->getCornerIndices()),
      f_(mesh->getFaceIndices()) {
    // Initialize quantities
    fairnessNormalization = 0.5;
    geometry->requireEdgeLengths();
    geometry->requireVertexGaussianCurvatures();
    geometry->requireVertexAngleSums();
}

// ======================= Subdivision code ====================================
std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry>>
EmbeddingOptimization::initializeSubdivision(int N) {
    n = N;

    // Initialize union find data structure
    top = vector<int>(n * n * nCorners);
    for (int i = 0; i < n * n * nCorners; i++)
        top[i] = i;
    next = vector<int>(n * n * nCorners, -1);

    buildEquivalenceClasses();
    cout << "built equivalence classes" << endl;
    nSubdividedVertices = buildFinalIndices();
    cout << "built final indices" << endl;

    // Initialize submesh and subgeometry
    buildSubdivision();
    cout << "built subdivision" << endl;

    // initialize inital guess for x
    initialSolution = Eigen::VectorXd::Zero(3 * nSubdividedVertices);
    for (Corner c : mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int X = 0; X < n; X++) {
            for (int Y = 0; Y < n; Y++) {
                Vector3 pos = bary(c, X, Y);
                size_t startIndex = 3 * finalIndices[cOffset + X + n * Y];
                initialSolution[startIndex] = pos.x;
                initialSolution[startIndex + 1] = pos.y;
                initialSolution[startIndex + 2] = pos.z;
            }
        }
    }
    buildIntrinsicCheckerboard();
    cout << "built intrinsic lengths" << endl;

    intrinsicQuantitiesInitialized = true;

    // polyscope::show();

    return {submesh, subgeometry};
}
// TODO: double check indices here
void EmbeddingOptimization::buildEquivalenceClasses() {
    // Given a triangle mesh, our goal is subdivide it into a quad mesh
    // by first splitting each triangle into 3 quads (by connecting the
    // barycenter of each triangle with the midpoints of each side), then
    // subdividing each quad into an n by n grid. To do so, we associate each
    // quad obtained in the first step with its associated triangle corner,
    // indexing within this quad as follows:
    //        (0,n-1) -> ... -> (n-1,n-1)
    //        /                    /
    //      .........................
    //     /                      /
    //    (0,1)               (n-1,1)
    //   /                      /
    // (0,0) -> (1,0) -> ... ->(n-1,0)
    //
    //  Flattening this coordinate system, vertex (x,y) on the cth corner is
    //  then assigned index c*n^2 + n*x + y.
    //
    //  Note that some vertices created in this process are redundant.
    //  In order to assemble our final quad mesh, we need to identify
    //  vertices that are double counted with this index system. Working it out,
    //  the correspondence is as given below, suppressing the corner-based
    //  offsets on each quad:
    //
    //        n*(n-1) -> ... -> n^2 - 1
    //        /                    /
    //      .........................
    //     /                      /
    //    n               2*n - 1 <-> n^2-n+1
    //   /                      /
    //   0 -> 1 -> ... -> n - 1  <-> n^2 - n
    //   ^    ^
    //   |    |
    //   v    v
    //   0    n
    for (Corner c : mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        // merge the "interior edge" shared with the next corner
        size_t cNext = c_[c.halfedge().next().corner()] * n * n;
        for (int i = 0; i < n; i++) {
            merge(cOffset + n * i + n - 1, cNext + n * n - n + i);
        }
        // merge with the corner on the same edge
        if (!c.halfedge().edge().isBoundary()) {
            size_t cOpp = c_[c.halfedge().twin().next().corner()] * n * n;
            for (int i = 0; i < n; i++) {
                merge(cOffset + i, cOpp + n * i);
            }
        }
    }
}
// Helper functions for buildEquivalenceClasses()
void EmbeddingOptimization::merge(int a, int b) {
    // Union find-like method that collapes the equivalence classes of a and b
    // Somewhat important detail: always chooses the lower index as the new
    // representative of the equivalence class
    if (top[a] > top[b])
        std::swap(a, b);
    if (top[a] == top[b])
        return;
    b = top[b];
    // Insert b and it's descendants in between a and aNext
    int aNext = next[a];
    next[a] = b;
    // find everything pointed to by b and make it point to a
    while (next[b] != -1) {
        top[b] = top[a];
        b = next[b];
    }
    // finalize changes
    top[b] = top[a];
    next[b] = aNext;
}

// returns canonical index of [a]
int EmbeddingOptimization::find(int a) { return top[a]; }

// returns the final number of equivalence classes.
int EmbeddingOptimization::buildFinalIndices() {
    // reindexing so all our final vertex indices are contiguous,
    // which is required by geometry-central.
    // finalIndices[i] gives the canonical compressed index of vertex i.
    std::map<int, int> reindex;
    size_t count = 0;
    for (int i = 0; i < nCorners * n * n; i++) {
        if (!reindex.count(find(i))) {
            reindex[find(i)] = count;
            count++;
        }
        finalIndices.push_back(reindex[find(i)]);
    }
    return count;
}

void EmbeddingOptimization::buildSubdivision() {
    // now all the indexing garbage is done...let's construct the mesh and setup
    // the energy
    vector<vector<size_t>> polygons;
    for (Corner c : mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                // each quad within a corner quad is given by
                // {(i,j), (i+1,j), (i+1,j+1), (i,j+1)} in local coordinates
                size_t i = finalIndices[cOffset + x + n * y];
                size_t j = finalIndices[cOffset + (x + 1) + n * y];
                size_t k = finalIndices[cOffset + (x + 1) + n * (y + 1)];
                size_t l = finalIndices[cOffset + x + n * (y + 1)];
                polygons.push_back({i, j, k, l});
            }
        }
    }
    submesh =
        std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(polygons));
    VertexData<Vector3> positions(*submesh);
    for (Corner c : mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                positions[finalIndices[cOffset + x + n * y]] = bary(c, x, y);
            }
        }
    }
    subgeometry = std::unique_ptr<VertexPositionGeometry>(
        new VertexPositionGeometry(*submesh, positions));
    writeSurfaceMesh(*submesh, *subgeometry, "my_mesh.obj");
    // TODO: change this call to solve()
    polyscope::registerSurfaceMesh("New mesh", positions,
                                   submesh->getFaceVertexList());
    return;
}

// writes the expected edge lengths/angles associated to each quad to c_iso_i
// for i = 0,1,2. Also initializes the vectors quads and fairVertices, which
// store the indices/permutations over which the energy is taken c_iso_0 =
// |ik|^2, c_iso_1 = |jl|^2, c_iso_2 = < ki, lj >
//
//    l --------k
//   /          /
//  /          /
// i -------- j

void EmbeddingOptimization::buildIntrinsicCheckerboard() {
    // ==============Isometry energy stuff======================
    // Allocate space for target edge lengths. Note that these are lengths per
    // quad, hence the count is different from the number of vertices number of
    // quads in the subdivided mesh
    size_t nQuads = nCorners * (n - 1) * (n - 1);
    c_iso_0 = vector<double>(nQuads);
    c_iso_1 = vector<double>(nQuads);
    c_iso_2 = vector<double>(nQuads);
    // storing the associated permutation of indices associated to each quad
    quads = vector<tuple<size_t, size_t, size_t, size_t>>(nQuads);

    for (Corner c : mesh->corners()) {
        // initialize the target lengths per quad as well as the appropriate
        // permutation
        size_t cQuadOffset =
            c_[c] * (n - 1) * (n - 1); // base index for coarse quad faces
        size_t cVertexOffset =
            c_[c] * n * n; // base index for coarse quad vertices

        // grab face coordinates
        Halfedge IJ = c.halfedge();
        Vector3 I = geometry->inputVertexPositions[IJ.vertex()];
        Halfedge JK = IJ.next();
        Vector3 J = geometry->inputVertexPositions[JK.vertex()];
        Halfedge KI = JK.next();
        Vector3 K = geometry->inputVertexPositions[KI.vertex()];

        // solve for coefficients associated to input data
        double bij = beta[IJ.corner()];
        double bjk = beta[JK.corner()];
        double bki = beta[KI.corner()];
        BezierTriangle T = Coefficients(I, J, K, beta[IJ.corner()],
                                        beta[JK.corner()], beta[KI.corner()]);
        // iterate over fine quads
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                // get face index for current fine quad
                size_t quadIndex = cQuadOffset + x + (n - 1) * y;

                // get vertex indices for current fine quad ijkl
                size_t iIndex = finalIndices[cVertexOffset + x + (n)*y];
                size_t jIndex = finalIndices[cVertexOffset + (x + 1) + (n)*y];
                size_t kIndex =
                    finalIndices[cVertexOffset + (x + 1) + n * (y + 1)];
                size_t lIndex = finalIndices[cVertexOffset + x + n * (y + 1)];

                // store the permutation
                quads[quadIndex] = {iIndex, jIndex, kIndex, lIndex};

                // read in target edge lengths
                // For each quad ijkl, the energy is
                //    l --------k
                //   /          /
                //  /          /
                // i -------- j
                //
                // c_iso_0 = |v_i - v_k|^2 - Lik^2,
                // c_iso_1 = |v_j - v_l|^2 - Ljl^2,
                // c_iso_2 = < v_k - v_i, v_l - v_j > - θijkl

                Vector2 i = RationalBezierTriangle(T, baryCoords(x, y));
                Vector2 j = RationalBezierTriangle(T, baryCoords(x + 1, y));
                Vector2 k = RationalBezierTriangle(T, baryCoords(x + 1, y + 1));
                Vector2 l = RationalBezierTriangle(T, baryCoords(x, y + 1));
                Vector2 e20 = i - k;
                Vector2 e31 = j - l;
                c_iso_0[quadIndex] = e20.norm2() * 1;
                c_iso_1[quadIndex] = e31.norm2() * 1;
                c_iso_2[quadIndex] = dot(e20, e31) * 1;
            }
        }
    }
    // ===============================Fairness energy
    // stuff================================== number of vertices for which we
    // calculate the fairness term allocate space
    fairVertices.clear();

    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    // "graph laplacian" type energy
    /*
    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Vertex i: submesh->vertices()) {
        vector<size_t> neighbors = {subVertexIndices[i]};
        for (Vertex j: i.adjacentVertices()) {
            neighbors.push_back(subVertexIndices[j]);
        }
        fairVertices.push_back(neighbors);
    }
    return;
    */
    for (Corner c : mesh->corners()) {
        // Recall that the subdivided grid looks like the following:
        //        (0,n-1) -> ... -> (n-1,n-1)
        //        /                    /
        //      .........................
        //     /                      /
        //    (0,1)               (n-1,1)
        //   /                      /
        // (0,0) -> (1,0) -> ... ->(n-1,0)
        //
        // For every vertex i with two neighbors j and k following it in the
        // associated row/column, we will add a fairness term ||i - 2j + k|| to
        // the energy To avoid double counting edges shared between adjacent
        // quads, for such shared edges, we only count them once when they
        // appear as a "vertical" edge, which corresponds to only having y range
        // over [1,n-1) in the first loop below. The total number of resulting
        // terms in the energy per corner is 3*((n-2)*(n-2) + n*(n-2))

        size_t cVertexOffset =
            c_[c] * n * n; // base index for coarse quad vertices
        // regularization term - horizontal
        for (int x = 0; x < n - 2; x++) {
            for (int y = 1; y < n - 1; y++) {
                size_t i = finalIndices[cVertexOffset + x + n * y];
                size_t j = finalIndices[cVertexOffset + (x + 1) + n * y];
                size_t k = finalIndices[cVertexOffset + (x + 2) + n * y];
                fairVertices.push_back({j, i, k});
            }
        }
        // regularization term - vertical
        // Note x is allowed to take the value n-1 here
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n - 2; y++) {
                size_t i = finalIndices[cVertexOffset + x + n * y];
                size_t j = finalIndices[cVertexOffset + x + n * (y + 1)];
                size_t k = finalIndices[cVertexOffset + x + n * (y + 2)];
                fairVertices.push_back({j, i, k});
            }
        }
    }
    // regularizer for the center vertex on each edge
    // I do some probably unnecessarily complicated stuff here to account for
    // boundary; I should check if it's necessary.
    for (Edge e : mesh->edges()) {
        // grab the two associated corners
        Corner c1 = e.halfedge().isInterior() ? e.halfedge().corner()
                                              : e.halfedge().twin().corner();
        Corner c2 = c1.halfedge().next().corner();
        // find the index for the midpoint vertex
        size_t centerIndex = finalIndices[c_[c1] * n * n + (n - 1)];
        size_t leftIndex = finalIndices[c_[c1] * n * n + (n - 2)];
        size_t rightIndex = finalIndices[c_[c2] * n * n + n * (n - 2)];
        fairVertices.push_back({centerIndex, leftIndex, rightIndex});
    }
    // Regularizer for edges like the following:
    //        _____________________
    //       /                    /
    //      .....................
    //     /                    /
    //    /         (n-2) <-> (n-1,y) <-> (y, n-2) (on the other corner)
    //   /                    /
    //   ---------------------
    //
    for (Corner c1 : mesh->corners()) {
        Corner c2 = c1.halfedge().next().corner();
        for (int y = 1; y < n - 1; y++) {
            size_t centerIndex = finalIndices[c_[c1] * n * n + (n - 1) + n * y];
            size_t leftIndex = finalIndices[c_[c1] * n * n + (n - 2) + n * y];
            size_t rightIndex = finalIndices[c_[c2] * n * n + y + n * (n - 2)];
            fairVertices.push_back({centerIndex, leftIndex, rightIndex});
        }
    }

    // "triple point" regularizer
    // for each face, we add a regularization term for the center vertex and
    // it's 3 neighbors, which are all in distinct "corner quads"
    for (Face f : mesh->faces()) {
        Corner c = f.halfedge().corner();
        size_t centerVertexIndex =
            finalIndices[c_[c] * n * n + (n - 1) + (n - 1) * n];
        Vertex i = submesh->vertex(centerVertexIndex);
        vector<size_t> neighbors = {centerVertexIndex};
        for (Vertex j : i.adjacentVertices()) {
            neighbors.push_back(subVertexIndices[j]);
        }
        fairVertices.push_back(neighbors);
    }

    // regularizer for vertices of original mesh
    /*
    for (Vertex v0: mesh->vertices()) {
        size_t subdividedIndex = finalIndices[c_[v0.corner()] * n * n];
        Vertex i = submesh->vertex(subdividedIndex);
        vector<size_t> neighbors = {subdividedIndex};
        for (Vertex j: i.adjacentVertices()) {
            neighbors.push_back(subVertexIndices[j]);
        }
        fairVertices.push_back(neighbors);
    }
    */

    // assert(fairVertices.size() == nCorners * ((n-2)*(n-2) + n*(n-2)));
}
//=========================Indexing convenience functions=======================
// ======================= Basis Function Stuff=========================
double EmbeddingOptimization::Angle(Vector2 u, Vector2 v) {
    return atan2(cross(u, v), dot(u, v));
}

// solves for angle offsets alpha given Euclidean angles ti and beta values
std::tuple<double, double, double>
EmbeddingOptimization::bendAngles(double t1, double t2, double t3, double b1,
                                  double b2, double b3) {
    double aij = (b1 + b2 - b3 - t1 - t2 + t3) / 2;
    double ajk = (-b1 + b2 + b3 + t1 - t2 - t3) / 2;
    double aki = (b1 - b2 + b3 - t1 + t2 - t3) / 2;
    return {aij, ajk, aki};
}

// isometrically projects triangle i j k to the plane.
// for convenience, j is assumed to be on the x-axis, and i is set to 0
std::tuple<Vector2, Vector2, Vector2>
EmbeddingOptimization::projectToPlane(Vector3 i, Vector3 j, Vector3 k) {
    j -= i;
    k -= i;
    i = Vector3{0., 0., 0.};
    /*
       Vector3 U = j.normalize();
       Vector3 V = cross(j,cross(j,k)).normalize();

       return {{0.,0.}, {j.norm(),0}, {dot(U, k), dot(V,k)}};
       */
    double lij = j.norm();
    double ljk = (k - j).norm();
    double lki = k.norm();
    double thetai = acos((lij * lij - ljk * ljk + lki * lki) / (2 * lij * lki));
    Vector2 I = {0., 0.};
    Vector2 J = {j.norm(), 0.};
    Vector2 K = {lki, 0.};
    K = K.rotate(-thetai);
    assert(abs((K - J).norm() - (k - j).norm()) < 1e-6);
    return {I, J, K};
}

BezierTriangle EmbeddingOptimization::Coefficients(Vector3 I, Vector3 J,
                                                   Vector3 K, double Bi,
                                                   double Bj, double Bk) {
    auto [i, j, k] = projectToPlane(I, J, K);
    Vector2 eij = i - j;
    Vector2 ejk = j - k;
    Vector2 eki = k - i;
    Vector2 nij = eij.rotate90();
    Vector2 njk = ejk.rotate90();
    Vector2 nki = eki.rotate90();

    double thetai = Angle(-nij, nki);
    double thetaj = Angle(-njk, nij);
    double thetak = Angle(-nki, njk);

    auto [aij, ajk, aki] = bendAngles(thetai, thetaj, thetak, Bi, Bj, Bk);

    Vector2 mij = (i + j) / 2;
    Vector2 mjk = (j + k) / 2;
    Vector2 mki = (k + i) / 2;

    Vector2 p200 = i;
    Vector2 p020 = j;
    Vector2 p002 = k;
    double w200 = 1, w020 = 1, w002 = 1;

    Vector2 p110 = mij + tan(aij) * nij / 2;
    double w110 = cos(aij);

    Vector2 p011 = mjk + tan(ajk) * njk / 2;
    double w011 = cos(ajk);

    Vector2 p101 = mki + tan(aki) * nki / 2;
    double w101 = cos(aki);

    return {p200, p020, p002, p110, p011, p101,
            w200, w020, w002, w110, w011, w101};
}

BezierTriangle EmbeddingOptimization::rotIndices(BezierTriangle T) {
    return {T.p002, T.p200, T.p020, T.p101, T.p110, T.p011,
            T.w020, T.w002, T.w200, T.w101, T.w110, T.w011};
}

Vector2 EmbeddingOptimization::RationalBezierTriangle(
    BezierTriangle T, std::tuple<double, double, double> coords) {
    auto [t1, t2, t3] = coords;

    double B200 = t1 * t1;
    double B020 = t2 * t2;
    double B002 = t3 * t3;
    double B110 = 2 * t1 * t2;
    double B011 = 2 * t2 * t3;
    double B101 = 2 * t1 * t3;
    Vector2 y = B200 * T.w200 * T.p200 + B020 * T.w020 * T.p020 +
                B002 * T.w002 * T.p002 + B110 * T.w110 * T.p110 +
                B011 * T.w011 * T.p011 + B101 * T.w101 * T.p101;
    double h = B200 * T.w200 + B020 * T.w020 + B002 * T.w002 + B110 * T.w110 +
               B011 * T.w011 + B101 * T.w101;
    return y / h;
}

// ======================= Optimization Stuff==========================
// Given a corner c and integer weights X,Y
// representing local coordinates on c's quad, returns the associated point in
// R^3. Here the integers X,Y are assumed to be in the range [0,N-1], with the
// mapping to barycentric coordinates being given by the transformation induced
// by (0,0)->(1,0,0), (N-1,N-1)->(1/3,1/3,1/3), (0,N-1)->(1/2,0,1/2) and linear
// interpolation in between
Vector3 EmbeddingOptimization::bary(Corner c, int X, int Y) {
    // Normalize x and y
    double x = (double)X / (n - 1);
    double y = (double)Y / (n - 1);

    // grab face coordinates
    Halfedge ij = c.halfedge();
    Vector3 i = geometry->inputVertexPositions[ij.vertex()];
    Halfedge jk = ij.next();
    Vector3 j = geometry->inputVertexPositions[jk.vertex()];
    Halfedge ki = jk.next();
    Vector3 k = geometry->inputVertexPositions[ki.vertex()];

    double jWeight = -(y - 3) * x / 6;
    double kWeight = -(x - 3) * y / 6;
    return (1 - jWeight - kWeight) * i + jWeight * j + kWeight * k;
}

// a convenience function analagous to above that returns the barycentric
// coordinates associated to a point on the grid
std::tuple<double, double, double> EmbeddingOptimization::baryCoords(int X,
                                                                     int Y) {
    // Normalize x and y
    double x = (double)X / (n - 1);
    double y = (double)Y / (n - 1);

    double jWeight = -(y - 3) * x / 6;
    double kWeight = -(x - 3) * y / 6;
    return {(1 - jWeight - kWeight), jWeight, kWeight};
}

void EmbeddingOptimization::basisFunctionDebugging() {
    double N = 2 * n * (n - 1);
    // ==================================TEST============================================
    std::ofstream ss("test.svg", std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" " << endl
       << "<svg width=\"1000\" height=\"1000\" "
          "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >"
       //<< "<rect width=\"100%\" height =\"100%\" fill=\"#ffffff\"/>"
       << endl;
    BezierTriangle T0 = Coefficients({0, 0, 0}, {400, 0, 0}, {0, 400, 0},
                                     2 * PI / 3, PI / 4, PI / 4);
    for (int x = 0; x <= N; x++) {
        for (int y = 0; x + y <= N; y++) {
            Vector2 i =
                RationalBezierTriangle(T0, {x / N, y / N, (1 - x / N - y / N)});
            cout << i << endl;
            ss << "<circle cx=\"" << 500 + (i.x) << "\" cy=\"" << 500 + (i.y)
               << "\" r=\"3\" fill=\"red\"/>" << endl;
        }
    }

    vector<BezierTriangle> temp = {T0, rotIndices(T0),
                                   rotIndices(rotIndices(T0))};
    for (auto T : temp) {
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                Vector2 i = RationalBezierTriangle(T, baryCoords(x, y));
                // cout << i << endl;
                // Vector2 j = RationalBezierTriangle(T,baryCoords(x+1, y));
                // Vector2 k = RationalBezierTriangle(T,baryCoords(x+1, y+1));
                // Vector2 l = RationalBezierTriangle(T,baryCoords(x, y+1));
                ss << "<circle cx=\"" << 500 + (i.x) << "\" cy=\""
                   << 500 + (i.y) << "\" r=\"3\" fill=\"green\"/>" << endl;
            }
        }
    }

    ss << "</svg>";
}
