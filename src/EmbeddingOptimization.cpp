#include "EmbeddingOptimization.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


// ======================= Basis Function Stuff ==================================
double EmbeddingOptimization::Angle (Vector2 u, Vector2 v) {
    return atan2(cross(u,v), dot(u,v));
}


// solves for angle offsets alpha given Euclidean angles ti and beta values
std::tuple<double, double, double> EmbeddingOptimization::bendAngles(double t1, double t2, double t3, double b1, double b2, double b3) {
    double aij = (b1+b2-b3-t1-t2+t3)/2;
    double ajk = (-b1+b2+b3+t1-t2-t3)/2;
    double aki = (b1-b2+b3-t1+t2-t3)/2;
    return {aij, ajk, aki};
}

// isometrically projects triangle i j k to the plane.
// for convenience, j is assumed to be on the x-axis, and i is set to 0
std::tuple<Vector2,Vector2,Vector2> EmbeddingOptimization::projectToPlane(Vector3 i, Vector3 j, Vector3 k) {
    j -= i; k -= i;
    i = Vector3{0.,0.,0.};
    /*
       Vector3 U = j.normalize();
       Vector3 V = cross(j,cross(j,k)).normalize();

       return {{0.,0.}, {j.norm(),0}, {dot(U, k), dot(V,k)}};
       */
    double lij = j.norm();
    double ljk = (k-j).norm();
    double lki = k.norm();
    double thetai = acos((lij*lij - ljk*ljk + lki*lki)/(2*lij*lki));
    Vector2 I = {0.,0.};
    Vector2 J = {j.norm(),0.};
    Vector2 K = {lki,0.};
    K = K.rotate(-thetai);
    assert(abs((K-J).norm() - (k-j).norm()) < 1e-6);
    return {I,J,K};
}

BezierTriangle EmbeddingOptimization::Coefficients (Vector3 I, Vector3 J, Vector3 K, double Bi, double Bj, double Bk) {
    auto [i,j,k] = projectToPlane(I,J,K);
    Vector2 eij = i-j;
    Vector2 ejk = j-k;
    Vector2 eki = k-i;
    Vector2 nij = eij.rotate90();
    Vector2 njk = ejk.rotate90();
    Vector2 nki = eki.rotate90();

    double thetai = Angle(-nij,nki);
    double thetaj = Angle(-njk,nij);
    double thetak = Angle(-nki,njk);

    auto [aij,ajk,aki] = bendAngles(thetai,thetaj,thetak,Bi,Bj,Bk);

    Vector2 mij = (i+j)/2;
    Vector2 mjk = (j+k)/2;
    Vector2 mki = (k+i)/2;

    Vector2 p200 = i;
    Vector2 p020 = j;
    Vector2 p002 = k;
    double w200 = 1, w020 = 1, w002 = 1;

    Vector2 p110 = mij + tan(aij)*nij/2;
    double w110 = cos(aij);

    Vector2 p011 = mjk + tan(ajk)*njk/2;
    double w011 = cos(ajk);

    Vector2 p101 = mki + tan(aki)*nki/2;
    double w101 = cos(aki);

    return { p200, p020, p002, p110, p011, p101, w200, w020, w002, w110, w011, w101 };
}

BezierTriangle EmbeddingOptimization::rotIndices(BezierTriangle T) {
    return {T.p002, T.p200, T.p020, T.p101,T.p110,T.p011,T.w020, T.w002, T.w200, T.w101,T.w110,T.w011};
}

Vector2 EmbeddingOptimization::RationalBezierTriangle(BezierTriangle T, std::tuple<double,double,double> coords) {
    auto [t1,t2,t3] = coords;

    double B200 = t1 * t1;
    double B020 = t2 * t2;
    double B002 = t3 * t3;
    double B110 = 2 * t1 * t2;
    double B011 = 2 * t2 * t3;
    double B101 = 2 * t1 * t3;
    Vector2 y = 
        B200 * T.w200 * T.p200 +
        B020 * T.w020 * T.p020 +
        B002 * T.w002 * T.p002 +
        B110 * T.w110 * T.p110 +
        B011 * T.w011 * T.p011 +
        B101 * T.w101 * T.p101;
    double h = 
        B200 * T.w200 +
        B020 * T.w020 +
        B002 * T.w002 +
        B110 * T.w110 +
        B011 * T.w011 +
        B101 * T.w101;
    return y/h;
}


// ======================= Optimization Stuff ====================================
// Given a corner c and integer weights X,Y representing local coordinates 
// on c's quad, returns the associated point in R^3.
// Here the integers X,Y are assumed to be in the range
// [0,N-1], with the mapping to barycentric
// coordinates being given by the transformation induced by
// (0,0)->(1,0,0), (N-1,N-1)->(1/3,1/3,1/3), (0,N-1)->(1/2,0,1/2)
// and linear interpolation in between
Vector3 EmbeddingOptimization::bary(Corner c, int X, int Y) {
    // Normalize x and y
    double x = (double)X/(n-1);
    double y = (double)Y/(n-1);

    // grab face coordinates
    Halfedge ij = c.halfedge();
    Vector3 i = geometry->inputVertexPositions[ij.vertex()];
    Halfedge jk = ij.next();
    Vector3 j = geometry->inputVertexPositions[jk.vertex()];
    Halfedge ki = jk.next();
    Vector3 k = geometry->inputVertexPositions[ki.vertex()];


    double jWeight = -(y-3) * x/6;
    double kWeight = -(x-3) * y/6;
    return (1 - jWeight - kWeight) * i + jWeight * j + kWeight * k;
}

// a convenience function analagous to above that returns the barycentric coordinates associated to
// a point on the grid
std::tuple<double, double, double> EmbeddingOptimization::baryCoords(int X, int Y) {
    // Normalize x and y
    double x = (double)X/(n-1);
    double y = (double)Y/(n-1);

    double jWeight = -(y-3) * x/6;
    double kWeight = -(x-3) * y/6;
    return {(1 - jWeight - kWeight), jWeight, kWeight};
}

// constructor
EmbeddingOptimization::EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, 
        shared_ptr<VertexPositionGeometry> geometry, CornerData<double> beta) 
    : mesh(mesh), geometry(geometry), beta(beta) {
        // Initialize quantities
        intrinsicQuantitiesInitialized = false;
        fairnessNormalization = 0.5;
        LMInitialized = false;
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
           auto[i,j,k] = projectToPlane({0,0,0}, {3,0,0}, {0,4.,0});
           cout << i << j << k;
           */
    }


// ======================= Subdivision code ====================================
// Union find-like method that collapes the equivalence classes of a and b
// Somewhat important detail: always chooses the lower index as the new representative of the equivalence class
void EmbeddingOptimization::merge(int a, int b) {
    if (top[a] > top[b]) std::swap(a,b);
    if (top[a] == top[b]) return;
    b = top[b];
    int aNext = next[a];
    next[a] = b;
    while (next[b] != -1) {
        top[b] = top[a];
        b = next[b];
    }
    top[b] = top[a];
    next[b] = aNext;
}

// returns canonical index of [a]
int EmbeddingOptimization::find(int a) {
    return top[a];
}

void EmbeddingOptimization::buildEquivalenceClasses() {
    // Given a triangle mesh, our goal is subdivide it into a quad mesh
    // by first splitting each triangle into 3 quads (by connecting the barycenter
    // of each triangle with the midpoints of each side), then subdividing each quad 
    // into an n by n grid.
    // To do so, we associate each quad obtained in the first step with its associated
    // triangle corner, indexing within this quad as follows:
    //        (0,n-1) -> ... -> (n-1,n-1)
    //        /                    /
    //      .........................
    //     /                      /
    //    (0,1)               (n-1,1)
    //   /                      /
    // (0,0) -> (1,0) -> ... ->(n-1,0)
    //
    //  Flattening this coordinate system, vertex (x,y) on the cth corner is then 
    //  assigned index c*n^2 + n*x + y.
    //
    //  Note that some vertices created in this process are redundant.
    //  In order to assemble our final quad mesh, we hence need to identify vertices
    //  that are double counted with this index system.
    //  Working it out, the correspondence is as given below, suppressing 
    //  the corner-based offsets on each quad:
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
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        // merge the "interior edge" shared with the next corner
        size_t cNext = c_[c.halfedge().next().corner()] * n * n;
        for (int i = 0; i < n; i++) {
            merge(cOffset + n*i + n - 1, cNext + n*n -n + i);
        }
        // merge with the corner on the same edge
        if(!c.halfedge().edge().isBoundary()) {
            size_t cOpp = c_[c.halfedge().twin().next().corner()] * n * n;
            for (int i = 0; i < n; i++) {
                merge(cOffset + i, cOpp + n*i);
            }
        }
    }
}

// returns the final number of equivalence classes.
int EmbeddingOptimization::buildFinalIndices() {
    // reindexing so all our final vertex indices are contiguous,
    // which is required by geometry-central.
    // finalIndices[i] gives the canonical compressed index of vertex i.
    std::map<int,int> reindex;
    size_t temp = 0;
    for (int i = 0; i < nCorners * n * n; i++) {
        if (!reindex.count(find(i))){
            reindex[find(i)] = temp;
            temp++;
        }
        finalIndices.push_back(reindex[find(i)]);
    }
    return temp;
}


void EmbeddingOptimization::buildSubdivision(){
    // now all the indexing garbage is done...let's construct the mesh and setup the energy
    vector<vector<size_t>> polygons;
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                // each quad within a corner quad is given by
                // {(i,j), (i+1,j), (i+1,j+1), (i,j+1)} in local coordinates
                size_t i = finalIndices[cOffset + x + n * y];
                size_t j = finalIndices[cOffset + (x+1) + n * y];
                size_t k = finalIndices[cOffset + (x+1) + n * (y+1)];
                size_t l = finalIndices[cOffset + x + n * (y+1)];
                polygons.push_back({i,j,k,l});
            }
        }
    }
    submesh = std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(polygons));
    VertexData<Vector3> positions(*submesh);
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                positions[finalIndices[cOffset + x + n * y]] = bary(c, x, y);
            }
        }
    }
    subgeometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*submesh,positions));
    writeSurfaceMesh(*submesh, *subgeometry, "my_mesh.obj");
    // TODO: change this call to solve()
    polyscope::registerSurfaceMesh("New mesh", positions, submesh->getFaceVertexList());
    return;
}


// writes the expected edge lengths/angles associated to each quad to c_iso_i for
// i = 0,1,2. 
// Also initializes the vectors quads and fairVertices, which store the indices/permutations over which the energy
// is taken
// c_iso_0 = |ik|^2, 
// c_iso_1 = |jl|^2, 
// c_iso_2 = < ki, lj >
//
//    l --------k
//   /          /
//  /          /
// i -------- j

void EmbeddingOptimization::buildIntrinsicCheckerboard(){
    // ===============================Isometry energy stuff==================================
    // Allocate space for target edge lengths
    // Note that these are lengths per quad, hence the count is different from the number of vertices
    // number of quads in the subdivided mesh
    size_t nQuads = nCorners * (n-1) * (n-1);
    c_iso_0 = vector<double>(nQuads);
    c_iso_1 = vector<double>(nQuads);
    c_iso_2 = vector<double>(nQuads);
    // storing the associated permutation of indices associated to each quad
    quads = vector<tuple<size_t, size_t, size_t, size_t>>(nQuads);

    for (Corner c: mesh->corners()) {
        // initialize the target lengths per quad as well as the appropriate permutation
        size_t cQuadOffset = c_[c] * (n-1) * (n-1); // base index for coarse quad faces
        size_t cVertexOffset = c_[c] * n * n; // base index for coarse quad vertices

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
        BezierTriangle T = Coefficients(I,J,K,beta[IJ.corner()],beta[JK.corner()],beta[KI.corner()]);
        // iterate over fine quads
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                // get face index for current fine quad
                size_t quadIndex = cQuadOffset + x + (n-1)*y;

                // get vertex indices for current fine quad ijkl
                size_t iIndex = finalIndices[cVertexOffset + x + (n)*y];
                size_t jIndex = finalIndices[cVertexOffset + (x+1) + (n)*y];
                size_t kIndex = finalIndices[cVertexOffset + (x+1) + n*(y+1)];
                size_t lIndex = finalIndices[cVertexOffset + x + n*(y+1)];

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


                Vector2 i = RationalBezierTriangle(T,baryCoords(x, y));
                Vector2 j = RationalBezierTriangle(T,baryCoords(x+1, y));
                Vector2 k = RationalBezierTriangle(T,baryCoords(x+1, y+1));
                Vector2 l = RationalBezierTriangle(T,baryCoords(x, y+1));
                Vector2 e20 = i - k;
                Vector2 e31 = j - l;
                c_iso_0[quadIndex] = e20.norm2()/2;
                c_iso_1[quadIndex] = e31.norm2()/2;
                c_iso_2[quadIndex] = dot(e20,e31)/2;
            }
        }
    }
    // ===============================Fairness energy stuff==================================
    // number of vertices for which we calculate the fairness term
    // size_t nFairVertices = nCorners * ((n-2)*(n-2) + n*(n-2));
    // allocate space
    fairVertices.clear();

    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Vertex i: submesh->vertices()) {
        vector<size_t> neighbors = {subVertexIndices[i]};
        for (Vertex j: i.adjacentVertices()) {
            neighbors.push_back(subVertexIndices[j]);
        }
        fairVertices.push_back(neighbors);
    }
    /*
     for (Corner c: mesh->corners()) {
    // Recall that the subdivided grid looks like the following:
    //        (0,n-1) -> ... -> (n-1,n-1)
    //        /                    /
        //      .........................
        //     /                      /
        //    (0,1)               (n-1,1)
        //   /                      /
        // (0,0) -> (1,0) -> ... ->(n-1,0)
        //
        // For every vertex i with two neighbors j and k following it in the associated row/column, 
        // we will add a fairness term ||i - 2j + k|| to the energy
        // To avoid double counting edges shared between adjacent quads, for such shared edges, we only count 
        // them once when they appear as a "vertical" edge, which corresponds to only having y range over [1,n-1)
        // in the first loop below.
        // The total number of resulting terms in the energy per corner is 3*((n-2)*(n-2) + n*(n-2))

        size_t cVertexOffset = c_[c] * n * n; // base index for coarse quad vertices
        size_t cOffset = 3 * c_[c] * ((n-2)*(n-2) + n*(n-2)); // local offset in this energy term
                                                              // regularization term - horizontal
        for (int x = 0; x < n - 2; x++) {
            for (int y = 1; y < n - 1; y++) {
                //size_t energyIndex = cOffset + 3*(x + (n-2) * (y-1));
                size_t i = finalIndices[cVertexOffset + x + n*y];
                size_t j = finalIndices[cVertexOffset + (x+1) + n*y];
                size_t k = finalIndices[cVertexOffset + (x+2) + n*y];
                fairVertices.push_back({j, i, k}); 
            }
        }
        // regularization term - vertical
        // Note x is allowed to take the value n-1 here
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n - 2; y++) {
                //size_t energyIndex = cOffset + 3 * (n-2) * (n-2) + 3*(x + n * y);
                size_t i = finalIndices[cVertexOffset + x + n*y];
                size_t j = finalIndices[cVertexOffset + x + n*(y+1)];
                size_t k = finalIndices[cVertexOffset + x + n*(y+2)];
                fairVertices.push_back({j, i, k}); 
            }
        }
    }
    //assert(fairVertices.size() == nCorners * ((n-2)*(n-2) + n*(n-2)));
    */
}

// ======================= Energy evaluation code ====================================
inline double EmbeddingOptimization::sqr(double x) { return x * x; }
// Given indices for vertices i, j
// adds the term (|i-j|^2 - target) and it's corresponding gradient to the energy
inline void EmbeddingOptimization::addLengthTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, size_t energyIndex, size_t iIndex, size_t jIndex, double target) {
    Vector3 vi = {v[3*iIndex], v[3*iIndex+1],v[3*iIndex+2]};
    Vector3 vj = {v[3*jIndex], v[3*jIndex+1],v[3*jIndex+2]};

    const double normalization = 1000;
    energy[energyIndex] = ((vi-vj).norm2() - target)*normalization;
}

inline void EmbeddingOptimization::addLengthGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, size_t energyIndex, size_t iIndex, size_t jIndex, double target) {
    typedef Eigen::Triplet<double> T;
    Vector3 i = {v[3*iIndex], v[3*iIndex+1],v[3*iIndex+2]};
    Vector3 j = {v[3*jIndex], v[3*jIndex+1],v[3*jIndex+2]};

    const double normalization = 1000;
    // i gradients
    tripletList.push_back(T(energyIndex,3*iIndex,   2 *(i.x - j.x) * normalization));
    tripletList.push_back(T(energyIndex,3*iIndex+1, 2 *(i.y - j.y) * normalization));
    tripletList.push_back(T(energyIndex,3*iIndex+2, 2 *(i.z - j.z) * normalization));

    // j gradients
    tripletList.push_back(T(energyIndex,3*jIndex,   2 *(j.x - i.x) * normalization));
    tripletList.push_back(T(energyIndex,3*jIndex+1, 2 *(j.y - i.y) * normalization));
    tripletList.push_back(T(energyIndex,3*jIndex+2, 2 *(j.z - i.z) * normalization));
}

inline void EmbeddingOptimization::addAngleTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, size_t energyIndex, size_t iIndex, size_t jIndex, size_t kIndex, size_t lIndex, double target) {
    Vector3 vi = {v[3*iIndex], v[3*iIndex+1], v[3*iIndex+2]};
    Vector3 vj = {v[3*jIndex], v[3*jIndex+1], v[3*jIndex+2]};
    Vector3 vk = {v[3*kIndex], v[3*kIndex+1], v[3*kIndex+2]};
    Vector3 vl = {v[3*lIndex], v[3*lIndex+1], v[3*lIndex+2]};
    const double normalization = 1000;
    energy[energyIndex] = (dot(vi-vk,vj-vl) - target) * normalization;
}

inline void EmbeddingOptimization::addAngleGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, size_t energyIndex, size_t iIndex, size_t jIndex, size_t kIndex, size_t lIndex, double target) {
    typedef Eigen::Triplet<double> T;
    Vector3 i = {v[3*iIndex], v[3*iIndex+1], v[3*iIndex+2]};
    Vector3 j = {v[3*jIndex], v[3*jIndex+1], v[3*jIndex+2]};
    Vector3 k = {v[3*kIndex], v[3*kIndex+1], v[3*kIndex+2]};
    Vector3 l = {v[3*lIndex], v[3*lIndex+1], v[3*lIndex+2]};

    const double normalization = 1000;
    // i partial: j - l
    tripletList.push_back(T(energyIndex,3*iIndex,   j.x - l.x * normalization));
    tripletList.push_back(T(energyIndex,3*iIndex+1, j.y - l.y * normalization));
    tripletList.push_back(T(energyIndex,3*iIndex+2, j.z - l.z * normalization));

    // k partial: l - j
    tripletList.push_back(T(energyIndex,3*kIndex,   l.x - j.x * normalization));
    tripletList.push_back(T(energyIndex,3*kIndex+1, l.y - j.y * normalization));
    tripletList.push_back(T(energyIndex,3*kIndex+2, l.z - j.z * normalization));

    // j partial: i - k
    tripletList.push_back(T(energyIndex,3*jIndex,   i.x - k.x * normalization));
    tripletList.push_back(T(energyIndex,3*jIndex+1, i.y - k.y * normalization));
    tripletList.push_back(T(energyIndex,3*jIndex+2, i.z - k.z * normalization));

    // l partial: k - i
    tripletList.push_back(T(energyIndex,3*lIndex,   k.x - i.x * normalization));
    tripletList.push_back(T(energyIndex,3*lIndex+1, k.y - i.y * normalization));
    tripletList.push_back(T(energyIndex,3*lIndex+2, k.z - i.z * normalization));
}

// indices is a vector consisting of [center term, neighbors].
// The corresponding contribution to the energy is || # neighbors * center position - Σ neighbors ||^2 
inline void EmbeddingOptimization::addCenterTerm(Eigen::VectorXd& energy, const Eigen::VectorXd& v, 
        size_t energyIndex, vector<size_t>& indices) {
    size_t nNeighbors = indices.size() - 1;
    size_t centerIndex = indices[0];
    for (int i = 0; i < 3; i++) {
        energy[energyIndex + i] -= nNeighbors * v[3*centerIndex + i] * fairnessNormalization;
    }

    for (int n = 1; n < indices.size(); n++) {
        for (int i = 0; i < 3; i++) {
            energy[energyIndex + i] += v[3*indices[n] + i] * fairnessNormalization;
        }
    }
}

inline void EmbeddingOptimization::addCenterGradient(vector<Eigen::Triplet<double>>& tripletList, const Eigen::VectorXd& v, 
        size_t energyIndex, vector<size_t>& indices) {
    typedef Eigen::Triplet<double> T;
    size_t nNeighbors = indices.size() - 1;
    size_t centerIndex = indices[0];
    for (int i = 0; i < 3; i++) {
        tripletList.push_back(T(energyIndex + i, 3 * centerIndex + i, -nNeighbors * fairnessNormalization));
    }

    for (int n = 1; n < indices.size(); n++) {
        for (int i = 0; i < 3; i++) {
            tripletList.push_back(T(energyIndex + i, 3 * centerIndex + i, fairnessNormalization));
        }
    }
}


inline double cube(double x) {
    return x * x * x;
}
void EmbeddingOptimization::LMOneStep(int MAX_ITERS) {
    if (!LMInitialized) {
        throw std::runtime_error("Error: LM stuff not initialized yet");
        return;
    }

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorType;
    // initialize parameters
    double mu = 1e-6;
    double nu = 2.;
    const double tau = 1e-3;
    const double eps1 = 1e-8;
    const double eps2 = 1e-8;
    size_t k = 0;

    // initialize the Jacobian
    Eigen::SparseMatrix<double> J(LMValues, LMInputs);
    J.setZero();
    evaluateJacobian(currentSolution, J);
    // compute the relevant things
    Eigen::SparseMatrix<double> A = J.transpose() * J;
    Eigen::VectorXd f = Eigen::VectorXd::Zero(LMValues);
    Eigen::VectorXd fNew = Eigen::VectorXd::Zero(LMValues);
    evaluateEnergy(currentSolution, f);
    cout << "INITIAL ENERGY" << f.norm() << endl;
    Eigen::VectorXd g = J.transpose() * f;
    bool found = (g.norm() < eps1);
    // initialize our solution currentSolution to the current vertex positions x
    VectorType xStar = currentSolution;

    while (!found && k < MAX_ITERS) {
        if (!(k % 50)) {
            cout << "ITERATION" << k << endl;
            cout << "New ENERGY" << f.norm() << endl;
        }
        k += 1;
        Eigen::SparseMatrix<double> DescentMatrix(LMInputs, LMInputs);
        DescentMatrix.setIdentity();
        DescentMatrix *= mu;
        DescentMatrix += A;
        //cout << "DescentMatrix norm: " << DescentMatrix.norm() << endl;
        Eigen::VectorXd hLM = -solvePositiveDefinite(DescentMatrix,g);
        // debug
        if (hLM.norm() <= eps2 * (xStar.norm() + eps2)) {
            cout << "TERMINATION CONDITION MET ON ITERATION " << k << endl;
            cout << "Final ENERGY" << f.norm() << endl;
            cout << hLM.norm() << endl;
            cout << xStar.norm() << endl;
            return;
        } else {
            xStar = currentSolution + hLM;
            // evaluate the energy after taking a descent step
            fNew.setZero();
            evaluateEnergy(xStar, fNew);
            double numerator = f.squaredNorm() - fNew.squaredNorm();
            double denominator = (hLM.dot(mu * hLM - g))/2.;
            double rho = numerator/denominator;
            if (rho > 0) {
                currentSolution = xStar;
                J.setZero();
                evaluateJacobian(currentSolution, J);
                // compute the relevant things
                A = J.transpose() * J;
                f = fNew;
                g = J.transpose() * f;
                found = (g.norm() < eps1);
                mu = mu * std::max(1./3., 1 - cube(2 * rho - 1));
                nu = 2;
            } else {
                mu = mu * nu;
                nu = 2 * nu;
            }
        }
    }
    // Solve the optimization problem
    cout << "Finished on iteration " << k << endl;

}

void EmbeddingOptimization::basisFunctionDebugging() {
    double N = 2*n*(n-1);
    // ==================================TEST============================================
    std::ofstream ss("test.svg", std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        << endl
        << "<svg width=\"1000\" height=\"1000\" "
        "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >"
        //<< "<rect width=\"100%\" height =\"100%\" fill=\"#ffffff\"/>"
        << endl;
    BezierTriangle T0 = Coefficients({0,0,0},{400,0,0},{0,400,0},2*PI/3,PI/4,PI/4);
    for (int x = 0; x <= N; x++) {
        for (int y = 0; x+y <= N; y++) {
            Vector2 i = RationalBezierTriangle(T0, {x/N, y/N,(1-x/N-y/N)});
            cout << i << endl;
            ss << "<circle cx=\"" << 500+(i.x) << "\" cy=\""
                << 500+(i.y) << "\" r=\"3\" fill=\"red\"/>" << endl;
        }
    }

    vector<BezierTriangle> temp = {T0, rotIndices(T0), rotIndices(rotIndices(T0))};
    for (auto T: temp) {
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                Vector2 i = RationalBezierTriangle(T,baryCoords(x, y));
                //cout << i << endl;
                //Vector2 j = RationalBezierTriangle(T,baryCoords(x+1, y));
                //Vector2 k = RationalBezierTriangle(T,baryCoords(x+1, y+1));
                //Vector2 l = RationalBezierTriangle(T,baryCoords(x, y+1));
                ss << "<circle cx=\"" << 500+(i.x) << "\" cy=\""
                    << 500+(i.y) << "\" r=\"3\" fill=\"green\"/>" << endl;
            }
        }
    }

    ss << "</svg>";

}

void EmbeddingOptimization::initializeLM() {
    if (!intrinsicQuantitiesInitialized) {
        throw std::runtime_error("Error: intrinsic checkboard/mesh not initialized yet");
        return;
    }
    currentSolution = initialSolution;

    int nCoordinates = 3*nSubdividedVertices;

    // number of quads in the fine mesh
    size_t nQuads = nCorners * (n-1) * (n-1);
    size_t nFairVertices = fairVertices.size();

    LMValues = 3*nQuads + 3 * nFairVertices;
    /*
    // the energy functor calls this->evaluateEnergy to evaluate the objective
    // and this->evaluateJacobian for the Jacobian
    EnergyFunctor EF(nCoordinates, values, this);
    */


    // make sure that the size of the vertex position vector x agrees with the
    // number of fine vertices in the mesh times three
    assert(currentSolution.size() == nCoordinates);
    LMInputs = nCoordinates;
    LMInitialized = true;
    return;
}
// Call the optimization procedure for one step, then update the mesh in
// polyscope with the new vertex positions
void EmbeddingOptimization::optimizeOneStep(int MAX_ITERS) {
    LMOneStep(MAX_ITERS);
    VertexData<Vector3> positions(*submesh);
    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Vertex v: submesh->vertices()) {
        size_t i = subVertexIndices[v];
        positions[v] = {currentSolution[3*i],currentSolution[3*i+1],currentSolution[3*i+2]};
        //cout << positions[v] << endl;
    }

    polyscope::registerSurfaceMesh("Optimized Mesh", positions, submesh->getFaceVertexList());
}

// ======================================= LM Stuff ==========================================
void EmbeddingOptimization::evaluateEnergy(
        const Eigen::VectorXd& v, // vertex positions
        Eigen::VectorXd& energy) // each entry of this vector is the (square root of) a term in the energy summand
{

    size_t nQuads = nCorners * (n-1) * (n-1);
    // make sure we were given the right number of vertex coordinates
    assert(v.size() == 3 * nSubdividedVertices);

    // make sure the given energy summand vector is the right size; there are 3 * nQuads
    assert(energy.size() == 3 * nQuads + 3 * fairVertices.size());

    // fill the energy summand vector with the individual terms
    for (size_t quadIndex = 0; quadIndex < nQuads; quadIndex++) {
        // vertex indices
        auto [i, j, k, l] = quads[quadIndex];
        addLengthTerm(energy, v, quadIndex, i, k, c_iso_0[quadIndex]);
        addLengthTerm(energy, v, quadIndex + nQuads, j, l, c_iso_1[quadIndex]);
        addAngleTerm(energy, v, quadIndex + 2*nQuads, i, j, k, l, c_iso_2[quadIndex]);
    }

    for (int i = 0; i < fairVertices.size(); i++) {
        addCenterTerm(energy, v, 3 * nQuads + 3 * i, fairVertices[i]);
    }
    

}
void EmbeddingOptimization::evaluateJacobian(const Eigen::VectorXd& v, Eigen::SparseMatrix<double>& J) {
    // number of fine quads
    //cout << "evaluating jacobian" << endl;
    size_t nQuads = nCorners * (n-1) * (n-1);
    vector<Eigen::Triplet<double>> tripletList;

    // isometry term entries
    for (size_t quadIndex = 0; quadIndex < nQuads; quadIndex++) {
        // vertex indices
        auto [i, j, k, l] = quads[quadIndex];
        addLengthGradient(tripletList, v, quadIndex, i, k, c_iso_0[quadIndex]);
        addLengthGradient(tripletList, v, quadIndex + nQuads, j, l, c_iso_1[quadIndex]);
        addAngleGradient(tripletList, v, quadIndex + 2*nQuads, i, j, k, l, c_iso_2[quadIndex]);
    }
    for (int i = 0; i < fairVertices.size(); i++) {
        addCenterGradient(tripletList, v, 3 * nQuads + 3 * i, fairVertices[i]);
    }
    // build Jacobian matrix from triplets
    J.setFromTriplets(tripletList.begin(), tripletList.end());
}

std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry> >EmbeddingOptimization::initializeSubdivision(int N) {
    n = N;

    // Initialize union find data structure
    top = vector<int> (n*n*nCorners);
    for (int i = 0; i < n * n * nCorners; i++) top[i] = i;
    next = vector<int> (n*n*nCorners, -1);

    buildEquivalenceClasses();
    cout << "built equivalence classes" << endl;
    nSubdividedVertices = buildFinalIndices();
    cout << "built final indices" << endl;

    // Initialize submesh and subgeometry
    buildSubdivision();
    cout << "built subdivision" << endl;
    buildIntrinsicCheckerboard();
    cout << "built intrinsic lengths" << endl;

    //===============DEBUG===================
    //basisFunctionDebugging();
    //return {submesh, subgeometry};
    // initialize inital guess for x
    initialSolution = Eigen::VectorXd::Zero(3*nSubdividedVertices);
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int X = 0; X < n; X++) {
            for (int Y = 0; Y < n; Y++) {
                Vector3 pos = bary(c, X, Y);
                size_t startIndex = 3*finalIndices[cOffset + X + n * Y];
                initialSolution[startIndex] = pos.x;
                initialSolution[startIndex + 1] = pos.y;
                initialSolution[startIndex + 2] = pos.z;
            }
        }
    }
    intrinsicQuantitiesInitialized = true;

    //polyscope::show();

    return {submesh, subgeometry};
}

