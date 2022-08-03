#include "fusion.h"
#include "EmbeddingOptimization.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace mosek::fusion;
using namespace monty;


// ======================= Basis Function Stuff ==================================
double Angle (Vector2 u, Vector2 v) {
    return atan2(cross(u,v), dot(u,v));
}
// solves for angle offsets alpha given Euclidean angles ti and beta values
std::tuple<double, double, double> bendAngles(double t1, double t2, double t3, double b1, double b2, double b3) {
    double aij = (b1+b2-b3-t1-t2+t3)/2;
    double ajk = (-b1+b2+b3+t1-t2-t3)/2;
    double aki = (b1-b2+b3-t1+t2-t3)/2;
    return {aij, ajk, aki};
}
// isometrically projects triangle i j k to the plane.
// for convenience, j is assumed to be on the x-axis, and i is set to 0
std::tuple<Vector2,Vector2,Vector2> projectToPlane(Vector3 i, Vector3 j, Vector3 k) {
    j -= i; k -= i;
    i -= i;
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

BezierTriangle Coefficients (Vector3 I, Vector3 J, Vector3 K, double Bi, double Bj, double Bk) {
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
    //cout << "midpoints" << mij << " " << mjk << " " << mki << endl;

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
BezierTriangle rotIndices(BezierTriangle T) {
    return {T.p002, T.p200, T.p020, T.p101,T.p110,T.p011,T.w020, T.w002, T.w200, T.w101,T.w110,T.w011};
}
Vector2 RationalBezierTriangle(BezierTriangle T, std::tuple<double,double,double> coords) {
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
EmbeddingOptimization::EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry, CornerData<double> beta) : mesh(mesh), geometry(geometry), beta(beta) {
    // Initialize quantities
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

// convenience function for making a sparse matrix
monty::rc_ptr<mosek::fusion::Matrix> EmbeddingOptimization::sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values ) {
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto res =  Matrix::sparse(m,n,r,c,v);
    return res;
}

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
    // To do so, we associate each quad obtained in the first step with it's associated
    // corner, indexing within this quad as follows:
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
    //  In order to assemble our final quad mesh though, we need to identify vertices
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
                //debugging code
                /*
                   if(positions[finalIndices[cOffset + x + n * y]].norm() != 0) {
                   double err = (positions[finalIndices[cOffset + x + n * y]] 
                   - bary(c, x, y)).norm();
                   if (err > 1e-6) {
                   cout << "reindex error:" << err << " from " 
                   << cOffset + x + n * y << endl;
                   }
                   } 
                   */
                positions[finalIndices[cOffset + x + n * y]] = bary(c, x, y);
            }
        }
    }
    subgeometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*submesh,positions));

    // TODO: change this call to solve()
    polyscope::registerSurfaceMesh("New mesh", positions, submesh->getFaceVertexList());
    return;
}


// writes the expected edge lengths/angles associated to each quad to c_iso_i for
// i = 0,1,2. 
// c_iso_0 = |ik|^2, 
// c_iso_1 = |jl|^2, 
// c_iso_2 = < ki, lj >
//
//    l --------k
//   /          /
//  /          /
// i -------- j

void EmbeddingOptimization::buildIntrinsicCheckerboard(){
    // initialize vectors.
    // Note that these are lenghts per quad, hence the count is different from the number of vertices
    c_iso_0 = vector<double>(nCorners * (n-1) * (n-1));
    c_iso_1 = vector<double>(nCorners * (n-1) * (n-1));
    c_iso_2 = vector<double>(nCorners * (n-1) * (n-1));
    VertexData<size_t> vMap = submesh->getVertexIndices();
    for (Corner c: mesh->corners()) {
        // see above
        size_t cQuadOffset = c_[c] * (n-1) * (n-1);
        // grab face coordinates
        Halfedge IJ = c.halfedge();
        Vector3 I = geometry->inputVertexPositions[IJ.vertex()];
        Halfedge JK = IJ.next();
        Vector3 J = geometry->inputVertexPositions[JK.vertex()];
        Halfedge KI = JK.next();
        Vector3 K = geometry->inputVertexPositions[KI.vertex()];

        // solve for coefficients associated to input data
        // DEBUG
        double bij = beta[IJ.corner()];
        double bjk = beta[JK.corner()];
        double bki = beta[KI.corner()];
        BezierTriangle T = Coefficients(I,J,K,beta[IJ.corner()],beta[JK.corner()],beta[KI.corner()]);
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                size_t index = cQuadOffset + x + (n-1)*y;
                Vector2 i = RationalBezierTriangle(T,baryCoords(x, y));
                Vector2 j = RationalBezierTriangle(T,baryCoords(x+1, y));
                Vector2 k = RationalBezierTriangle(T,baryCoords(x+1, y+1));
                Vector2 l = RationalBezierTriangle(T,baryCoords(x, y+1));
                Vector2 e20 = i - k;
                Vector2 e31 = j - l;
                c_iso_0[index] = e20.norm2();
                c_iso_1[index] = e31.norm2();
                c_iso_2[index] = dot(e20,e31);
            }
        }
    }
}

inline double sqr(double x) { return x * x; }

void EmbeddingOptimization::evaluateEnergy(double& energy, const Eigen::VectorXd& v){
    // For each quad ijkl, the energy is
    //    l --------k
    //   /          /
    //  /          /
    // i -------- j
    //  
    // c_iso_0 = |ik|^2, 
    // c_iso_1 = |jl|^2, 
    // c_iso_2 = < ki, lj >

    energy = 0.;
    for (Corner c: mesh->corners()) {
        // offset for the 3 stored intrinsic lengths
        size_t cQuadOffset = c_[c] * (n-1) * (n-1);
        size_t cVertexOffset = c_[c] * n * n;
        // grab face coordinates
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                size_t iIndex = finalIndices[cVertexOffset + x + (n)*y];
                Vector3 i = {v[3*iIndex], v[3*iIndex+1],v[3*iIndex+2]};

                size_t jIndex = finalIndices[cVertexOffset + (x+1) + (n)*y];
                Vector3 j = {v[3*jIndex], v[3*jIndex+1],v[3*jIndex+2]};

                size_t kIndex = finalIndices[cVertexOffset + (x+1) + n*(y+1)];
                Vector3 k = {v[3*kIndex], v[3*kIndex+1],v[3*kIndex+2]};

                size_t lIndex = finalIndices[cVertexOffset + x + n*(y+1)];
                Vector3 l = {v[3*lIndex], v[3*lIndex+1],v[3*lIndex+2]};

                size_t quadIndex = cQuadOffset + x + (n-1)*y;

                Vector3 e20 = i - k;
                // c_iso_0 term
                energy += (sqr(e20.norm2() - c_iso_0[quadIndex]))/2;
                //DEBUG
                //cout << "i : " << i << " k: " << k << "norm: " << e20.norm() << "should: "<< sqrt(c_iso_0[quadIndex]) << endl;
                Vector3 e31 = j - l;

                //cout << "actual: " <<  e20.norm2() << " expected:  " << c_iso_0[quadIndex] << endl;

                // c_iso_1 term
                energy += (sqr(e31.norm2() - c_iso_1[quadIndex]))/2;
                // c_iso_2 term
                energy += (sqr(dot(e20,e31) - c_iso_2[quadIndex]))/2;
            }
        }
        // regularization term
    }
}

void EmbeddingOptimization::evaluateGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& v) {
    for (Corner c: mesh->corners()) {
        for (int i = 0; i < 3 * nSubdividedVertices; i++) {
            gradient[i] = 0.;
        }
        size_t cQuadOffset = c_[c] * (n-1) * (n-1);
        size_t cVertexOffset = c_[c] * n * n;
        // grab face coordinates
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                size_t iIndex = finalIndices[cVertexOffset + x + (n)*y];
                Vector3 i = {v[3*iIndex], v[3*iIndex+1],v[3*iIndex+2]};

                size_t jIndex = finalIndices[cVertexOffset + (x+1) + (n)*y];
                Vector3 j = {v[3*jIndex], v[3*jIndex+1],v[3*jIndex+2]};

                size_t kIndex = finalIndices[cQuadOffset + (x+1) + n*(y+1)];
                Vector3 k = {v[3*kIndex], v[3*kIndex+1],v[3*kIndex+2]};

                size_t lIndex = finalIndices[cQuadOffset + x + n*(y+1)];
                Vector3 l = {v[3*lIndex], v[3*lIndex+1],v[3*lIndex+2]};

                size_t quadIndex = cQuadOffset + x + (n-1)*y;

                Vector3 e20 = i - k;
                Vector3 e31 = j - l;
                //================== c_iso_0 gradients===================
                double c_iso_0_weight = (e20.norm2() - c_iso_0[quadIndex]);
                // partials with respect to i are
                // (|e20|^2 - |e20'|^2|) * 2 * (i0 - k0)
                // i gradients
                gradient[3*iIndex] += c_iso_0_weight * 2 *(i.x - k.x);
                gradient[3*iIndex+1] += c_iso_0_weight * 2 *(i.y - k.y);
                gradient[3*iIndex+2] += c_iso_0_weight * 2 *(i.z - k.z);

                // k gradients
                gradient[3*kIndex] += c_iso_0_weight * 2 *(k.x - i.x);
                gradient[3*kIndex+1] += c_iso_0_weight * 2 *(k.y - i.y);
                gradient[3*kIndex+2] += c_iso_0_weight * 2 *(k.z - i.z);


                //===================== c_iso_1 gradients===========================
                double c_iso_1_weight = (e31.norm2() - c_iso_1[quadIndex]);
                // j gradients
                gradient[3*jIndex] += c_iso_1_weight * 2 *(j.x - l.x);
                gradient[3*jIndex+1] += c_iso_1_weight * 2 *(j.y - l.y);
                gradient[3*jIndex+2] += c_iso_1_weight * 2 *(j.z - l.z);

                // l gradients
                gradient[3*lIndex] += c_iso_1_weight * 2 *(l.x - j.x);
                gradient[3*lIndex+1] += c_iso_1_weight * 2 *(l.y - j.y);
                gradient[3*lIndex+2] += c_iso_1_weight * 2 *(l.z - j.z);

                //===================== c_iso_2 gradients===========================
                double c_iso_2_weight = (dot(e20,e31) - c_iso_2[quadIndex]);

                // i partial: j - l
                gradient[3*iIndex] += c_iso_2_weight *(j.x - l.x);
                gradient[3*iIndex+1] += c_iso_2_weight *(j.y - l.y);
                gradient[3*iIndex+2] += c_iso_2_weight *(j.z - l.z);

                // k partial: l - j
                gradient[3*kIndex] += c_iso_2_weight * (l.x - j.x);
                gradient[3*kIndex+1] += c_iso_2_weight * (l.y - j.y);
                gradient[3*kIndex+2] += c_iso_2_weight * (l.z - j.z);

                // j partial: i - k
                gradient[3*jIndex] += c_iso_2_weight * (i.x - k.x);
                gradient[3*jIndex+1] += c_iso_2_weight * (i.y - k.y);
                gradient[3*jIndex+2] += c_iso_2_weight * (i.z - k.z);

                // l partial: k - i
                gradient[3*lIndex] += c_iso_2_weight * (k.x - i.x);
                gradient[3*lIndex+1] += c_iso_2_weight * (k.y - i.y);
                gradient[3*lIndex+2] += c_iso_2_weight * (k.z - i.z);
            }
        }
    }
}

Eigen::VectorXd EmbeddingOptimization::gradientDescent() {
    double BETA = 0.5;
    double EPSILON = 1e-7;
    double ALPHA = 0.25;
    int MAX_ITERS=1000;
    int k = 1;
    Eigen::VectorXd v = Eigen::VectorXd::Zero(3*nSubdividedVertices);
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int X = 0; X < n; X++) {
            for (int Y = 0; Y < n; Y++) {
                Vector3 pos = bary(c, X, Y);
                size_t startIndex = 3*finalIndices[cOffset + X + n * Y];
                v[startIndex] = pos.x;
                v[startIndex + 1] = pos.y;
                v[startIndex + 2] = pos.z;
            }
        }
    }
    Eigen::VectorXd v_new = Eigen::VectorXd::Zero(3*nSubdividedVertices);

    while (true) {
        double energy = 0.0;
        evaluateEnergy(energy, v);

        // compute update direction
        double energy_new = 0.;
        Eigen::VectorXd g(3*nSubdividedVertices);
        evaluateGradient(g, v);
        double grad_size = sqrt(g.dot(g));

        // compute step size
        double t = 1.;
        v_new = v - t * g;
        evaluateEnergy(energy_new, v_new);
        while (energy_new > energy - ALPHA * t * grad_size) {
            t = BETA*t;
            v_new = v - t * g;
            evaluateEnergy(energy_new, v_new);
        }

        // update
        k++;

        // check termination condition
        cout << "ENERGY: " << energy << endl;
        cout << grad_size << endl;
        if (grad_size < EPSILON || k > MAX_ITERS) break;
        if (k % 100 == 0) cout << "iters:" << k << endl;
    }

    std::cout << "Iterations: " << k << std::endl;
    return v;
}

// check that our derivative code agrees with numerical
// derivatives, up to a small error
void EmbeddingOptimization::testFlatteningDerivatives() {
    // perturbation size
    const double eps = 1e-7;

    // grab current configuration
    Eigen::VectorXd x = Eigen::VectorXd::Zero(3*nSubdividedVertices);
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int X = 0; X < n; X++) {
            for (int Y = 0; Y < n; Y++) {
                Vector3 pos = bary(c, X, Y);
                size_t startIndex = 3*finalIndices[cOffset + X + n * Y];
                x[startIndex] = pos.x;
                x[startIndex + 1] = pos.y;
                x[startIndex + 2] = pos.z;
            }
        }
    }

    // compute energy and its derivatives
    double E0 = 0.;
    evaluateEnergy(E0, x);
    Eigen::VectorXd dE(3*nSubdividedVertices); // differential of energy with respect to the variable
    evaluateGradient(dE, x);

    double worstError = 0.;
    double averageError = 0.;

    for (int i = 0; i < 3 * nSubdividedVertices; i++) {
        double dEi;  // numerical derivative with respect to index i

        x[i] += eps;  // perturb coordinate k
        double E = 0.;
        evaluateEnergy(E,x);
        dEi = (E - E0) / eps;
        x[i] -= eps;  // restore it

        // compare to expression computed in code
        double error = fabs(dEi - dE[i])/(1 + dEi);
        if (error > 1) {
            cout << dEi << " " << dE[i] << " error: " << error << endl;
            cout << "i = " << i/3 << endl;
        }

        worstError = std::max(worstError, error);
        averageError += error;
    }

    averageError /= (3 * nSubdividedVertices);

    std::cout << "Worst error: " << worstError << std::endl;
    std::cout << "Average error: " << averageError << std::endl;
}

void EmbeddingOptimization::basisFunctionDebugging() {
    double N = 2*n*(n-1);
    // ==================================TEST============================================
    std::ofstream ss("test.svg", std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
          "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
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
std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry> >EmbeddingOptimization::solve(int N) {
    n = N;
    //===============DEBUG===================
    //basisFunctionDebugging();
    //testFlatteningDerivatives();
    //return {submesh, subgeometry};
    //=======================================

    // Initialize union find data structure
    top = vector<int> (n*n*nCorners);
    for (int i = 0; i < n * n * nCorners; i++) top[i] = i;
    next = vector<int> (n*n*nCorners, -1);

    buildEquivalenceClasses();
    nSubdividedVertices = buildFinalIndices();

    // Initialize submesh and subgeometry
    buildSubdivision();
    buildIntrinsicCheckerboard();



    Eigen::VectorXd x = gradientDescent();

    VertexData<Vector3> positions(*submesh);
    for (int i = 0; i < nSubdividedVertices; i++) {
        positions[i] = {x[3*i],x[3*i+1],x[3*i+2]};
    }

    // TODO: why doesn't this work?
    polyscope::registerSurfaceMesh("Mesh", positions, submesh->getFaceVertexList());
    polyscope::show();
    /*
    */

    return {submesh, subgeometry};
}

