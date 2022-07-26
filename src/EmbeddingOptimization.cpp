#include "fusion.h"
#include "EmbeddingOptimization.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace mosek::fusion;
using namespace monty;

// Given a corner c and integer weights X,Y representing local coordinates 
// on c's quad, returns the associated point in R^3.
// Here the integers X,Y are assumed to be in the range
// [0,N-1], with the mapping to barycentric
// coordinates being given by the transformation induced by
// (0,0)->(1,0,0), (N-1,N-1)->(1/3,1/3,1/3), (0,N-1)->(1/2,0,1/2)
// and linear interpolation in between
Vector3 EmbeddingOptimization::bary(Corner c, int X, int Y)
{
    // Normalize x and y
    double x = (double)X/(n-1);
    double y = (double)Y/(n-1);

    // grab face coordinates
    auto it = c.halfedge();
    Vector3 i = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 j = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 k = geometry->inputVertexPositions[it.vertex()];


    double jWeight = -(y-3) * x/6;
    double kWeight = -(x-3) * y/6;
    return (1 - jWeight - kWeight) * i + jWeight * j + kWeight * k;
}
EmbeddingOptimization::EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry) : mesh(mesh), geometry(geometry) {
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
                size_t i = finalIndices[find(cOffset + x + n * y)];
                size_t j = finalIndices[find(cOffset + (x+1) + n * y)];
                size_t k = finalIndices[find(cOffset + (x+1) + n * (y+1))];
                size_t l = finalIndices[find(cOffset + x + n * (y+1))];
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
                //debuging code
                /*
                if(positions[finalIndices[find(cOffset + x + n * y)]].norm() != 0) {
                    double err = (positions[finalIndices[find(cOffset + x + n * y)]] 
                                - bary(c, x, y)).norm();
                    if (err > 1e-6) {
                        cout << "reindex error:" << err << " from " 
                            << cOffset + x + n * y << endl;
                    }
                } 
                */
                positions[finalIndices[find(cOffset + x + n * y)]] = bary(c, x, y);
            }
        }
    }
    subgeometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*submesh,positions));
    return;
}


void EmbeddingOptimization::evaluateEnergy(double& energy, const Eigen::VectorXd& v){
    // For each quad ijkl, the energy is
    //    l --------k
    //   /          /
    //  /          /
    // i -------- j
    //
}

void EmbeddingOptimization::evaluateGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& v) {
}

Eigen::VectorXd EmbeddingOptimization::gradientDescent() {
    double beta = 0.5;
    double EPSILON = 1e-7;
    int MAX_ITERS=10000;
    int k = 1;
    double f = 0.0, tp = 1.0;
    vector<double> obj;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(3*n);
    evaluateEnergy(f, x);
    obj.push_back(f);
    Eigen::VectorXd xp = Eigen::VectorXd::Zero(3*n);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(3*n);

    while (true) {
        // compute momentum term
        v = x;
        if (k > 1) v += (k-2)*(x - xp)/(k+1);
        
        // compute update direction
        Eigen::VectorXd g(3*n);
        evaluateGradient(g, v);
        
        // compute step size
        double t = tp;
        double fp = 0.0;
        Eigen::VectorXd xn = v - t*g;
        Eigen::VectorXd xnv = xn - v;
        evaluateEnergy(fp, v);
        evaluateEnergy(f, xn);
        while (f > fp + g.dot(xnv) + xnv.dot(xnv)/(2*t)) {
            t = beta*t;
            xn = v - t*g;
            xnv = xn - v;
            evaluateEnergy(f, xn);
        }
    
        // update
        tp = t;
        xp = x;
        x  = xn;
        obj.push_back(f);
        k++;

        // check termination condition
        if (fabs(f - fp) < EPSILON || k > MAX_ITERS) break;
    }

    std::cout << "Iterations: " << k << std::endl;
    return x;
}

std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry> >EmbeddingOptimization::solve(int N) {
    n = N;
    // Initialize union find data structure
    top = vector<int> (n*n*nCorners);
    for (int i = 0; i < n * n * nCorners; i++) top[i] = i;
    next = vector<int> (n*n*nCorners, -1);

    buildEquivalenceClasses();
    buildFinalIndices();

    // Initialize submesh and subgeometry
    buildSubdivision();

    polyscope::registerSurfaceMesh("New mesh", subgeometry->vertexPositions, submesh->getFaceVertexList());
    polyscope::show();

    return {submesh, subgeometry};
}
