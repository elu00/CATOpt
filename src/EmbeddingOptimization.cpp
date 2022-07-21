#include "fusion.h"
#include "EmbeddingOptimization.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

using namespace mosek::fusion;
using namespace monty;

Vector3 EmbeddingOptimization::bary(Corner c, int X, int Y) {
    double x = X;
    double y = Y;
   auto it = c.halfedge();
   Vector3 i = geometry->inputVertexPositions[it.vertex()];
   it = it.next();
   Vector3 j = geometry->inputVertexPositions[it.vertex()];
   it = it.next();
   Vector3 k = geometry->inputVertexPositions[it.vertex()];
   double jWeight = x * (3*n-3-2*y) * (3*n-3-y)/(2*(n-1)*(9-18*n + 9 * n * n + 6*y -6*n*y-x*y+y*y));
   double kWeight = y * (3*n-3-2*x) * (3*n-3-x)/(2*(n-1)*(9-18*n + 9 * n * n + 6*x -6*n*x-x*y+x*x));
   return (1 - jWeight - kWeight) * i + jWeight * j + kWeight * k;
   }
EmbeddingOptimization::EmbeddingOptimization(shared_ptr<ManifoldSurfaceMesh> mesh, shared_ptr<VertexPositionGeometry> geometry):
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
            
    }


// convenience function for making a sparse matrix
monty::rc_ptr<mosek::fusion::Matrix> EmbeddingOptimization::sMatrix(int m, int n, vector<int>& rows, vector<int>& cols, vector<double>& values ) {
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto res =  Matrix::sparse(m,n,r,c,v);
    rows.clear();
    cols.clear();
    values.clear();
    return res;
}

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
int EmbeddingOptimization::find(int a) {
    return top[a];
}
 std::pair<shared_ptr<ManifoldSurfaceMesh>, shared_ptr<VertexPositionGeometry> >EmbeddingOptimization::solve(int N) {
    n = N;
    top = vector<int> (n*n*nCorners);
    for (int i = 0; i < n * n * nCorners; i++) top[i] = i;
    next = vector<int> (n*n*nCorners, -1);

    /*
     *        n*(n-1) -> .. -> n^2 - 1
     *        /                    /
     *      .........................
     *     /                      /
     *    n               2*n - 1 <-> n^2-n+1
     *   /                      /
     *   0 -> 1 -> ... -> n - 1  <-> n^2 - n
     *   ^    ^
     *   |    |
     *   v    v
     *   0    n
     */
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
    for (int i = 0; i < n * n * nCorners; i++) cout << find(i) << " ";
    cout << endl;

    // reindexing so all our final vertex indices are contiguous
    std::map<int,int> reindex;
    size_t temp = 0;
    for (int i = 0; i < nCorners * n * n; i++) {
        if (!reindex.count(find(i))){
            reindex[find(i)] = temp;
            temp++;
        }
        finalIndices.push_back(reindex[find(i)]);
    }

    for (int i = 0; i < n * n * nCorners; i++) cout << finalIndices[i] << " ";
    cout << endl;

    // now all the indexing garbage is done...let's construct the mesh and setup the energy
    vector< vector<size_t> > polygons;
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n - 1; x++) {
            for (int y = 0; y < n - 1; y++) {
                size_t i = finalIndices[find(cOffset + x + n * y)];
                size_t j = finalIndices[find(cOffset + (x+1) + n * y)];
                size_t k = finalIndices[find(cOffset + (x+1) + n * (y+1))];
                size_t l = finalIndices[find(cOffset + x + n * (y+1))];
                polygons.push_back({i,j,k,l});
                cout << i << " " << j << " " << k << " " << l << "wat" << endl;
            }
        }
    }
    submesh = std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(polygons));
    VertexData<Vector3> positions(*submesh);
    for (Corner c: mesh->corners()) {
        size_t cOffset = c_[c] * n * n;
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                positions[finalIndices[find(cOffset + x + n * y)]] = bary(c, x, y);
                //cout << (bary(c,x,y));
                //cout << find(cOffset + x + n * y) << endl;
            }
        }
    }
    for (Vertex v: submesh->vertices()) cout << positions[v] << "wah" << endl;
    subgeometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*submesh,positions));
    polyscope::registerSurfaceMesh("New mesh", positions, submesh->getFaceVertexList());
    polyscope::show();
    return {submesh, subgeometry};
    }
