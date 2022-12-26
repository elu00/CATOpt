#include "PrescribeCurvature.h"
// Input: A CAT surface with edge lengths l : E -> R>0
// and interior angles 𝛽 : C -> R, together with target Gaussian curvatures 
// Ω𝑖, Ω𝑖𝑗 at each interior vertex i and edge ij
// (resp.), and target geodesic curvatures k𝑖, k𝑖𝑗 at each boundary vertex
// i and edge ij (resp.). Curvatures must satisfy the CAT Gauss-Bonnet condition from Equation 10.
// Output: A CAT surface (T, \tilde l, \tilde 𝛽) with the prescribed curvatures. For 
// convenience, we also return the corresponding edge bend angles \tilde \alpha
void PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l, CornerData<double> beta, 
        VertexData<double> VertexCurvatures, EdgeData<double> EdgeCurvatures) {
    VertexData<double> VertexCurvatureStar(*mesh);
    for (Vertex i: mesh->vertices()) {
        if (i.isBoundary()) {
        }
    }
    // calculate intersection angles from CAS
    EdgeData<double> eta(*mesh);
    // calculate 
    //EdgeData<double> lStar = MinimizeLengthEnergy(mesh);
}
void LayoutMesh() {

}
