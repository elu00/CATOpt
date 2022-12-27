#include "PrescribeCurvature.h"
// Input: A CAT surface with edge lengths l : E -> R>0
// and interior angles 𝛽 : C -> R, together with target Gaussian curvatures 
// Ωi, Ωij at each interior vertex i and edge ij
// (resp.), and target geodesic curvatures k𝑖, k𝑖𝑗 at each boundary vertex
// i and edge ij (resp.). Curvatures must satisfy the CAT Gauss-Bonnet condition from Equation 10.
// Output: A CAT surface (T, \tilde l, \tilde 𝛽) with the prescribed curvatures. For 
// convenience, we also return the corresponding edge bend angles \tilde \alpha
void PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l, CornerData<double> beta, 
        VertexData<double> VertexCurvatures, EdgeData<double> EdgeCurvatures) {

    VertexData<double> VertexCurvatureStar(*mesh);
    for (Vertex i: mesh->vertices()) {
        if (i.isBoundary()) {
            // get the two boundary halfedges into and out of i
            Edge e1, e2;
            for (Halfedge h: i.outgoingHalfedges()) {
                if (!h.isInterior()) {
                    e1 = h.edge();
                    break;
                }
            }
            for (Halfedge h: i.incomingHalfedges()) {
                if (!h.isInterior()) {
                    e2 = h.edge();
                    break;
                }
            }
            assert(e1 != e2);
            VertexCurvatureStar[i] = VertexCurvatures[i] + (EdgeCurvatures[e1] + EdgeCurvatures[e2])/2.;
        } else {
            // calculate (Σ_ij Ωij)/2
            double accum = 0.;
            for (Edge ij: i.adjacentEdges()) {
                accum += EdgeCurvatures[ij]/2.;
            }
            VertexCurvatureStar[i] = VertexCurvatures[i] + accum;
        }
    }
    // TODO: add intrinsic flattening stuff in
    /*
    IntrinsicFlattening flattening(mesh); 
    // arguments are VertexData targetCurvatures, CornerData targetBetas
    // output is pair
    auto [A,B] = flattening.CoherentAngleSystem(
    */

    // calculate intersection angles from CAS
    EdgeData<double> eta(*mesh);

    //EdgeData<double> lStar = MinimizeLengthEnergy(mesh);
}
void LayoutMesh() {

}
