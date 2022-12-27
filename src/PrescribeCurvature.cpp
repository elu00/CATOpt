#include "PrescribeCurvature.h"
// Input: A CAT surface with edge lengths l : E -> R>0
// and interior angles 𝛽 : C -> R, together with target Gaussian curvatures 
// Ωi, Ωij at each interior vertex i and edge ij
// (resp.), and target geodesic curvatures k𝑖, k𝑖𝑗 at each boundary vertex
// i and edge ij (resp.). Curvatures must satisfy the CAT Gauss-Bonnet condition from Equation 10.
// Output: A CAT surface (T, \tilde l, \tilde 𝛽) with the prescribed curvatures. For 
// convenience, we also return the corresponding edge bend angles \tilde \alpha
pair<EdgeData<double>, CornerData<double>> PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, 
        EdgeData<double> l, CornerData<double> beta, 
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
    IntrinsicFlattening flattening(mesh, l); 
    // arguments are VertexData targetCurvatures, CornerData targetBetas
    auto [thetaHat, betaHat] = flattening.CoherentAngleSystem(VertexCurvatureStar, beta);

    // calculate intersection angles from CAS
    EdgeData<double> eta(*mesh);
    for (Edge e: mesh->edges()) {
        if (e.isBoundary()) {
            // calculate which of the halfedges is interior
            Halfedge ij = e.halfedge().isInterior() ? e.halfedge() : e.halfedge().twin();
            eta[e] = M_PI - thetaHat[ij.next().next()];
        } else {
            Halfedge ij = e.halfedge();
            Halfedge ji = e.halfedge().twin();
            eta[e] = M_PI - thetaHat[ij.next().next()] - thetaHat[ji.next().next()];
        }
    }

    EdgeLengthOptimization E(mesh, l, eta);
    Eigen::VectorXd temp = E->MinimizeLengthEnergy();
    EdgeData<double> lHat(*mesh);
    EdgeData<size_t> eIndices = mesh->getEdgeIndices();
    for (Edge e: mesh->edges()) {
        lHat[e] = temp[eIndices[e]];
    }

    // TODO: double check why lines 7-10 of the pseudocode are there...
    return {lHat, betaHat};
}
void LayoutMesh() {

}
