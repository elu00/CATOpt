#include "PrescribeCurvature.h"
// Input: A CAT surface with edge lengths l : E -> R>0
// and interior angles 𝛽 : C -> R, together with target Gaussian curvatures
// Ωi, Ωij at each interior vertex i and edge ij
// (resp.), and target geodesic curvatures ki, kij at each boundary vertex
// i and edge ij (resp.). Curvatures must satisfy the CAT Gauss-Bonnet condition
// from Equation 10.
// Output: A CAT surface (T, \tilde l, \tilde 𝛽) with the
// prescribed curvatures. For convenience, we also return the corresponding edge
// bend angles \tilde \alpha
std::pair<EdgeData<double>, CornerData<double>>
PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l,
                   CornerData<double> beta, VertexData<double> VertexCurvatures,
                   EdgeData<double> EdgeCurvatures) {

    VertexData<double> VertexCurvatureStar(*mesh);
    for (Vertex i : mesh->vertices()) {
        if (i.isBoundary()) {
            // get the two boundary halfedges into and out of i
            Edge e1, e2;
            for (Halfedge h : i.outgoingHalfedges()) {
                if (!h.isInterior()) {
                    e1 = h.edge();
                    break;
                }
            }
            for (Halfedge h : i.incomingHalfedges()) {
                if (!h.isInterior()) {
                    e2 = h.edge();
                    break;
                }
            }
            assert(e1 != e2);
            VertexCurvatureStar[i] =
                VertexCurvatures[i] +
                (EdgeCurvatures[e1] + EdgeCurvatures[e2]) / 2.;
        } else {
            // calculate (Σ_ij Ωij)/2
            double accum = 0.;
            for (Edge ij : i.adjacentEdges()) {
                accum += EdgeCurvatures[ij] / 2.;
            }
            VertexCurvatureStar[i] = VertexCurvatures[i] + accum;
        }
    }
    IntrinsicFlattening flattening(mesh, l);
    // arguments are VertexData targetCurvatures, CornerData targetBetas
    auto [thetaHat, betaHat] =
        flattening.CoherentAngleSystem(VertexCurvatureStar, beta);

    // calculate intersection angles from CAS
    EdgeData<double> eta(*mesh);
    for (Edge e : mesh->edges()) {
        if (e.isBoundary()) {
            // calculate which of the halfedges is interior
            Halfedge ij =
                e.halfedge().isInterior() ? e.halfedge() : e.halfedge().twin();
            eta[e] = M_PI - thetaHat[ij.next().next().corner()];
        } else {
            Halfedge ij = e.halfedge();
            Halfedge ji = e.halfedge().twin();
            eta[e] = M_PI - thetaHat[ij.next().next().corner()] -
                     thetaHat[ji.next().next().corner()];
        }
    }

    EdgeLengthOptimization E(mesh, l, eta);
    Eigen::VectorXd temp = E.MinimizeLengthEnergy();
    EdgeData<double> lHat(*mesh);
    EdgeData<size_t> eIndices = mesh->getEdgeIndices();
    for (Edge e : mesh->edges()) {
        lHat[e] = temp[eIndices[e]];
    }

    // TODO: double check why lines 7-10 of the pseudocode are there...
    return {lHat, betaHat};
}

Vector2 Rotate(Vector2 u, double phi) {
    return {cos(phi) * u.x - sin(phi) * u.y, sin(phi) * u.x + cos(phi) * u.y};
}
Vector2 FinishTriangle(Vector2 x_i, Vector2 x_j, double l_ij, double l_jk,
                       double l_ki) {
    double theta_ijk =
        acos((l_ij * l_ij - l_jk * l_jk + l_ki * l_ki) / (2. * l_ij * l_ki));
    return x_i + Rotate(x_j - x_i, theta_ijk) * l_ki / l_ij;
}

inline std::tuple<Corner, Corner, Corner> Corners(Halfedge h) {
    return {h.corner(), h.next().corner(), h.next().next().corner()};
}
CornerData<Vector2> LayoutMesh(shared_ptr<ManifoldSurfaceMesh> mesh,
                               EdgeData<double> l, EdgeData<bool> S) {
    CornerData<Vector2> x(*mesh);

    FaceData<bool> visited(*mesh, false);
    std::queue<Face> Q;
    Face f = mesh->face(0);
    // assign coordinates to face f
    Halfedge ij = f.halfedge();
    auto [i, j, k] = Corners(ij);
    x[i] = {0, 0};
    x[j] = {l[ij.edge()], 0};
    x[k] = FinishTriangle(x[i], x[j], l[ij.edge()], l[ij.next().edge()],
                          l[ij.next().next().edge()]);
    Q.push(f);
    while (!Q.empty()) {
        Face ijk = Q.front();
        Q.pop();
        // grab adjacent faces
        for (Halfedge ij : ijk.adjacentHalfedges()) {
            Halfedge ji = ij.twin();
            Face jil = ji.face();
            if (!visited[jil] && !S[ij.edge()]) {

                auto [i_jk, j_ki, k_ij] = Corners(ij);
                auto [j_il, i_lj, l_ji] = Corners(ji);

                x[i_lj] = x[i_jk];
                x[j_il] = x[j_ki];
                x[l_ji] = FinishTriangle(x[j_il], x[i_lj], l[ji.edge()],
                                         l[ji.next().edge()],
                                         l[ji.next().next().edge()]);

                visited[jil] = true;
                Q.push(jil);
            }
        }
    }
    return x;
}

// Input: A CAT surface with corner positions p: C -> R^2 and edge
// and interior angles 𝛽 : C -> R
// Output: an SVG with the given corner positions and edge bend angles written
// to the specified filename
void CATToSVG(shared_ptr<ManifoldSurfaceMesh> mesh, CornerData<Vector2> p,
              CornerData<double> beta, std::string filename) {

    CornerData<Vector2> offsetPositions(*mesh);

    double minX = p[mesh->corner(0)].x;
    double maxX = p[mesh->corner(0)].x;
    double minY = p[mesh->corner(0)].y;
    double maxY = p[mesh->corner(0)].y;

    // calculate min/max of corner positions
    for (Corner c : mesh->corners()) {
        minX = std::min(minX, p[c].x);
        maxX = std::max(maxX, p[c].x);
        minY = std::min(minY, p[c].y);
        maxY = std::max(maxY, p[c].y);
    }
    double deltaX = maxX - minX;
    double deltaY = maxY - minY;
    // int width = ceil(deltaX);
    // int height = ceil(deltaY);
    int width = 100, height = 100;

    // do renormalization
    for (Corner c : mesh->corners()) {
        Vector2 pos = p[c];
        offsetPositions[c] = {(pos.x - minX) / deltaX * width,
                              (pos.y - minY) / deltaY * height};
    }
    // initialize file I/O
    std::ofstream ss(filename, std::ofstream::out);

    // SVG Header
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
          "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
       << endl
       << "<svg width=\"" << width << "\" height=\"" << height
       << "2000\" "
          "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >"
       << "<rect width=\"100%\" height =\"100%\" fill=\"#ffffff\"/>" << endl;

    // point labels
    for (Corner c : mesh->corners()) {
        ss << "<circle cx=\"" << offsetPositions[c].x << "\" cy=\""
           << offsetPositions[c].y << "\" r=\"1\"/>" << endl;
    }

    // circular arcs
    for (Halfedge h : mesh->halfedges()) {
        Vector2 i = offsetPositions[h.corner()];
        Vector2 j = offsetPositions[h.next().corner()];
        // calculate arc length/radius
        // DEBUG: fix this
        double angle = 0;
        // double angle = beta[h.edge()];
        double radius = (i - j).norm() / abs(2 * sin(angle));
        if (abs(angle) > 1e-7) {
            // this flag is 0 if the drawn arc should be the "minor" arc
            // and 1 if it should be the major arc
            std::string largeArcFlag = std::abs(2 * angle) <= M_PI ? "0" : "1";
            // sweep flag is 1 if going outward, 0 if going inward
            // visually check this at some point
            std::string sweepFlag = angle >= 0 ? "1" : "0";
            ss << "<path d=\"M" << i.x << "," << i.y << " A" << radius << ","
               << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
               << j.x << "," << j.y
               << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />"
               << endl;
        } else {
            // flat line case
            ss << "<line x1=\"" << i.x << "\" x2=\"" << j.x << "\" y1=\"" << i.y
               << "\" y2=\"" << j.y << "\" stroke=\"blue\" stroke-width=\"2\"/>"
               << endl;
        }
    }
    // footer
    ss << "</svg>";
}
