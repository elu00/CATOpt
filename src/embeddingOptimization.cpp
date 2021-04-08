#include "CatOpt.h"
#include "CAT.h"
Vector3 CatOpt::bary(Face f, double a, double b, double c) {
    auto it = f.halfedge();
    Vector3 i = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 j = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 k = geometry->inputVertexPositions[it.vertex()];
    return a * i + b * j + c * k;
}

// initializes CATMesh, handles all subdivision energy initialization
void CatOpt::subdivision() {
    EdgeData<size_t> eStart(*mesh, -1);
    VertexData<size_t> v(*mesh, -1);
    FaceData<std::map<size_t,size_t>> indexing(*mesh);
    FaceData<CAT> CATs(*mesh);
    FaceData<vector<double>> f_points(*mesh);
    // initializing indexing map because I'm lazy
    cout << "starting subdivision" << endl;
    size_t index = 0;
    vector<vector<size_t>> polygons;
    // initialization of initial guesses
    for (Face f : mesh->faces()) {
        double ij, jk, ki, a_ij, a_jk, a_ki;
        Halfedge h = f.halfedge();
        ij = geometry->edgeLength(h.edge());
        a_ij = (sol)[eInd[h.edge()]];
        h = h.next();
        jk = geometry->edgeLength(h.edge());
        a_jk = (sol)[eInd[h.edge()]];
        h = h.next();
        ki = geometry->edgeLength(h.edge());
        a_ki = (sol)[eInd[h.edge()]];
        CAT triangle = CAT(ij, jk, ki, a_ij, a_jk, a_ki);
        vector<double> finalPoints;
        vector<int> triangles;
        std::tie(finalPoints, triangles) = triangle.triangulation(subdiv_level);
        std::map<size_t, size_t> toFin;
        CATs[f] = triangle;
        f_points[f] = finalPoints;


        auto it = f.halfedge();
        Vertex i = it.vertex();
        if (v[i] == -1) {
            v[i] = index;
            toFin[0] = index;
            index++;
        } else {
            toFin[0] = v[i];
        }
        it = it.next();

        Vertex j = it.vertex();
        if (v[j] == -1) {
            v[j] = index;
            toFin[subdiv_level] = index;
            index++;
        } else {
            toFin[subdiv_level] = v[j];
        }
        it = it.next();

        Vertex k = it.vertex();
        if (v[k] == -1) {
            v[k] = index;
            toFin[2*subdiv_level] = index;
            index++;
        } else {
            toFin[2*subdiv_level] = v[k];
        }
        it = it.next();

        // edge ij
        size_t start;
        //which_e = eInd[it.edge()];
        if (eStart[it.edge()] == -1) {
            eStart[it.edge()] = index;
            start = index;
            index += subdiv_level - 1;
        } else {
            start = eStart[it.edge()];
        }
        if (it.edge().halfedge() == it) {
            // i is canonical vertex
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[pos] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[subdiv_level - pos] = start;
                start++;
            }
        }
        it = it.next();

        // edge jk
        //which_e = eInd[it.edge()];
        if (eStart[it.edge()] == -1) {
            eStart[it.edge()] = index;
            start = index;
            index += subdiv_level - 1;
        } else {
            start = eStart[it.edge()];
        }
        if (it.edge().halfedge() == it) {
            // j is canonical vertex
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[subdiv_level + pos] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[2*subdiv_level - pos] = start;
                start++;
            }
        }
        it = it.next();

        // edge ki
        //which_e = eInd[it.edge()];
        if (eStart[it.edge()] == -1) {
            eStart[it.edge()] = index;
            start = index;
            index += subdiv_level - 1;
        } else {
            start = eStart[it.edge()];
        }
        if (it.edge().halfedge() == it) {
            // j is canonical vertex
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[2*subdiv_level + pos] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                toFin[3*subdiv_level - pos] = start;
                start++;
            }
        }

        for (int i = 3*subdiv_level; i < finalPoints.size()/2; i++) {
            toFin[i] = index;
            index++;
        }
        
        for (int i = 0; i < triangles.size()/3; i++) {
            int a = triangles[3*i], b = triangles[3*i + 1], c = triangles[3*i + 2];
            // update connectivity
            polygons.push_back({toFin[a], toFin[b], toFin[c]});
            // add correct distances
            correct_dist.push_back(make_tuple(toFin[a], toFin[b], 
            sqr(finalPoints[2*a] - finalPoints[2*b]) + sqr(finalPoints[2*a +1] - finalPoints[2*b + 1])));
            correct_dist.push_back(make_tuple(toFin[b], toFin[c], 
            sqr(finalPoints[2*b] - finalPoints[2*c]) + sqr(finalPoints[2*b +1] - finalPoints[2*c + 1])));
            correct_dist.push_back(make_tuple(toFin[c], toFin[a], 
            sqr(finalPoints[2*c] - finalPoints[2*a]) + sqr(finalPoints[2*c +1] - finalPoints[2*a + 1])));
        }
        indexing[f] = toFin;
    }
    CATmesh = std::unique_ptr<ManifoldSurfaceMesh>(new ManifoldSurfaceMesh(polygons));
    std::map<size_t, std::map<size_t, Edge>> intrinsicEdgeMap;
    VertexData<size_t> vMap = CATmesh->getVertexIndices();
    for (Edge e : CATmesh->edges()) {
        size_t v1 = vMap[e.halfedge().vertex()];
        size_t v2 = vMap[e.halfedge().twin().vertex()];
        intrinsicEdgeMap[v1][v2] = e;
        intrinsicEdgeMap[v2][v1] = e;
    }
    subdiv_points = vector<Vector3>(index, Vector3{0, 0, 0});
    for (Face f : mesh->faces()) {
        auto toFin = indexing[f];
        auto triangle = CATs[f];
        auto finalPoints = f_points[f];
        // initialize initial guess
        for (int i = 0; i < finalPoints.size()/2; i++) {
            vec3 stuff = triangle.planeToBary(vec2(finalPoints[2*i], finalPoints[2*i + 1]));
            subdiv_points[toFin[i]] = bary(f,stuff.x, stuff.y, stuff.z);
        }
    }
    EdgeData<double> edgeLengths(*CATmesh);
    size_t i1, i2;
    for (auto t : correct_dist) {
        i1 = std::get<0>(t);
        i2 = std::get<1>(t);
        double distsq = std::get<2>(t);
        edgeLengths[intrinsicEdgeMap[i1][i2]] = sqrt(distsq);
    }
    intrinsicGeometry = std::unique_ptr<EdgeLengthGeometry>(new EdgeLengthGeometry(*CATmesh, edgeLengths));
    // Building bending energy
    
    intrinsicGeometry->requireCotanLaplacian();
    intrinsicGeometry->requireVertexLumpedMassMatrix();
    Eigen::SparseMatrix<double> L = intrinsicGeometry->cotanLaplacian;

    //cout << "Laplacian: " << L.norm();
    Eigen::SparseMatrix<double> M = intrinsicGeometry->vertexLumpedMassMatrix.cwiseInverse();
    //cout << "Mass: " << M.norm();
    bendingMatrix = L.transpose() * M * L;
    //cout << "bending matrix made" << endl;

}



// registers subdiv_points to the polyscope view
void CatOpt::buildNewMesh() {
    CATpsMesh = polyscope::registerSurfaceMesh(
        "CAT Mesh",
        subdiv_points, CATmesh->getFaceVertexList(),
        polyscopePermutations(*CATmesh));
}
// returns \sum (dist^2 - dist_actual^2) + bending energy
double CatOpt::objective(const VectorXd &x1, const VectorXd &x2, const VectorXd &x3) {
    double result = 0.;
    int i1, i2;
    double distsq, actualDistsq;
    for (auto t : correct_dist) {
        i1 = std::get<0>(t);
        i2 = std::get<1>(t);
        distsq = std::get<2>(t);
        actualDistsq = (sqr(x1[i1] - x1[i2])+
                    sqr(x2[i1] - x2[i2])  +
                    sqr(x3[i1] - x3[i2]));
        result += sqr(log(actualDistsq/distsq))/2;
    }
    // bending energy
    result += 0.5 * bendingWeight * (x1.transpose() * bendingMatrix * x1)(0, 0);
    result += 0.5 * bendingWeight * (x2.transpose() * bendingMatrix * x2)(0, 0);
    result += 0.5 * bendingWeight * (x3.transpose() * bendingMatrix * x3)(0, 0);
    return result;
}
// gradient of metric embedding + bending energy
tuple<VectorXd, VectorXd, VectorXd> CatOpt::gradient(const VectorXd &x1, const VectorXd &x2, const VectorXd &x3, VectorXd &grad1, VectorXd &grad2, VectorXd &grad3) {
    int i1, i2;
    double distsq, actualDistsq, diff;
    grad1.setZero(subdiv_points.size());
    grad2.setZero(subdiv_points.size());
    grad3.setZero(subdiv_points.size()); 
    for (auto t : correct_dist) {
        i1 = std::get<0>(t);
        i2 = std::get<1>(t);
        distsq = std::get<2>(t);
        actualDistsq = (sqr(x1[i1] - x1[i2])+
                    sqr(x2[i1] - x2[i2])  +
                    sqr(x3[i1] - x3[i2]));
        diff = 2 * log(sqrt(actualDistsq)/sqrt(distsq))/actualDistsq;
        // update gradient
        grad1[i1] += diff * (x1[i1] - x1[i2]);
        grad1[i2] += diff * (x1[i2] - x1[i1]);
        grad2[i1] += diff * (x2[i1] - x2[i2]);
        grad2[i2] += diff * (x2[i2] - x2[i1]);
        grad3[i1] += diff * (x3[i1] - x3[i2]);
        grad3[i2] += diff * (x3[i2] - x3[i1]);
    }
    // bending energy
    grad1 += bendingWeight * bendingMatrix * x1;
    grad2 += bendingWeight * bendingMatrix * x2;
    grad3 += bendingWeight * bendingMatrix * x3;
    return make_tuple(grad1, grad2, grad3);
}

// runs the actual optimization, updating subdiv_points and the polyscope view
void CatOpt::step(int n) {
    //for (auto& v : subdiv_points) v *= 1.1;
    cout << "Starting descent" << endl;
    VectorXd grad1 = VectorXd::Zero(subdiv_points.size());
    VectorXd grad2 = VectorXd::Zero(subdiv_points.size());
    VectorXd grad3 = VectorXd::Zero(subdiv_points.size());
    double grad_size = 1.;
    VectorXd x1(subdiv_points.size());
    VectorXd x2(subdiv_points.size());
    VectorXd x3(subdiv_points.size());
    for (int i = 0; i < subdiv_points.size(); i++) {
        x1[i] = subdiv_points[i].x;
        x2[i] = subdiv_points[i].y;
        x3[i] = subdiv_points[i].z;
    }
    VectorXd x1_new = x1;
    VectorXd x2_new = x2;
    VectorXd x3_new = x3;
    for (int m = 0; m < n; m++){
    //while (sqrt(grad_size) > ep) {
        double result = objective(x1, x2, x3);
        gradient(x1, x2, x3, grad1, grad2, grad3);
        grad_size = grad_norm_sq(grad1, grad2, grad3);
        double t = 1.;
        x1_new = x1 - t * grad1;
        x2_new = x2 - t * grad2;
        x3_new = x3 - t * grad3;
        while (objective(x1_new, x2_new, x3_new) > result - alpha * t * grad_size) {
            t = beta * t;
            x1_new = x1 - t * grad1;
            x2_new = x2 - t * grad2;
            x3_new = x3 - t * grad3;
        }
        x1 = x1_new;
        x2 = x2_new;
        x3 = x3_new;
        if (iter % 100 == 0) {
            
            cout << "Starting iteration " << iter << endl;
            cout << "grad size squared:" << grad_size << endl;
            cout << "objective:" << result << endl;
            /*
           cout << "updating";
            for (int i = 0; i < subdiv_points.size(); i++) {
                subdiv_points[i].x = x1[i];
                subdiv_points[i].y = x2[i];
                subdiv_points[i].z = x3[i];
            }
            //return pos;
            CATpsMesh->updateVertexPositions(subdiv_points);
            */
            }
        iter++;
    }
    //cout << "iteration count: " << iter << endl;
    //cout << "grad size:" << sqrt(grad_size) << endl;
    //cout << "objective:" << objective(x1, x2, x3) << endl;
    //vector<Vector3> pos = subdiv_points;
    VertexData<Vector3> positions(*CATmesh);
    for (int i = 0; i < subdiv_points.size(); i++) {
        subdiv_points[i].x = x1[i];
        subdiv_points[i].y = x2[i];
        subdiv_points[i].z = x3[i];
        positions[i] = subdiv_points[i];
    }
    //return pos;
    CATpsMesh->updateVertexPositions(subdiv_points);
    //auto embedded = VertexPositionGeometry(*CATmesh, positions);
    auto embedded = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*CATmesh, positions));
    writeSurfaceMesh(*CATmesh, *embedded, "embedded_opt.obj"); 
}

// imgui wrapper
void CatOpt::myCallback() {
  if (ImGui::Button("do work")) {
    step(100);
  }
}