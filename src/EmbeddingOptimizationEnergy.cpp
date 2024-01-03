#include "EmbeddingOptimization.h"

#include "polyscope/surface_mesh.h"
// ======================= Optimization code==========================
void EmbeddingOptimization::initializeLM() {
    if (!intrinsicQuantitiesInitialized) {
        throw std::runtime_error(
            "Error: intrinsic checkboard/mesh not initialized yet");
        return;
    }
    currentSolution = initialSolution;

    int nCoordinates = 3 * nSubdividedVertices;

    // number of quads in the fine mesh
    size_t nQuads = nCorners * (n - 1) * (n - 1);
    size_t nFairVertices = fairVertices.size();

    LMValues = 3 * nQuads + nFairVertices + 1;

    targetVolume = 10 * calculateVolume(currentSolution);
    cout << "target volume is " << targetVolume << endl;

    // make sure that the size of the vertex position vector x agrees with the
    // number of fine vertices in the mesh times three
    assert(currentSolution.size() == nCoordinates);
    LMInputs = nCoordinates;
    LMInitialized = true;
    return;
}
// One step of Levenberg-Marquardt.
// For reference, see
// https://www.dmg.tuwien.ac.at/geom/ig/publications/isoforfab/isoforfab.pdf#cite.madsen04
// and http://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf for pseudocode
void EmbeddingOptimization::LMOneStep(int MAX_ITERS) {
    if (!LMInitialized) {
        throw std::runtime_error("Error: LM stuff not initialized yet");
        return;
    }

    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorType;
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
        // cout << "DescentMatrix norm: " << DescentMatrix.norm() << endl;
        Eigen::VectorXd hLM = -solvePositiveDefinite(DescentMatrix, g);
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
            double denominator = (hLM.dot(mu * hLM - g)) / 2.;
            double rho = numerator / denominator;
            if (rho > 0) {
                currentSolution = xStar;
                J.setZero();
                evaluateJacobian(currentSolution, J);
                // compute the relevant things
                A = J.transpose() * J;
                f = fNew;
                g = J.transpose() * f;
                found = (g.norm() < eps1);
                mu = mu * std::max(1. / 3., 1 - cube(2 * rho - 1));
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

// Call the optimization procedure for one step, then update the mesh in
// polyscope with the new vertex positions
void EmbeddingOptimization::optimizeOneStep(int MAX_ITERS) {
    LMOneStep(MAX_ITERS);
    VertexData<Vector3> positions(*submesh);
    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Vertex v : submesh->vertices()) {
        size_t i = subVertexIndices[v];
        positions[v] = {currentSolution[3 * i], currentSolution[3 * i + 1],
                        currentSolution[3 * i + 2]};
        // cout << positions[v] << endl;
    }
    // visualize fairness term

    vector<Vector3> nodes;
    vector<std::array<size_t, 2>> edges;
    for (auto indices : fairVertices) {
        size_t centerIndex = indices[0];

        size_t centerNodeIndex = nodes.size();
        Vector3 center = {currentSolution[3 * centerIndex],
                          currentSolution[3 * centerIndex + 1],
                          currentSolution[3 * centerIndex + 2]};
        nodes.push_back(center);

        for (int n = 1; n < indices.size(); n++) {
            size_t index = indices[n];

            size_t neighborNodeIndex = nodes.size();
            Vector3 neighbor = {currentSolution[3 * index],
                                currentSolution[3 * index + 1],
                                currentSolution[3 * index + 2]};

            nodes.push_back(2 * center / 3 + neighbor / 3);
            edges.push_back({centerNodeIndex, neighborNodeIndex});
        }
    }
    polyscope::registerCurveNetwork("fairness terms", nodes, edges);
    // calculate relative error
    vector<double> error;
    for (size_t quadIndex = 0; quadIndex < quads.size(); quadIndex++) {
        // vertex indices
        auto [i, j, k, l] = quads[quadIndex];
        error.push_back(abs(
            (positions[i] - positions[k]).norm2() / c_iso_0[quadIndex] - 1));
        error.push_back(abs(
            (positions[j] - positions[l]).norm2() / c_iso_1[quadIndex] - 1));
    }
    cout << "average relative error: "
         << std::reduce(error.begin(), error.end()) / error.size() << endl;
    cout << "max relative error: "
         << (*std::max_element(error.begin(), error.end())) << endl;
    cout << "min relative error: "
         << (*std::min_element(error.begin(), error.end())) << endl;

    polyscope::registerSurfaceMesh("Optimized Mesh", positions,
                                   submesh->getFaceVertexList());
}

// ======================= Energy evaluation code==========================
void EmbeddingOptimization::evaluateEnergy(
    const Eigen::VectorXd &v, // vertex positions
    Eigen::VectorXd &energy)  // each entry of this vector is the (square root
                              // of) a term in the energy summand
{

    size_t nQuads = nCorners * (n - 1) * (n - 1);
    // make sure we were given the right number of vertex coordinates
    assert(v.size() == 3 * nSubdividedVertices);
    // make sure the given energy summand vector is the right size; there are 3
    // * nQuads
    assert(energy.size() == 3 * nQuads + fairVertices.size());
    // fill the energy summand vector with the individual terms
    for (size_t quadIndex = 0; quadIndex < nQuads; quadIndex++) {
        // vertex indices
        auto [i, j, k, l] = quads[quadIndex];
        addLengthTerm(energy, v, quadIndex, i, k, c_iso_0[quadIndex]);
        addLengthTerm(energy, v, quadIndex + nQuads, j, l, c_iso_1[quadIndex]);
        addAngleTerm(energy, v, quadIndex + 2 * nQuads, i, j, k, l,
                     c_iso_2[quadIndex]);
    }
    // center regularizer
    for (int i = 0; i < fairVertices.size(); i++) {
        addCenterTerm(energy, v, 3 * nQuads + i, fairVertices[i]);
    }
    // volume term
    addVolumeTerm(energy, v, 3 * nQuads + fairVertices.size(), targetVolume);
}
void EmbeddingOptimization::evaluateJacobian(const Eigen::VectorXd &v,
                                             Eigen::SparseMatrix<double> &J) {
    // number of fine quads
    size_t nQuads = nCorners * (n - 1) * (n - 1);
    vector<Eigen::Triplet<double>> tripletList;

    // isometry term entries
    for (size_t quadIndex = 0; quadIndex < nQuads; quadIndex++) {
        // vertex indices
        auto [i, j, k, l] = quads[quadIndex];
        addLengthGradient(tripletList, v, quadIndex, i, k, c_iso_0[quadIndex]);
        addLengthGradient(tripletList, v, quadIndex + nQuads, j, l,
                          c_iso_1[quadIndex]);
        addAngleGradient(tripletList, v, quadIndex + 2 * nQuads, i, j, k, l,
                         c_iso_2[quadIndex]);
    }
    // regularization term
    for (int i = 0; i < fairVertices.size(); i++) {
        addCenterGradient(tripletList, v, 3 * nQuads + i, fairVertices[i]);
    }
    // volume term
    addVolumeGradient(tripletList, v, 3 * nQuads + fairVertices.size());
    // build Jacobian matrix from triplets
    J.setFromTriplets(tripletList.begin(), tripletList.end());
}

// ======================= Energy components==========================
// Given indices for vertices i, j
// adds the term (|i-j|^2 - target) and it's corresponding gradient to the
// energy
inline void EmbeddingOptimization::addLengthTerm(Eigen::VectorXd &energy,
                                                 const Eigen::VectorXd &v,
                                                 size_t energyIndex,
                                                 size_t iIndex, size_t jIndex,
                                                 double target) {
    Vector3 vi = indexToVector(iIndex, v);
    Vector3 vj = indexToVector(jIndex, v);

    const double normalization = 1;
    energy[energyIndex] = ((vi - vj).norm2() - target) * normalization;
}

inline void EmbeddingOptimization::addLengthGradient(
    vector<Eigen::Triplet<double>> &tripletList, const Eigen::VectorXd &v,
    size_t energyIndex, size_t iIndex, size_t jIndex, double target) {
    typedef Eigen::Triplet<double> T;

    const double normalization = 1;
    for (int d = 0; d < 3; d++) {
        // i gradients
        tripletList.push_back(
            T(energyIndex, 3 * iIndex + d,
              2 * (v[3 * iIndex + d] - v[3 * jIndex + d]) * normalization));
        // j gradients
        tripletList.push_back(
            T(energyIndex, 3 * jIndex + d,
              2 * (v[3 * jIndex + d] - v[3 * iIndex + d]) * normalization));
    }
}

inline void EmbeddingOptimization::addAngleTerm(
    Eigen::VectorXd &energy, const Eigen::VectorXd &v, size_t energyIndex,
    size_t iIndex, size_t jIndex, size_t kIndex, size_t lIndex, double target) {
    Vector3 vi = indexToVector(iIndex, v);
    Vector3 vj = indexToVector(jIndex, v);
    Vector3 vk = indexToVector(kIndex, v);
    Vector3 vl = indexToVector(lIndex, v);

    const double normalization = 1;
    energy[energyIndex] = (dot(vi - vk, vj - vl) - target) * normalization;
}

inline void EmbeddingOptimization::addAngleGradient(
    vector<Eigen::Triplet<double>> &tripletList, const Eigen::VectorXd &v,
    size_t energyIndex, size_t iIndex, size_t jIndex, size_t kIndex,
    size_t lIndex, double target) {
    typedef Eigen::Triplet<double> T;

    const double normalization = 1;
    // i partial: j - l
    for (int d = 0; d < 3; d++) {
        tripletList.push_back(
            T(energyIndex, 3 * iIndex,
              v[3 * jIndex + d] - v[3 * lIndex + d] * normalization));

        // k partial: l - j
        tripletList.push_back(
            T(energyIndex, 3 * kIndex,
              v[3 * lIndex + d] - v[3 * jIndex + d] * normalization));

        // j partial: i - k
        tripletList.push_back(
            T(energyIndex, 3 * jIndex,
              v[3 * iIndex + d] - v[3 * kIndex + d] * normalization));

        // l partial: k - i
        tripletList.push_back(
            T(energyIndex, 3 * lIndex,
              v[3 * kIndex + d] - v[3 * iIndex + d] * normalization));
    }
}

// indices is a vector consisting of [center term, neighbors].
// The corresponding contribution to the energy is || # neighbors * center
// position - Σ neighbors ||^2
inline void EmbeddingOptimization::addCenterTerm(Eigen::VectorXd &energy,
                                                 const Eigen::VectorXd &v,
                                                 size_t energyIndex,
                                                 vector<size_t> &indices) {
    size_t nNeighbors = indices.size() - 1;
    size_t centerIndex = indices[0];

    // the position of the desired center
    Vector3 center = indexToVector(centerIndex, v);
    // the center of neighbors
    Vector3 average = {0, 0, 0};

    for (int n = 1; n < indices.size(); n++) {
        size_t index = indices[n];
        average += indexToVector(index, v);
    }
    average /= nNeighbors;
    Vector3 difference = center - average;
    energy[energyIndex] = ((center - average).norm2()) * fairnessNormalization;
}

inline void EmbeddingOptimization::addCenterGradient(
    vector<Eigen::Triplet<double>> &tripletList, const Eigen::VectorXd &v,
    size_t energyIndex, vector<size_t> &indices) {
    typedef Eigen::Triplet<double> T;

    size_t nNeighbors = indices.size() - 1;
    size_t centerIndex = indices[0];

    // the position of the desired center
    Vector3 center = indexToVector(centerIndex, v);
    // the center of neighbors
    Vector3 average = {0, 0, 0};

    for (int n = 1; n < indices.size(); n++) {
        size_t neighborIndex = indices[n];
        average += indexToVector(neighborIndex, v);
    }
    average /= nNeighbors;
    Vector3 difference = center - average;
    // center index gradients
    tripletList.push_back(T(energyIndex, 3 * centerIndex,
                            2 * (difference.x) * fairnessNormalization));
    tripletList.push_back(T(energyIndex, 3 * centerIndex + 1,
                            2 * (difference.y) * fairnessNormalization));
    tripletList.push_back(T(energyIndex, 3 * centerIndex + 2,
                            2 * (difference.z) * fairnessNormalization));

    // neighbor gradients
    for (int n = 1; n < indices.size(); n++) {
        size_t neighborIndex = indices[n];
        tripletList.push_back(
            T(energyIndex, 3 * neighborIndex,
              -2 * (difference.x) * fairnessNormalization / nNeighbors));
        tripletList.push_back(
            T(energyIndex, 3 * neighborIndex + 1,
              -2 * (difference.y) * fairnessNormalization / nNeighbors));
        tripletList.push_back(
            T(energyIndex, 3 * neighborIndex + 2,
              -2 * (difference.z) * fairnessNormalization / nNeighbors));
    }
}

inline double EmbeddingOptimization::calculateVolume(const Eigen::VectorXd &v) {
    double volume = 0.0;
    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Face f : submesh->faces()) {
        vector<size_t> vertexIndices = {
            subVertexIndices[f.halfedge().vertex()],
            subVertexIndices[f.halfedge().next().vertex()],
            subVertexIndices[f.halfedge().next().next().vertex()],
            subVertexIndices[f.halfedge().next().next().next().vertex()]};
        vector<Vector3> V;
        for (auto ind : vertexIndices) {
            V.push_back(indexToVector(ind, v));
        }
        for (int i = 0; i < 4; i++) {
            volume += tripProd(V[i], V[(i + 1) % 4], V[(i + 2) % 4]) / 12.;
        }
    }
    cout << "current volume is" << volume << endl;
    return volume;
}
inline void EmbeddingOptimization::addVolumeTerm(Eigen::VectorXd &energy,
                                                 const Eigen::VectorXd &v,
                                                 size_t energyIndex,
                                                 double targetVolume) {
    energy[energyIndex] = (calculateVolume(v) - targetVolume) / 10;
}

inline void EmbeddingOptimization::addVolumeGradient(
    vector<Eigen::Triplet<double>> &tripletList, const Eigen::VectorXd &v,
    size_t energyIndex) {
    typedef Eigen::Triplet<double> T;
    VertexData<size_t> subVertexIndices = submesh->getVertexIndices();
    for (Face f : submesh->faces()) {
        vector<size_t> vertexIndices = {
            subVertexIndices[f.halfedge().vertex()],
            subVertexIndices[f.halfedge().next().vertex()],
            subVertexIndices[f.halfedge().next().next().vertex()],
            subVertexIndices[f.halfedge().next().next().next().vertex()]};
        vector<Vector3> V;
        for (auto ind : vertexIndices) {
            V.push_back(indexToVector(ind, v));
        }
        for (int i = 0; i < 4; i++) {
            Vector3 grad = cross(V[(i + 1) % 4], V[(i + 2) % 4]) +
                           cross(V[(i + 2) % 4], V[(i + 3) % 4]) +
                           cross(V[(i + 1) % 4], V[(i + 3) % 4]);
            tripletList.push_back(
                T(energyIndex, 3 * vertexIndices[i] + 0, grad.x / 120.));
            tripletList.push_back(
                T(energyIndex, 3 * vertexIndices[i] + 1, grad.y / 120.));
            tripletList.push_back(
                T(energyIndex, 3 * vertexIndices[i] + 2, grad.z / 120.));
        }
    }
}

// =================== Utility==========================================

// grabs the 3 indices associated to index
inline Vector3 EmbeddingOptimization::indexToVector(size_t index,
                                                    const Eigen::VectorXd &v) {
    return {v[3 * index], v[3 * index + 1], v[3 * index + 2]};
}

inline double EmbeddingOptimization::sqr(double x) { return x * x; }
inline double EmbeddingOptimization::cube(double x) { return x * x * x; }
// volume calculations
inline double EmbeddingOptimization::tripProd(Vector3 i, Vector3 j, Vector3 k) {
    return dot(i, cross(j, k));
}
