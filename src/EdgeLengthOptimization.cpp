#define EPSILON 1e-11
#define MAX_ITERS 20000

#include "EdgeLengthOptimization.h"
EdgeLengthOptimization::EdgeLengthOptimization(
    shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> eta)
    : mesh(mesh), eta(eta) {
    // Initialize quantities
    nEdges = mesh->nEdges();
    e_ = mesh->getEdgeIndices();
    lambda0 = Eigen::VectorXd::Zero(nEdges);
}
EdgeLengthOptimization::EdgeLengthOptimization(
    shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l,
    EdgeData<double> eta)
    : mesh(mesh), eta(eta) {
    // Initialize quantities
    nEdges = mesh->nEdges();
    e_ = mesh->getEdgeIndices();
    lambda0 = Eigen::VectorXd::Zero(nEdges);
    for (Edge e : mesh->edges()) {
        lambda0[e_[e]] = log(l[e]);
    }
}

// clausen integral
inline double EdgeLengthOptimization::Cl2(double x) {
    if (x == 0.0)
        return 0.0;
    x = std::remainder(x, 2 * M_PI);
    if (x == 0.0)
        return 0.0;

    if (fabs(x) <= 2.0944) {
        double xx = x * x;
        return ((((((((((((2.3257441143020875e-22 * xx +
                           1.0887357368300848e-20) *
                              xx +
                          5.178258806090624e-19) *
                             xx +
                         2.5105444608999545e-17) *
                            xx +
                        1.2462059912950672e-15) *
                           xx +
                       6.372636443183181e-14) *
                          xx +
                      3.387301370953521e-12) *
                         xx +
                     1.8978869988971e-10) *
                        xx +
                    1.1482216343327455e-8) *
                       xx +
                   7.873519778281683e-7) *
                      xx +
                  0.00006944444444444444) *
                     xx +
                 0.013888888888888888) *
                    xx -
                log(fabs(x)) + 1.0) *
               x;
    }

    x += ((x > 0.0) ? -M_PI : M_PI);
    double xx = x * x;
    return ((((((((((((3.901950904063069e-15 * xx + 4.566487567193635e-14) *
                          xx +
                      5.429792727596476e-13) *
                         xx +
                     6.5812165661369675e-12) *
                        xx +
                    8.167010963952222e-11) *
                       xx +
                   1.0440290284867003e-9) *
                      xx +
                  1.3870999114054669e-8) *
                     xx +
                 1.941538399871733e-7) *
                    xx +
                2.927965167548501e-6) *
                   xx +
               0.0000496031746031746) *
                  xx +
              0.0010416666666666667) *
                 xx +
             0.041666666666666664) *
                xx +
            log(0.5)) *
           x;
}

// returns the angle at corner i^jk given the edge lengths l_ij, l_jk, l_ki
double EdgeLengthOptimization::CornerAngle(double l_ij, double l_jk,
                                           double l_ki) {
    // cout << l_ij << " " << l_jk << " " << l_ki << endl;
    if (l_ij > l_jk + l_ki || l_ki > l_ij + l_jk)
        return 0.0;
    if (l_jk > l_ki + l_ij)
        return M_PI;
    return acos((l_ij * l_ij - l_jk * l_jk + l_ki * l_ki) / (2 * l_ij * l_ki));
}
// returns the angle at corner i^jk given the logarithmic edge lengths λ_ij,
// λ_jk, λ_ki
double EdgeLengthOptimization::ExpCornerAngle(double lambda_ij,
                                              double lambda_jk,
                                              double lambda_ki) {
    return CornerAngle(exp(lambda_ij), exp(lambda_jk), exp(lambda_ki));
}
// given logarithmic edge lengths lambda, returns the corresponding interior
// angles at each corner
CornerData<double>
EdgeLengthOptimization::InitializeCornerAngles(Eigen::VectorXd &lambda) {
    CornerData<double> theta(*mesh);
    for (Corner c : mesh->corners()) {
        // index into the appropriate edge positions
        double lambda_ij = lambda[e_[c.halfedge().edge()]];
        double lambda_jk = lambda[e_[c.halfedge().next().edge()]];
        double lambda_ki = lambda[e_[c.halfedge().next().next().edge()]];
        theta[c] = ExpCornerAngle(lambda_ij, lambda_jk, lambda_ki);
    }
    return theta;
}
// returns the edge length energy associated to the logarithmic edge lengths λ,
// which is defined via E= Σ_{a^bc} θ_a^bc λ_bc+ Cl(2θ_a^bc)/2 -Σ_{ij} λ_ij
// (𝜋-η_ij) where Cl denotes the Clausen integral
double EdgeLengthOptimization::LengthEnergy(Eigen::VectorXd &lambda) {
    // initialize energy
    double energy = 0.0;
    // initialize all corner angles
    CornerData<double> theta = InitializeCornerAngles(lambda);
    for (Corner c : mesh->corners()) {
        double lambda_jk = lambda[e_[c.halfedge().next().edge()]];
        energy += theta[c] * lambda_jk + Cl2(2 * theta[c]) / 2;
    }
    for (Edge e : mesh->edges()) {
        energy -= lambda[e_[e]] * (M_PI - eta[e]);
    }
    return energy;
}

Eigen::VectorXd
EdgeLengthOptimization::LengthEnergyGradient(Eigen::VectorXd &lambda) {
    Eigen::VectorXd g = Eigen::VectorXd::Zero(nEdges);
    // initialize corner angles
    CornerData<double> theta = InitializeCornerAngles(lambda);
    for (Edge e : mesh->edges()) {
        size_t eIndex = e_[e];
        // the gradient at each edge is given by
        // the sum of the two opposite interior angles
        // plus the intersection angle minus pi
        //            k
        //           /1\
        //          /   \
        //         i ----j
        //          \   /
        //           \2/
        //            l
        g[eIndex] = eta[e] - M_PI;
        for (Halfedge h : {e.halfedge(), e.halfedge().twin()}) {
            if (h.isInterior()) {
                g[eIndex] += theta[h.next().next().corner()];
            }
        }
    }
    return g;
}
Eigen::SparseMatrix<double>
EdgeLengthOptimization::LengthEnergyHessian(Eigen::VectorXd &lambda) {
    // initialize sparse matrix and triplet list
    typedef Eigen::Triplet<double> T;
    Eigen::SparseMatrix<double> H(nEdges, nEdges);
    vector<Eigen::Triplet<double>> tripletList;
    // initialize corner angles corresponding to edge lengths
    CornerData<double> theta = InitializeCornerAngles(lambda);
    for (Face f : mesh->faces()) {
        Halfedge h = f.halfedge();
        size_t ij = e_[h.edge()];
        size_t jk = e_[h.next().edge()];
        size_t ki = e_[h.next().next().edge()];
        double cot_ijk = 1 / tan(theta[h.corner()]);
        double cot_jki = 1 / tan(theta[h.next().corner()]);
        double cot_kij = 1 / tan(theta[h.next().next().corner()]);
        cout << cot_ijk << " " << cot_jki << " " << cot_kij << endl;
        // diagonal entries
        tripletList.push_back(T(ij, ij, cot_ijk + cot_jki));
        tripletList.push_back(T(jk, jk, cot_jki + cot_kij));
        tripletList.push_back(T(ki, ki, cot_kij + cot_ijk));
        // off diagonal entries
        tripletList.push_back(T(ij, jk, -cot_jki));
        tripletList.push_back(T(jk, ki, -cot_kij));
        tripletList.push_back(T(ki, ij, -cot_ijk));
        // symmetric off diagonal entries
        tripletList.push_back(T(jk, ij, -cot_jki));
        tripletList.push_back(T(ki, jk, -cot_kij));
        tripletList.push_back(T(ij, ki, -cot_ijk));
    }
    H.setFromTriplets(tripletList.begin(), tripletList.end());
    return H;
}
// Input: A triangulation with edge lengths l : E -> R>0,
// and circumcircle intersection angles 𝜂ij that are compatible
// with some coherent angle system. A parameter 𝜀 > 0 determines
// the stopping tolerance, and parameters 𝛼 ∈ (0, 1/2),
// 𝛽 ∈ (0, 1) control line search
// Output: New edge lengths l : 𝐸 → R that exhibit the prescribed
// intersection angles.
Eigen::VectorXd EdgeLengthOptimization::MinimizeLengthEnergy() {
    cout << "Starting length energy minimization" << endl;
    int k = 0;
    Eigen::VectorXd x = lambda0;

    const double alpha = 0.5;
    const double beta = 0.3;
    while (true) {
        // compute update direction
        // debug stuff
        cout << "iteration " << k << endl;
        for (int i = 0; i < x.size(); i++) {
            cout << x[i] << endl;
        }
        Eigen::VectorXd g = -LengthEnergyGradient(x);
        Eigen::SparseMatrix<double> H = LengthEnergyHessian(x);
        cout << "g size " << g.norm() << endl;
        cout << "H size " << H.norm() << endl;

        Eigen::VectorXd mu = solvePositiveDefinite(H, g);
        // debug
        cout << "step size " << mu.norm() << endl;
        if (abs(mu.dot(g)) <= 2 * EPSILON)
            break;

        // compute step size
        double t = 1.0;
        double initialEnergy = LengthEnergy(x);
        Eigen::VectorXd xp = x + t * mu;
        double currentEnergy = LengthEnergy(xp);
        while (currentEnergy > initialEnergy + alpha * t * g.dot(mu)) {
            t = beta * t;
            xp = x + t * mu;
            currentEnergy = LengthEnergy(xp);
        }

        // update
        x += t * mu;
        k++;

        // check termination condition
        if (fabs(initialEnergy - LengthEnergy(x)) < EPSILON || k > MAX_ITERS)
            break;
    }
    cout << "Iterations: " << k << endl;
    for (int i = 0; i < nEdges; i++) {
        x[i] = exp(x[i]);
    }
    return x;
}
