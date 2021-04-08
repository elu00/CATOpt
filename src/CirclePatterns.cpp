#include "CirclePatterns.h"

CirclePatterns::CirclePatterns(std::shared_ptr<ManifoldSurfaceMesh> mesh0, int optScheme0, vector<double>& solve,
EdgeData<size_t> eInd, VertexData<size_t> vInd, CornerData<size_t> cInd, FaceData<size_t> fInd):
angles(mesh->nHalfedges()),
thetas(mesh->nEdges()),
radii(mesh->nFaces()-1),
eIntIndices(mesh->nEdges()),
imaginaryHe(0),
eInd(eInd),
vInd(vInd),
cInd(cInd),
fInd(fInd),
OptScheme(optScheme0)
{
    mesh = mesh0;
    solver.n = mesh->nFaces()-1;
    sol = sol;
}


void CirclePatterns::setThetas()
{
    // TODODOJDOJODODJODJOJDOJDOJOD
    // set opposite angles
    /*
    int a = 0;
    for (Halfedge he : mesh->halfedges()) {
        if (!he.onBoundary) angles[he.getIndex()] = he->angle() + mosekSolver.xx[a++];
        else angles[he.getIndex()] = 0.0;
    }
    
    // set thetas
    for (Edge e = mesh->edges()) {
        HalfEdgeCIter he = e->he;
        thetas[eInd[e]] = M_PI - angles[he.getIndex()] - angles[he.twin().getIndex()];
    }
    */
}

bool CirclePatterns::computeAngles()
{
    int variables = 3 * (int)(mesh->nFaces() - mesh->nExteriorHalfedges()); 
    int constraints = (int)(mesh->nVertices() + mesh->nEdges() + mesh->nFaces() -
                            imaginaryHe - mesh->nExteriorHalfedges());
    int numanz = 3 * variables - imaginaryHe;
    int numqnz = variables;
    setThetas();

    return true;
}
inline double Cl2(double x)
{
    if (x == 0.0) return 0.0;
    x = std::remainder(x, 2*M_PI);
    if (x == 0.0) return 0.0;
    
    if (fabs(x) <= 2.0944) {
        double xx = x * x;
        return ((((((((((((2.3257441143020875e-22 * xx
                           + 1.0887357368300848e-20) * xx
                           + 5.178258806090624e-19) * xx
                           + 2.5105444608999545e-17) * xx
                           + 1.2462059912950672e-15) * xx
                           + 6.372636443183181e-14) * xx
                           + 3.387301370953521e-12) * xx
                           + 1.8978869988971e-10) * xx
                           + 1.1482216343327455e-8) * xx
                           + 7.873519778281683e-7) * xx
                           + 0.00006944444444444444) * xx
                           + 0.013888888888888888) * xx
                           - log(fabs(x)) + 1.0) * x;
    }
    
    x += ((x > 0.0) ? - M_PI : M_PI);
    double xx = x * x;
    return ((((((((((((3.901950904063069e-15 * xx
                       + 4.566487567193635e-14) * xx
                       + 5.429792727596476e-13) * xx
                       + 6.5812165661369675e-12) * xx
                       + 8.167010963952222e-11) * xx
                       + 1.0440290284867003e-9) * xx
                       + 1.3870999114054669e-8) * xx
                       + 1.941538399871733e-7) * xx
                       + 2.927965167548501e-6) * xx
                       + 0.0000496031746031746) * xx
                       + 0.0010416666666666667) * xx
                       + 0.041666666666666664) * xx
                       + log(0.5)) * x;
}

double ImLi2Sum(double dp, double theta)
{
    double tStar = M_PI - theta;
    double x = 2*atan(tanh(0.5*dp) * tan(0.5*tStar));
    
    return x*dp + Cl2(x + tStar) + Cl2(-x + tStar) - Cl2(2.0*tStar);
}

double fe(double dp, double theta)
{
    return atan2(sin(theta), exp(dp) - cos(theta));
}

void CirclePatterns::computeEnergy(double& energy, const Eigen::VectorXd& rho)
{
    energy = 0.0;
    
    // sum over edges
    for (Edge e : mesh->edges()) {
        int fk = e.halfedge().face().getIndex();
        
        if (e.isBoundary()) {
            energy -= 2*(M_PI - thetas[e.getIndex()]) * rho[fk];
            
        } else {
            int fl = e.halfedge().twin().face().getIndex();
            energy += ImLi2Sum(rho[fk] - rho[fl], thetas[e.getIndex()]) -
                      (M_PI - thetas[e.getIndex()])*(rho[fk] + rho[fl]);
        }
    }
    
    // sum over faces
    for (Face f: mesh->faces()) {
        if (!f.isBoundaryLoop()) energy += 2*M_PI*rho[fInd[f]];
    }
}

void CirclePatterns::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho)
{
    // loop over faces
    for (Face f : mesh->faces()) {
        if (!f.isBoundaryLoop()) {
            int fk = fInd[f];
            gradient[fk] = 2*M_PI;
            
            // sum of adjacent edges
            Halfedge he = f.halfedge();
            do {
                Edge e = he.edge();
                if (e.isBoundary()) {
                    gradient[fk] -= 2*(M_PI - thetas[eInd[e]]);
                    
                } else {
                    Halfedge h = e.halfedge();
                    int fl = fk == (int)fInd[h.face()] ? fInd[h.twin().face()] : fInd[h.face()];
                    gradient[fk] -= 2*fe(rho[fk] - rho[fl], thetas[e.getIndex()]);
                }
                
                he = he.next();
            } while (he != f.halfedge());
        }
    }
}

void CirclePatterns::computeHessian(Eigen::SparseMatrix<double>& hessian, const Eigen::VectorXd& rho)
{
    std::vector<Eigen::Triplet<double>> HTriplets;
    
    for (Edge e : mesh->edges()) {
        if (!e.isBoundary()) {
            int fk = fInd[e.halfedge().face()];
            int fl = fInd[e.halfedge().twin().face()];
                        
            double hessval = sin(thetas[eInd[e]]) / (cosh(rho(fk) - rho(fl)) - cos(thetas[eInd[e]]));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fk, hessval + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fl, hessval + 1e-8));
            HTriplets.push_back(Eigen::Triplet<double>(fk, fl, -hessval));
            HTriplets.push_back(Eigen::Triplet<double>(fl, fk, -hessval));
        }
    }
    
    hessian.setFromTriplets(HTriplets.begin(), HTriplets.end());
}

void CirclePatterns::setRadii()
{
    for (Face f : mesh->faces()) {
        if (!f.isBoundaryLoop()) radii[fInd[f]] = exp(solver.x[fInd[f]]);
    }
}

bool CirclePatterns::computeRadii()
{
    MeshHandle handle;
    handle.computeEnergy = std::bind(&CirclePatterns::computeEnergy, this, _1, _2);
    handle.computeGradient = std::bind(&CirclePatterns::computeGradient, this, _1, _2);
    handle.computeHessian = std::bind(&CirclePatterns::computeHessian, this, _1, _2);
    
    solver.handle = &handle;
    if (OptScheme == GRAD_DESCENT) solver.gradientDescent();
    else if (OptScheme == NEWTON) solver.newton();
    else if (OptScheme == TRUST_REGION) solver.trustRegion();
    else solver.lbfgs();
    
    // set radii
    setRadii();
    
    return true;
}

void CirclePatterns::computeAnglesAndEdgeLengths(Eigen::VectorXd& lengths)
{
    for (Edge e : mesh->edges()) {
        Halfedge h1 = e.halfedge();
        
        if (e.isBoundary()) {
            angles[h1.getIndex()] = M_PI - thetas[eInd[e]];
        
        } else {
            Halfedge h2 = h1.twin();
            double dp = log(radii[h1.face().getIndex()]) - log(radii[h2.face().getIndex()]);
            angles[h1.getIndex()] = fe(dp, thetas[e.getIndex()]);
            angles[h2.getIndex()] = fe(-dp, thetas[e.getIndex()]);
        }
        
        lengths[eInd[e]] = 2.0*radii[h1.face().getIndex()]*sin(angles[h1.getIndex()]);
    }
}

void CirclePatterns::performFaceLayout(Halfedge he, const Eigen::Vector2d& dir,
                                       Eigen::VectorXd& lengths, std::unordered_map<int, bool>& visited,
                                       std::stack<Edge>& stack)
{
    if (he.isInterior()) {
        int fIdx = he.face().getIndex();
        if (visited.find(fIdx) == visited.end()) {
            Halfedge next = he.next();
            Halfedge prev = he.next().next();
            
            // compute new uv position
            double angle = angles[next.getIndex()];
            Eigen::Vector2d newDir = {cos(angle)*dir[0] - sin(angle)*dir[1],
                                      sin(angle)*dir[0] + cos(angle)*dir[1]};
            
            // TODODOODOD
            //prev->vertex->uv = he->vertex->uv + newDir*lengths[eInd[prev.edge()]];
            
            // mark face as visited
            visited[fIdx] = true;
            
            // push edges onto stack
            stack.push(next.edge());
            stack.push(prev.edge());
        }
    }
}

void CirclePatterns::setUVs()
{
    // compute edge lengths
    Eigen::VectorXd lengths(mesh->nEdges());
    computeAnglesAndEdgeLengths(lengths);
    
    // push any edge
    std::stack<Edge> stack;
    Edge e = *(mesh->edges().begin());
    stack.push(e);

    //TODODODODOD
    //e.halfedge().vertex()->uv = Eigen::Vector2d::Zero();
    //e.halfedge().next().vertex()->uv = Eigen::Vector2d(lengths[eInd[e]], 0);
    
    // perform layout
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        Edge e = stack.top();
        stack.pop();
        
        Halfedge h1 = e.halfedge();
        Halfedge h2 = h1.twin();
        
        // compute edge vector
        //TODODODOD
        //Eigen::Vector2d dir = h2->vertex->uv - h1->vertex->uv;
        // delete this
        Eigen::Vector2d dir = Eigen::Vector2d::Zero();
        //dir.normalize();
        
        performFaceLayout(h1, dir, lengths, visited, stack);
        performFaceLayout(h2, -dir, lengths, visited, stack);
    }
    
    //normalize();
}

void CirclePatterns::parameterize()
{
    // set interior edge indices
    int eIdx = 0;
    for (Edge e : mesh->edges()) {
        if (!e.isBoundary()) eIntIndices[eInd[e]] = eIdx++;
        else {
            eIntIndices[eInd[e]] = -1;
            imaginaryHe++;
        }
    }
    
    // compute angles
    if (!computeAngles()) {
        std::cout << "Unable to compute angles" << std::endl;
        return;
    }
    
    // compute radii
    if (!computeRadii()) {
        std::cout << "Unable to compute radii" << std::endl;
        return;
    }
    
    // set uvs
    setUVs();
}
