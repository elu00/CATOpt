#include "CirclePatterns.h"
#include <Eigen/SparseQR>

CirclePatterns::CirclePatterns(shared_ptr<ManifoldSurfaceMesh> mesh0, Vertex infVertex, 
    EdgeData<bool> eMask, EdgeData<bool> eBdry,FaceData<bool> fMask, int optScheme0,
     Eigen::VectorXd thetas):
mesh(mesh0),
infVertex(infVertex),
eMask(eMask),
eBdry(eBdry),
fMask(fMask),
angles(mesh->nHalfedges()),
thetas(thetas),
radii(mesh->nFaces()),
eIntIndices(mesh->nEdges()),
imaginaryHe(0),
OptScheme(optScheme0)
{
    // I added a plus 1 here and at radii; should figure why I need to
    solver.n = mesh->nFaces();
    uv = VertexData<Eigen::Vector2d> (*mesh);
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    fInd = mesh->getFaceIndices();
    hInd = mesh->getHalfedgeIndices();
}

inline double Cl2(double x) {
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

double ImLi2Sum(double dp, double theta) {
    double tStar = M_PI - theta;
    double x = 2*atan(tanh(0.5*dp) * tan(0.5*tStar));
    
    return x*dp + Cl2(x + tStar) + Cl2(-x + tStar) - Cl2(2.0*tStar);
}

double fe(double dp, double theta) {
    return atan2(sin(theta), exp(dp) - cos(theta));
}

void CirclePatterns::computeEnergy(double& energy, const Eigen::VectorXd& rho)
{
    energy = 0.0;

    // sum over edges
    for (Edge e : mesh->edges()) {
        if (eMask[e]) {
            int fk = fInd[e.halfedge().face()];

            if (eBdry[e]) {
                energy -= 2 * (M_PI - thetas[eInd[e]]) * rho[fk];

            } else {
                int fl = fInd[e.halfedge().twin().face()];
                energy += ImLi2Sum(rho[fk] - rho[fl], thetas[eInd[e]]) -
                          (M_PI - thetas[eInd[e]]) * (rho[fk] + rho[fl]);
            }
        }
    }

    // sum over faces
    for (Face f: mesh->faces()) {
        if (!f.isBoundaryLoop() && fMask[f]) energy += 2*M_PI*rho[fInd[f]];
    }
}

void CirclePatterns::computeGradient(Eigen::VectorXd& gradient, const Eigen::VectorXd& rho)
{
    // loop over faces
    for (Face f : mesh->faces()) {
        if (!f.isBoundaryLoop() && fMask[f]) {
            int fk = fInd[f];
            gradient[fk] = 2*M_PI;
            
            // sum of adjacent edges
            Halfedge he = f.halfedge();
            do {
                Edge e = he.edge();
                if (eBdry[e]) {
                    gradient[fk] -= 2*(M_PI - thetas[eInd[e]]);
                    
                } else {
                    Halfedge h = e.halfedge();
                    int fl = fk == (int)fInd[h.face()] ? fInd[h.twin().face()] : fInd[h.face()];
                    gradient[fk] -= 2*fe(rho[fk] - rho[fl], thetas[eInd[e]]);
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
        if (!eBdry[e] && eMask[e]) {
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
        if (!f.isBoundaryLoop() && fMask[f]) radii[fInd[f]] = exp(solver.x[fInd[f]]);
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
        if(eMask[e]) {
            Halfedge h1 = e.halfedge();
        
            if (eBdry[e]) {
                angles[hInd[h1]] = M_PI - thetas[eInd[e]];
            
            } else {
                Halfedge h2 = h1.twin();
                double dp = log(radii[fInd[h1.face()]]) - log(radii[fInd[h2.face()]]);
                angles[hInd[h1]] = fe(dp, thetas[eInd[e]]);
                angles[hInd[h2]] = fe(-dp, thetas[eInd[e]]);
            }
        
            lengths[eInd[e]] = 2.0*radii[fInd[h1.face()]]*sin(angles[hInd[h1]]);
        }
    }
}

void CirclePatterns::performFaceLayout(Halfedge he, const Eigen::Vector2d& dir,
                                       Eigen::VectorXd& lengths, std::unordered_map<int, bool>& visited,
                                       std::stack<Edge>& stack)
{
    if (he.isInterior() && fMask[he.face()]) {
        int fIdx = fInd[he.face()];
        if (visited.find(fIdx) == visited.end()) {
            Halfedge next = he.next();
            Halfedge prev = he.next().next();
            
            // compute new uv position
            double angle = angles[hInd[next]];
            Eigen::Vector2d newDir = {cos(angle)*dir[0] - sin(angle)*dir[1],
                                      sin(angle)*dir[0] + cos(angle)*dir[1]};
            
            uv[prev.vertex()] = uv[he.vertex()] + newDir*lengths[eInd[prev.edge()]];
            
            // mark face as visited
            visited[fIdx] = true;
            
            // push edges onto stack
            if(eMask[next.edge()]) stack.push(next.edge());
            if(eMask[prev.edge()]) stack.push(prev.edge());
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
    Edge e0;
    for (Edge e :mesh->edges()) {
        if(eMask[e]) {
            e0 = e;
            break;
        }
    };
    stack.push(e0);

    uv[e0.halfedge().vertex()] = Eigen::Vector2d::Zero();
    uv[e0.halfedge().next().vertex()] = Eigen::Vector2d(lengths[eInd[e0]], 0);
    
    // perform layout
    std::unordered_map<int, bool> visited;
    while (!stack.empty()) {
        Edge e = stack.top();
        stack.pop();
        
        Halfedge h1 = e.halfedge();
        Halfedge h2 = h1.twin();
        
        // compute edge vector

        Eigen::Vector2d dir = uv[h2.vertex()] - uv[h1.vertex()];

        dir.normalize();
        // boundary edges
        performFaceLayout(h1, dir, lengths, visited, stack);
        performFaceLayout(h2, -dir, lengths, visited, stack);
    }
    
    normalize();
}

VertexData<Eigen::Vector2d> CirclePatterns::parameterize() {
    // set interior edge indices
    int eIdx = 0;
    for (Edge e : mesh->edges()) {
        if (eMask[e]) {
            if (!eBdry[e])
                eIntIndices[eInd[e]] = eIdx++;
            else {
                eIntIndices[eInd[e]] = -1;
                imaginaryHe++;
            }
        }
    }

    // compute radii
    if (!computeRadii()) {
        std::cout << "Unable to compute radii" << std::endl;
        return VertexData<Eigen::Vector2d>();
    }
    
    // set uvs
    setUVs();
    return uv;
}
double CirclePatterns::uvArea(Face f) {
    if (f.isBoundaryLoop() || !fMask[f]) {
        return 0;
    }
    
    const Eigen::Vector2d& a(uv[f.halfedge().vertex()]);
    const Eigen::Vector2d& b(uv[f.halfedge().next().vertex()]);
    const Eigen::Vector2d& c(uv[f.halfedge().next().next().vertex()]);
    
    const Eigen::Vector2d u = b - a;
    const Eigen::Vector2d v = c - a;
    
    return 0.5 * (u.x()*v.y() - v.x()*u.y());
}

Eigen::Vector2d CirclePatterns::uvBarycenter(Face f) {
    if (f.isBoundaryLoop() || !fMask[f]) {
        return Eigen::Vector2d::Zero();
    }
    
    const Eigen::Vector2d& a(uv[f.halfedge().vertex()]);
    const Eigen::Vector2d& b(uv[f.halfedge().next().vertex()]);
    const Eigen::Vector2d& c(uv[f.halfedge().next().next().vertex()]);
    
    return (a + b + c) / 3.0;
}
void CirclePatterns::normalize() {
    // compute center
    double totalArea = 0;
    Eigen::Vector2d center = Eigen::Vector2d::Zero();
    //uv[vInd[infVertex]] = Eigen::Vector2d::Zero();
    // DEBUG:
    for (Vertex v: mesh->vertices()){
        if(!std::isfinite(uv[v].x())) {
            uv[v] = Eigen::Vector2d::Zero();
            cout << "nananana" << endl;
        }
        //cout << uv[v].x() << " " << uv[v].y() << endl;
    }
    /*
    uv[infVertex.getIndex()].x() = -8;
    uv[infVertex.getIndex()].y() = 5;
    */
    for (Face f : mesh->faces()) {
        if (fMask[f]){
            double area = uvArea(f);
            center += area * uvBarycenter(f);
            totalArea += area;
        }
    }
    center /= totalArea;
    
    // shift
    double r = 0.0;
    for (Vertex v : mesh->vertices()) {
        //if (v != infVertex) {
            
            uv[v] -= center;
            r = std::max(r, uv[v].squaredNorm());
        //}
    }
    
    // scale
    r = sqrt(r);
    for (Vertex v : mesh->vertices()) {
        uv[v] /= r;
    }
}

inline double shift(double c) {
    return (c + 2.) * 500;
}