#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

#include "Bff.h"
#include "MeshIO.h"
#include "HoleFiller.h"
#include "Generators.h"
#include "ConePlacement.h"
#include "Cutter.h"

#include "args/args.hxx"
#include "imgui.h"

#include "fusion.h"
#include <algorithm>
#include <tuple>
#include <map>

#include "util.cpp"

#include "glm/vec3.hpp"

using namespace geometrycentral;
using namespace geometrycentral::surface;
using std::cout;
using std::endl;
using std::max;
using std::unique_ptr;

using namespace mosek::fusion;
using namespace monty;
using Eigen::VectorXd;
using std::make_tuple;
using std::tuple;
using std::vector;

// == Geometry-central data
string inputMeshPath;
unique_ptr<ManifoldSurfaceMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

unique_ptr<EdgeLengthGeometry> intrinsicGeometry;
unique_ptr<ManifoldSurfaceMesh> CATmesh;
//unique_ptr<VertexPositionGeometry> CATgeometry;

// Mesh data
size_t nVertices;
size_t nFaces;
size_t nEdges;
size_t nCorners;
EdgeData<size_t> eInd;
VertexData<size_t> vInd;
CornerData<size_t> cInd;
FaceData<size_t> fInd;
CornerData<double> cornerAngles;
VertexData<double> angleDefects;

// Optimization Stuff
// Vector<double> x_init;
vector<double> sol;
vector<double> rhs;
vector<double> ineqRHS0;
vector<double> ineqRHS1;

size_t iter = 0;
// specifies number of 
size_t subdiv_level = 5;
vector<Vector3> subdiv_points;
//vector<Vector3> grad;
vector<tuple<size_t, size_t, double>> correct_dist;
Eigen::SparseMatrix<double> bendingMatrix;
// tuning parameters for gradient descent
double alpha = 0.1;
double beta = 0.5;
double ep = 1e-3;
// weight for the bending energy
double bendingWeight = 1e-8;

//stuff for bff
bff::Mesh bffMesh;
bff::Model model;
vector<Vector2> flattened;
vector<double> alphas;
CornerData<double> targetAngles;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;
polyscope::SurfaceMesh *CATpsMesh;

inline double sqr(double x) { return x*x; }

void initializeQuantities() {
    // Initialization
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    nFaces = mesh->nFaces();
    eInd = mesh->getEdgeIndices();
    vInd = mesh->getVertexIndices();
    cInd = mesh->getCornerIndices();
    fInd = mesh->getFaceIndices();
    geometry->requireEdgeLengths();
    geometry->requireVertexGaussianCurvatures();
    geometry->requireVertexAngleSums();
    angleDefects = geometry->vertexGaussianCurvatures;
    sol = vector<double>(nEdges);
}

// calls mosek to generate optimal alphas, and places them into sol
void generateConstraints() {
    // Model initialization
    // TODO: check for optimality
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", nEdges, Domain::inRange(-2 * PI, 2 * PI));

    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    // Equality constraint initialization
    rhs = vector<double>(nVertices);
    for (size_t i = 0; i < nVertices; i++) {
        if (mesh-> vertex(i).isBoundary()) {
            rhs[i] = PI - geometry->vertexAngleSums[mesh->vertex(i)];
        } else {
            rhs[i] =  angleDefects[mesh->vertex(i)];                
        }
    }
    for (Edge e : mesh->edges()) {
        double weight = e.isBoundary() ? 1. : 2. ;
        rows.emplace_back(vInd[e.halfedge().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(weight);
        rows.emplace_back(vInd[e.halfedge().twin().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(weight);
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto Meq = Matrix::sparse(nVertices, nEdges, r, c, v);
    auto eqRHS = new_array_ptr(rhs);
    M->constraint("eq constraints", Expr::mul(Meq, x), Domain::equalsTo(eqRHS));
    //cout << "eq generated" << endl;
    // inequality constraints
    ineqRHS0 = vector<double>(nCorners);
    ineqRHS1 = vector<double>(nCorners);
    rows.clear();
    cols.clear();
    values.clear();
    for (Corner c : mesh->corners()) {
        ineqRHS0[cInd[c]] = -geometry->cornerAngle(c);
        ineqRHS1[cInd[c]] = 2 * PI - geometry->cornerAngle(c);
        Halfedge h = c.halfedge();
        rows.emplace_back(cInd[c]);
        cols.emplace_back(eInd[h.edge()]);
        values.emplace_back(1.);
        rows.emplace_back(cInd[h.next().corner()]);
        cols.emplace_back(eInd[h.edge()]);
        values.emplace_back(1.);
    }
    r = new_array_ptr<int>(rows);
    c = new_array_ptr<int>(cols);
    v = new_array_ptr<double>(values);
    auto Mineq = Matrix::sparse(nCorners, nEdges, r, c, v);
    auto inRHS0 = new_array_ptr(ineqRHS0);
    auto inRHS1 = new_array_ptr(ineqRHS1);
    // sum is greater than 0
    M->constraint("ineq0 constraints", Expr::mul(Mineq, x), Domain::greaterThan(inRHS0));
    // sum is less than 2pi
    M->constraint("ineq1 constraints", Expr::mul(Mineq, x), Domain::lessThan(inRHS1));

    auto ones = std::make_shared<ndarray<double, 1>>(shape(nEdges), 1.);
    //cout << "ineq generated" << endl;
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::vstack(t, x), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << M->getProblemStatus() << endl;
    //cout << x->level() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = x->getSize();
    auto xVal = x->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    for (int i = 0; i < xsize; ++i) {
        sol[i] = (*xVal)[i];
    }
}

inline Vector3 bary(Face f, double a, double b, double c) {
    //cout << a + b + c << endl;
    auto it = f.halfedge();
    Vector3 i = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 j = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 k = geometry->inputVertexPositions[it.vertex()];
    return a * i + b * j + c * k;
}

// initializes CATMesh, handles all subdivision energy initialization
void subdivision() {
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


        // DEBUG
        //std::ofstream s ("tri" + std::to_string(index) + ".svg", std::ofstream::out);
        //s << CAT(ij*500, jk*500, ki*500, a_ij, a_jk, a_ki, true).triangulation_svg(subdiv_level);
        //s.close();
        //int which_f = fInd[f], which_e;


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
void buildNewMesh() {
    //VertexData<Vector3> positions(*CATmesh);
    //for (size_t i = 0; i < subdiv_points.size(); i++) {
    //    positions[i] = subdiv_points[i];
    //}
    //CATgeometry = std::unique_ptr<VertexPositionGeometry>(new VertexPositionGeometry(*CATmesh, subdiv_points));
    CATpsMesh = polyscope::registerSurfaceMesh(
        "CAT Mesh",
        subdiv_points, CATmesh->getFaceVertexList(),
        polyscopePermutations(*CATmesh));
}
/*
 * Optimization code
 */
// returns \sum (dist^2 - dist_actual^2) + bending energy
double objective(const VectorXd &x1, const VectorXd &x2, const VectorXd &x3) {
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
tuple<VectorXd, VectorXd, VectorXd> gradient(const VectorXd &x1, const VectorXd &x2, const VectorXd &x3, VectorXd &grad1, VectorXd &grad2, VectorXd &grad3) {
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

// convenience function to return square of norm of gradient
double grad_norm_sq(const VectorXd &grad1, const VectorXd &grad2, const VectorXd &grad3) {
    return grad1.squaredNorm() + grad2.squaredNorm() + grad3.squaredNorm();
}
/*
void generateVisualization() {
    // Visualization
    psMesh->addEdgeScalarQuantity("Final Solution", sol);
    psMesh->addVertexScalarQuantity("curvature",
                                    geometry->vertexGaussianCurvatures);
    cout << "No of subdivision points:" << subdiv_points.size();
}
*/
// runs the actual optimization, updating subdiv_points and the polyscope view
void step(int n) {
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
void myCallback() {
  if (ImGui::Button("do work")) {
    step(100);
  }
}
void dbgOutput(string filename) {
    std::ofstream s(filename, std::ofstream::out);
    s << nVertices << " " << nEdges << endl;
    for (Vertex v : mesh->vertices()) {
        auto thing = flattened[vInd[v]];
        s << thing.x << " " << thing.y << "\n";
    }
    for (Edge e : mesh->edges()) {
        s << vInd[e.halfedge().vertex()] << " " << vInd[e.halfedge().twin().vertex()] << "\n";
    }
}
inline double shift(double c) {
    return (c + 2.) * 500;
}
void dbgSVG(string filename) {
    std::ofstream ss(filename, std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
       << "<svg width=\"2000\" height=\"2000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
    for (Vertex v : mesh->vertices()) {
        auto thing = flattened[vInd[v]];
        ss << "<circle cx=\"" << shift(thing.x) << "\" cy=\"" << shift(thing.y) << "\" r=\"1\"/>" << endl;
    }
    for (Edge e : mesh->edges()) {
        Vector2 i = flattened[vInd[e.halfedge().vertex()]];
        Vector2 j = flattened[vInd[e.halfedge().twin().vertex()]];
        double angle = alphas[eInd[e]];
        // FROM NORMALIZATION
        double radius = 500 * norm(i - j) / abs(2 * sin(angle));
        if (abs(angle) > 1e-7) {
            string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
            // sweep flag is 1 if going outward, 0 if going inward
            string sweepFlag = angle < 0 ? "0" : "1";
            ss << "<path d=\"M" << shift(i.x) << "," << shift(i.y) << " A" << radius << ","
            << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
            << shift(j.x) << "," << shift(j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
        } else {
            ss << "<line x1=\"" << shift(i.x) << "\" x2=\"" << shift(j.x) 
            << "\" y1=\"" << shift(i.y) <<"\" y2=\"" << shift(j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
        }
        
    }
    // footer
    ss << "</svg>";
}

double confObjective(vector<Vector2> &flattened, vector<double>& alphas) {
    double residual = 0.;
    bool first = true;
    for (Vertex v: mesh->vertices()) {
        double accum = 0.;
        double targetAccum = 0;
        //double vAngle = geometry->vertexAngleSums[v];
        for (Corner C : v.adjacentCorners()) {
            targetAccum += targetAngles[C];
            Halfedge h = C.halfedge();
            Halfedge ab = h;
            Halfedge bc = h.next();
            Halfedge ca = h.next().next();
            if (h.isInterior()) {
                const Vector2 &a = flattened[vInd[ab.vertex()]];
                const Vector2 &b = flattened[vInd[bc.vertex()]];
                const Vector2 &c = flattened[vInd[ca.vertex()]];
                auto u = b - a;
                auto v = c - a;
                double angle = orientedAngle(u,v);
                // alphas are offsets from canonical halfedge orientation
                if (h.edge().halfedge() == h) {
                    angle += alphas[eInd[h.edge()]];
                } else {
                    angle -= alphas[eInd[h.edge()]];
                }
                if (ca.edge().halfedge() == ca) {
                    angle += alphas[eInd[ca.edge()]];
                } else {
                    angle -= alphas[eInd[ca.edge()]];
                }

                accum += angle;
            }
            
        }
        if (first) {
            if (v.isBoundary()) {
                //cout << "Actual boundary sum: " << accum << "target: " << targetAccum << endl;
            } 
            first = false;
        }
    }
    //for (Corner C: mesh->corners()) cout << targetAngles[C] << endl;
    double result = 0.;
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        if (h.isInterior()) {
            const Vector2 &a = flattened[vInd[ab.vertex()]];
            const Vector2 &b = flattened[vInd[bc.vertex()]];
            const Vector2 &c = flattened[vInd[ca.vertex()]];
            auto u = b - a;
            auto v = c - a;
            double angle = orientedAngle(u,v);
            // alphas are offsets from canonical halfedge orientation
            if (h.edge().halfedge() == h) {
                angle += alphas[eInd[h.edge()]];
            } else {
                angle -= alphas[eInd[h.edge()]];
            }
            if (ca.edge().halfedge() == ca) {
                angle += alphas[eInd[ca.edge()]];
            } else {
                angle -= alphas[eInd[ca.edge()]];
            }
            //cout << "target:" << targetAngles[C] << "   " << "actual:" << angle << endl;
            // check: sign issue?
            if (ab.vertex().isBoundary()) {
                //cout << "BOUNDARY " << angle << " target " << targetAngles[C] << endl;
            }
            result += sqr(angle - targetAngles[C]);
        }
    }
    return result/2;
}
void confGradient(vector<Vector2> &x, vector<double> &alphas, vector<Vector2> &flattened_grad, vector<double>& alphas_grad) {
    for (int i = 0; i < nVertices; i++) {
        flattened_grad[i] = Vector2::zero();
    }
    fill(alphas_grad.begin(), alphas_grad.end(), 0);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        if (h.isInterior()) {
            const Vector2 &a = flattened[vInd[C.vertex()]];
            const Vector2 &b = flattened[vInd[bc.vertex()]];
            const Vector2 &c = flattened[vInd[ca.vertex()]];
            Vector2 u = b - a;
            Vector2 v = c - a;
            // this rotate is ccw
            Vector2 bGrad = -u.rotate90()/u.norm2();
            Vector2 cGrad = v.rotate90()/v.norm2();


            double angle = orientedAngle(u,v);
            //cout << "angle is" << angle << "\n";
            // alphas are offsets from canonical halfedge orientation
            if (h.edge().halfedge() == h) {
                angle += alphas[eInd[h.edge()]];
            } else {
                angle -= alphas[eInd[h.edge()]];
            }
            // fix second conditional in the objective
            if (ca.edge().halfedge() == ca) {
                angle += alphas[eInd[ca.edge()]];
            } else {
                angle -= alphas[eInd[ca.edge()]];
            }
            double diff = angle - targetAngles[C];
            flattened_grad[vInd[bc.vertex()]] += diff * bGrad;
            flattened_grad[vInd[ca.vertex()]] += diff * cGrad;
            // TODO double check this
            flattened_grad[vInd[C.vertex()]] -= diff*(bGrad + cGrad);
            if (h.edge().halfedge() == h) {
                alphas_grad[eInd[h.edge()]] += diff;
            } else {
                alphas_grad[eInd[h.edge()]] -= diff;
            }
            if (ca.edge().halfedge() == ca) {
                alphas_grad[eInd[ca.edge()]] += diff;
            } else {
                alphas_grad[eInd[ca.edge()]] -= diff;
            }
        }
    }
}
void confStep(int n) {
    //for (auto& v : subdiv_points) v *= 1.1;
    cout << "Starting planar optimization" << endl;
    vector<Vector2> flattened_new = flattened;
    vector<Vector2> flattened_grad(nVertices);
    vector<double> alphas_new = alphas;
    vector<double> alphas_grad(nEdges);
    double grad_size = 1.;

    for (int m = 0; m < n; m++){
    //while (sqrt(grad_size) > ep) {
        double result = confObjective(flattened, alphas);
        confGradient(flattened, alphas, flattened_grad, alphas_grad);
        grad_size = 0; // change this
        for (Vector2 a: flattened_grad) grad_size += a.norm2();
        for (double a: alphas_grad) grad_size += a*a;
        
        //double t = .00001;
        double t = 0.01;
        int steps = 0;
        while (confObjective(flattened_new, alphas_new) > result - alpha * t * grad_size) {
            steps++;
            t = beta * t;
            for (int i = 0; i < nVertices; i++) {
                flattened_new[i] = flattened[i] - t * flattened_grad[i];
                //cout << "flat" << flattened_grad[i].x << " " << flattened_grad[i].y << endl;
            }
            for (int i = 0; i < nEdges; i++) {
                alphas_new[i] = alphas[i] - t * alphas_grad[i];
                //cout << "a" << alphas_grad[i] << endl;
            }
            /*
            int i = 6;
            flattened_new[i] = flattened[i] + Vector2{t,0};
            //alphas_new[i] = alphas[i] + t;
            double new_result = confObjective(flattened_new, alphas_new);
            cout << "Objective change:" << new_result - result << endl;
            cout <<  "Estimated change:" << t *(flattened_grad[i].x) << endl;
            */
        }
        
        flattened = flattened_new;
        alphas = alphas_new;
        if (iter % 1000 == 0) {
            cout << "Starting iteration " << iter << endl;
            cout << "t: " << t << endl;
            cout << "grad" << grad_size << endl;
            cout << "steps: " << steps << endl;
            //cout << "grad size squared:" << grad_size << endl;
            cout << "objective:" << result << endl;
            dbgSVG("step" + std::to_string(iter + 1) + ".svg");
            //std::ofstream s ("chug" + std::to_string(iter), std::ofstream::out);
        }
        /*
        if (iter == 9999) {
            cout << "final";
            dbgOutput("final");
        }
        */
        iter++;
    }
    //cout << "iteration count: " << iter << endl;
    //cout << "grad size:" << sqrt(grad_size) << endl;
    //cout << "objective:" << objective(x1, x2, x3) << endl;
    //vector<Vector3> pos = subdiv_points;
}

void loadModel(const std::string& inputPath, bff::Model& model,
			   std::vector<bool>& surfaceIsClosed) {
	std::string error;
	if (bff::MeshIO::read(inputPath, model, error)) {
        bff::Mesh& mesh = model[0];
        int nBoundaries = (int)mesh.boundaries.size();
        if (nBoundaries >= 1) {
            // mesh has boundaries
            int eulerPlusBoundaries = mesh.eulerCharacteristic() + nBoundaries;
            if (eulerPlusBoundaries == 2) {
                // fill holes if mesh has more than 1 boundary
                if (nBoundaries > 1) {
                    if (bff::HoleFiller::fill(mesh)) {
                        // all holes were filled
                        surfaceIsClosed[0] = true;
                    }
                }
            } else {
                // mesh probably has holes and handles
                bff::HoleFiller::fill(mesh, true);
                bff::Generators::compute(mesh);
            }

        } else if (nBoundaries == 0) {
            if (mesh.eulerCharacteristic() == 2) {
                // mesh is closed
                surfaceIsClosed[0] = true;
            } else {
                // mesh has handles
                bff::Generators::compute(mesh);
            }
        }

	} else {
		std::cerr << "Unable to load file: " << inputPath << ". " << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

void flatten(bff::Model& model, const std::vector<bool>& surfaceIsClosed,
			 int nCones, bool flattenToDisk, bool mapToSphere) {
    bff::Mesh& mesh = model[0];
    bff::BFF bff(mesh);

    if (nCones > 0) {
        std::vector<bff::VertexIter> cones;
        bff::DenseMatrix coneAngles(bff.data->iN);
        int S = std::min(nCones, (int)mesh.vertices.size() - bff.data->bN);

        if (bff::ConePlacement::findConesAndPrescribeAngles(S, cones, coneAngles, bff.data, mesh)
            == bff::ConePlacement::ErrorCode::ok) {
            if (!surfaceIsClosed[0] || cones.size() > 0) {
                bff::Cutter::cut(cones, mesh);
                bff.flattenWithCones(coneAngles, true);
            }
        }
    } else {
        if (surfaceIsClosed[0]) {
                std::cerr << "Surface is closed. Either specify nCones or mapToSphere." << std::endl;
                exit(EXIT_FAILURE);
        } else {
                bff::DenseMatrix u(bff.data->bN);
                bff.flatten(u, true);
        }
    }
}


void conformalFlatten() {
	// parse command line options
	std::string inputPath = inputMeshPath;
	int nCones = 0;
	bool flattenToDisk = false;
	bool mapToSphere = false;
	bool normalizeUVs = false;
	// load model
	std::vector<bool> surfaceIsClosed(1);
	loadModel(inputPath, model, surfaceIsClosed);

	// set nCones to 8 for closed surfaces`
    if (surfaceIsClosed[0] && !mapToSphere && nCones < 3) {
        std::cout << "Setting nCones to 8." << std::endl;
        nCones = 8;
    }

    // flatten
    flatten(model, surfaceIsClosed, nCones, flattenToDisk, mapToSphere);
    for (bff::VertexCIter v = model[0].vertices.begin(); v != model[0].vertices.end(); v++) {
        flattened.push_back({v->wedge()->uv.x, v->wedge()->uv.y});
        //cout << v->wedge()->uv.x << endl;
        //cout << v->wedge()->uv.y << endl;
    }
    // DEBUG
    alphas = vector<double>(nEdges, 0);
    targetAngles = CornerData<double>(*mesh);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        if (h.isInterior()) {
            double angle = geometry->cornerAngle(C);
            angle += sol[eInd[h.edge()]];
            angle += sol[eInd[h.next().next().edge()]];
            targetAngles[C] = angle;
            //cout << "angle is" << angle;
        }
    }
    /*
    // DEBUG
    for (Vertex v: mesh->vertices()) {
        double accum = geometry->vertexAngleSums[v];
        for (Edge e: v.adjacentEdges()) {
            accum += (e.isBoundary() ? 1 : 2)*sol[eInd[e]];
        }
        if (v.isBoundary()) {
            cout << "boundary from solve" << accum << endl;
        } else {
            cout << "not boundary from solve" << accum << endl;
        }
    }
    
    // DEBUG
    for (Vertex v: mesh->vertices()) {
        double accum = 0;
        //double vAngle = geometry->vertexAngleSums[v];
        for (Corner C : v.adjacentCorners()) {
            accum += targetAngles[C];
        }
        if (v.isBoundary()) {
            cout << "boundary" << accum << endl;
        } else {
           cout << "not boundary" << accum << endl;
        }
    }
    */
    dbgSVG("step0.svg");
    //dbgOutput("chug0");
    confStep(50000);
    /*
    for (bff::WedgeCIter w: model[0].cutBoundary()) {
            w->uv;
    }
    */        

}
// just for validating the SVG formula I'm using
void testSVG() {
    std::ofstream ss("test.svg", std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
       << "<svg width=\"2000\" height=\"2000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
    ss << "<circle cx=\"" << 500 << "\" cy=\"" << (500) << "\" r=\"1\"/>" << endl;
    ss << "<circle cx=\"" << 500 << "\" cy=\"" << (1000) << "\" r=\"1\"/>" << endl;

    Vector2 i = {500, 500};
    Vector2 j = {1000, 500};
    double angle = 3*PI/4;
    double radius = norm(i - j) / abs(2 * sin(angle));
    string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
    // sweep flag is 1 if going outward, 0 if going inward
    string sweepFlag = angle < 0 ? "0" : "1";
    ss << "<path d=\"M" << (i.x) << "," << (i.y) << " A" << radius << ","
       << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
       << (j.x) << "," << (j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
    ss << "<line x1=\"" << (i.x) << "\" x2=\"" << (j.x)
       << "\" y1=\"" << (i.y) << "\" y2=\"" << (j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
    
    double wog = 1.3;
    j = {500 + 700*cos(wog), 500 + 800*sin(wog)};
    angle = -(PI - wog - angle);
    radius = norm(i - j) / abs(2 * sin(angle));
    largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
    // sweep flag is 1 if going outward, 0 if going inward
    sweepFlag = angle < 0 ? "0" : "1";
    ss << "<path d=\"M" << (i.x) << "," << (i.y) << " A" << radius << ","
       << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
       << (j.x) << "," << (j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
    ss << "<line x1=\"" << (i.x) << "\" x2=\"" << (j.x)
       << "\" y1=\"" << (i.y) << "\" y2=\"" << (j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
    ss << "</svg>";
}
int main(int argc, char **argv) {

    // Configure the argument parser
    args::ArgumentParser parser("");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a mesh name was given
    if (!inputFilename) {
        //inputMeshPath = "/home/elu/repos/catopt/meshes/spotwithhole.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/cube.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
        //std::cerr << "Please specify a mesh file as argument" << std::endl;
        //return EXIT_FAILURE;
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    // Initialize polyscope
    cout << "Initialized" << endl;
    //polyscope::init("openGL_mock");
    polyscope::init();

    polyscope::state::userCallback = myCallback;
    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        "main",
        geometry->inputVertexPositions, mesh->getFaceVertexList(),
        polyscopePermutations(*mesh));
    cout << "starting optimization" << endl;
    initializeQuantities();
    generateConstraints();

    subdivision();
    buildNewMesh();

    //fin = descent();

    //////////////////////////////////////////
    // BFF stuff
    conformalFlatten();
    //testSVG();

    //generateVisualization();
    // Give control to the polyscope gui
    //polyscope::show();

    return EXIT_SUCCESS;
}
