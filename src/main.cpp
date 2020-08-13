#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"

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
// SparseMatrix<double> constraints;
// Vector<double> x_init;
vector<double> sol;
vector<double> rhs;
vector<double> ineqRHS0;
vector<double> ineqRHS1;

size_t iter = 0;

size_t subdiv_level = 5;
vector<Vector3> subdiv_points;
//vector<Vector3> grad;
vector<tuple<size_t, size_t, double>> correct_dist;
Eigen::SparseMatrix<double> bendingMatrix;
double alpha = 0.1;
double beta = 0.5;
double ep = 1e-5;
double bendingWeight = 1e-8;

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
    angleDefects = geometry->vertexGaussianCurvatures;
    sol = vector<double>(nEdges);
}

void generateConstraints() {
    // Model initialization
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", nEdges, Domain::inRange(-2 * PI, 2 * PI));

    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    // Equality constraint initialization
    rhs = vector<double>(nVertices);
    for (size_t i = 0; i < nVertices; i++) {
        rhs[i] = angleDefects[mesh->vertex(i)] / 2.;
    }
    for (Edge e : mesh->edges()) {
        rows.emplace_back(vInd[e.halfedge().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(1.);
        rows.emplace_back(vInd[e.halfedge().twin().vertex()]);
        cols.emplace_back(eInd[e]);
        values.emplace_back(1.);
    }
    auto r = new_array_ptr<int>(rows);
    auto c = new_array_ptr<int>(cols);
    auto v = new_array_ptr<double>(values);
    auto Meq = Matrix::sparse(nVertices, nEdges, r, c, v);
    auto eqRHS = new_array_ptr(rhs);
    M->constraint("eq constraints", Expr::mul(Meq, x), Domain::equalsTo(eqRHS));
    cout << "eq generated" << endl;
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
    cout << "ineq generated" << endl;
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::vstack(t, x), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << x->level() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = x->getSize();
    auto xVal = x->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    for (int i = 0; i < xsize; ++i) {
        (sol)[i] = (*xVal)[i];
    }
    return;
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


        // DEBUGGGG
        std::ofstream s ("tri" + std::to_string(index) + ".svg", std::ofstream::out);
        s << CAT(ij*500, jk*500, ki*500, a_ij, a_jk, a_ki, true).triangulation_svg(subdiv_level);
        s.close();
        int which_f = fInd[f], which_e;



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
        which_e = eInd[it.edge()];
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
        which_e = eInd[it.edge()];
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
        which_e = eInd[it.edge()];
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

    cout << "Laplacian: " << L.norm();
    Eigen::SparseMatrix<double> M = intrinsicGeometry->vertexLumpedMassMatrix.cwiseInverse();
    cout << "Mass: " << M.norm();
    bendingMatrix = L.transpose() * M * L;
    cout << "bending matrix made" << endl;

}

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
tuple<VectorXd, VectorXd, VectorXd> gradient(const VectorXd &x1, const VectorXd &x2, const VectorXd &x3) {
    int i1, i2;
    double distsq, actualDistsq, diff;
    VectorXd grad1 = VectorXd::Zero(subdiv_points.size());
    VectorXd grad2 = VectorXd::Zero(subdiv_points.size());
    VectorXd grad3 = VectorXd::Zero(subdiv_points.size());
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

double grad_norm_sq(const VectorXd &grad1, const VectorXd &grad2, const VectorXd &grad3) {
    return grad1.squaredNorm() + grad2.squaredNorm() + grad3.squaredNorm();
}

void generateVisualization() {
    // Visualization
    psMesh->addEdgeScalarQuantity("Final Solution", sol);
    psMesh->addVertexScalarQuantity("curvature",
                                    geometry->vertexGaussianCurvatures);
    cout << "No of subdivision points:" << subdiv_points.size();
}

void step(int n) {
    //for (auto& v : subdiv_points) v *= 1.1;
    cout << "Starting descent" << endl;
    VectorXd grad1;
    VectorXd grad2;
    VectorXd grad3;
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
    //for (int m = 0; m < n; m++){
    while (sqrt(grad_size) > ep) {
        double result = objective(x1, x2, x3);
        std::tie(grad1, grad2, grad3) = gradient(x1, x2, x3);
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
    for (int i = 0; i < subdiv_points.size(); i++) {
        subdiv_points[i].x = x1[i];
        subdiv_points[i].y = x2[i];
        subdiv_points[i].z = x3[i];
    }
    //return pos;
    CATpsMesh->updateVertexPositions(subdiv_points);
}
void myCallback() {
  if (ImGui::Button("do work")) {
    step(100);
  }
}
int main(int argc, char **argv) {

    /*
    // Configure the argument parser
    args::ArgumentParser parser("geometry-central & Polyscope example project");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a mesh name was given
    if (!inputFilename)
    {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }
    */
    // Initialize polyscope
    cout << "Initialized" << endl;
    //polyscope::init("openGL_mock");
    polyscope::init();

    polyscope::state::userCallback = myCallback;
    // Load mesh
    std::tie(mesh, geometry) = readManifoldSurfaceMesh("/home/elu/repos/CATOpt/meshes/spot.obj");

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

    //generateVisualization();
    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
