#include "geometrycentral/surface/halfedge_mesh.h"
#include "geometrycentral/surface/meshio.h"
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
using std::unique_ptr;

using namespace mosek::fusion;
using namespace monty;
using std::vector;
using std::tuple;
using std::make_tuple;
typedef std::tuple<Face, size_t, size_t> bindex;

// == Geometry-central data
unique_ptr<HalfedgeMesh> mesh;
unique_ptr<VertexPositionGeometry> geometry;

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

size_t subdiv_level = 8;
vector<Vector3> subdiv_points;
//vector<Vector3> grad;
std::map<bindex, size_t> indexing;
vector<tuple<size_t, size_t, double>> correct_dist;
double alpha = 0.2;
double beta = 0.5;
double ep = 1e-5;
vector<Vector3> fin;

// Polyscope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh *psMesh;

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
    Model::t M = new Model(); auto _M = finally([&]() { M->dispose(); });
    Variable::t x = M->variable("x", nEdges, Domain::inRange(-2*PI, 2 * PI));

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
        ineqRHS0[cInd[c]] = - geometry->cornerAngle(c);
        ineqRHS1[cInd[c]] = 2*PI - geometry->cornerAngle(c);
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
    // sum is less thatn 2pi
    M->constraint("ineq1 constraints", Expr::mul(Mineq, x), Domain::lessThan(inRHS1));

    auto ones =  std::make_shared<ndarray<double,1>>(shape(nEdges),1.);
    cout << "ineq generated" << endl;
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint( Expr::vstack(t, x), Domain::inQCone() );
    M->objective( ObjectiveSense::Minimize, t );
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << x->level() << endl;
    cout << "Optimization Done" << endl;
    auto xsize = x->getSize();
    auto xVal = x->level();
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    for(int i = 0; i < xsize; ++i) {
        (sol)[i] =  (*xVal)[i];
    }
    return;
}


inline Vector3 bary (Face f, double a, double b, double c) {
    //cout << a + b + c << endl;
    auto it = f.halfedge();
    Vector3 i = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 j = geometry->inputVertexPositions[it.vertex()];
    it = it.next();
    Vector3 k = geometry->inputVertexPositions[it.vertex()];
    return a*i + b*j + c*k;
}

//returns -1 if invalid index
inline size_t get_index(Face f, int j, int k) {
    auto key = make_tuple(f, j, k);
    auto i = indexing.find(key);
    return i == indexing.end() ? -1 : i->second;
}
inline double i_coord(size_t j, size_t k) {
    return 1.-((double)(j+k))/subdiv_level;
}
inline double to_bary(size_t i) {
    return ((double)i)/subdiv_level;
}
void subdivision() {
    EdgeData<size_t> eStart(*mesh, -1);
    VertexData<size_t> v(*mesh, -1);
    // initializing indexing map because I'm lazy
    cout << "starting subdivision" << endl;
    size_t index = 0;
    // initialization of initial guesses
    for (Face f: mesh->faces()) {
        auto it = f.halfedge();

        Vertex i = it.vertex();
        auto T = make_tuple(f, 0, 0);
        if (v[i] == -1) {
            v[i] = index;
            indexing[T] = index;
            index++;
        } else {
            indexing[T] = v[i];
        }
        it = it.next();

        Vertex j = it.vertex();
        T = make_tuple(f, subdiv_level, 0);
        if (v[j] == -1) {
            v[j] = index;
            indexing[T] = index;
            index++;
        } else {
            indexing[T] = v[j];
        }
        it = it.next();

        Vertex k = it.vertex();
        T = make_tuple(f, 0, subdiv_level);
        if (v[k] == -1) {
            v[k] = index;
            indexing[T] = index;
            index++;
        } else {
            indexing[T] = v[k];
        }
        it = it.next();

        // edge ij
        size_t start;
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
                T = make_tuple(f, pos, 0);
                indexing[T] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                T = make_tuple(f, subdiv_level - pos, 0);
                indexing[T] = start;
                start++;
            }

        }
        it = it.next();

        // edge jk
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
                T = make_tuple(f, subdiv_level - pos, pos);
                indexing[T] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                T = make_tuple(f, pos, subdiv_level - pos);
                indexing[T] = start;
                start++;
            }

        }
        it = it.next();

        // edge ki
        if (eStart[it.edge()] == -1) {
            eStart[it.edge()] = index;
            start = index;
            index += subdiv_level - 1;
        } else {
            start = eStart[it.edge()];
        }
        if (it.edge().halfedge() == it) {
            // k is canonical vertex 
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                T = make_tuple(f, 0, subdiv_level - pos);
                indexing[T] = start;
                start++;
            }
        } else {
            for (size_t pos = 1; pos < subdiv_level; pos++) {
                T = make_tuple(f, 0, pos);
                indexing[T] = start;
                start++;
            }

        }
        it = it.next();

        // interior points
        for (size_t j = 1; j < subdiv_level; j++) {
            for (size_t k = 1; k < subdiv_level - j; k++) {
                T = make_tuple(f, j, k);
                indexing[T] = index;
                index++;
            }
        }

    }
    // initialize vector
    subdiv_points = vector<Vector3>(index, Vector3{0,0,0});
    for (Face f: mesh->faces()) {
        for (size_t j = 0; j <= subdiv_level; j++) {
            for (size_t k = 0; k <= subdiv_level - j; k++) {
                if (get_index(f, j, k) == 59) {
                    cout << "59 found:" << endl;
                    cout << fInd[f] << endl << j << endl << k << endl;
                }
                subdiv_points[get_index(f, j, k)] = bary(f,i_coord(j, k), to_bary(j) , to_bary(k));
            }
        }
    }

    // initialization of "correct" lengths
    for (Face f: mesh->faces())
    {
        size_t j_, k_;
        size_t index1, index2;
        double dist;
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
        // DEBUG
        a_ij = 0;
        a_jk = 0;
        a_ki = 0;

        for (size_t j = 0; j <= subdiv_level; j++) {
            for (size_t k = 0; k <= subdiv_level - j; k++) {
                index1 = get_index(f, j, k);
                // right and up 1
                j_ = j;
                k_ = k + 1;
                index2 = get_index(f, j_, k_);
                if (index2 != -1) {
                    dist = l2DistSquared(ij, jk, ki, a_ij, a_jk, a_ki, i_coord(j,k), to_bary(j), to_bary(k), i_coord(j_, k_), to_bary(j_), to_bary(k_));
                    /*
                    cout << "j:" << j << "k:" << k << endl;
                   cout << "j_:" << j_ << "k_:" << k_ << endl;
                       cout << dist << endl;
                       */
                    correct_dist.push_back(make_tuple(index1, index2, dist));
                }
                // right 1
                j_ = j + 1;
                k_ = k;
                index2 = get_index(f, j_, k_);
                if (index2 != -1) {
                    dist = l2DistSquared(ij, jk, ki, a_ij, a_jk, a_ki, i_coord(j,k), to_bary(j), to_bary(k), i_coord(j_, k_), to_bary(j_), to_bary(k_));
                    /*
                       cout << "j:" << j << "k:" << k << endl;
                       cout << "j_:" << j_ << "k_:" << k_ << endl;
                       cout << dist << endl;
                       */
                    correct_dist.push_back(make_tuple(index1, index2, dist));
                }
                // right and down
                j_ = j + 1;
                k_ = k - 1;
                index2 = get_index(f, j_, k_);
                if (index2 != -1) {
                    dist = l2DistSquared(ij, jk, ki, a_ij, a_jk, a_ki, i_coord(j,k), to_bary(j), to_bary(k), i_coord(j_, k_), to_bary(j_), to_bary(k_));
                    /*
                       cout << "j:" << j << "k:" << k << endl;
                       cout << "j_:" << j_ << "k_:" << k_ << endl;
                       cout << dist << endl;
                       */
                    correct_dist.push_back(make_tuple(index1, index2, dist));
                }
            }
        }
    }

}
/*
 * Optimization code
 */
double objective(vector<Vector3> x) {
    double result = 0.;
    int i1, i2;
    double distsq, diff;
    for (auto t: correct_dist) {
        i1 = std::get<0>(t);
        i2 = std::get<1>(t);
        distsq = std::get<2>(t);
        // update objective
        diff = (norm2(x[i1] - x[i2]) - distsq);
        result += diff * diff;
    }
    return result;
}
vector<Vector3> gradient(vector<Vector3> x) {
    int i1, i2;
    double distsq, diff;
    //zero out gradient
    vector<Vector3> grad (subdiv_points.size(), Vector3{0.,0.,0.});
    for (auto t: correct_dist) {
        i1 = std::get<0>(t);
        i2 = std::get<1>(t);
        distsq = std::get<2>(t);
        diff = (norm2(x[i1] - x[i2]) - distsq);
        // update gradient
        grad[i1] += 4. * diff * (x[i1] - x[i2]);
        grad[i2] += 4. * diff * (x[i2] - x[i1]);
    }
    return grad;
}

double grad_norm_sq(vector<Vector3> grad) {
    double result = 0.;
    for (auto v: grad) {
        result += norm2(v);
    }
    return result;
}

vector<Vector3> descent() {
    cout << "Starting descent" << endl;
    size_t iter = 0;
    double result;
    double t = 1.0;
    vector<Vector3> grad (subdiv_points.size(), Vector3{1.,0.,0.});
    double grad_size = 1.;
    vector<Vector3> x = subdiv_points;
    vector<Vector3> x_new = subdiv_points;
    while (sqrt(grad_size) > ep) {
        double result = objective(x);
        grad = gradient(x);
        grad_size = grad_norm_sq(grad);
        t = 1.;
        for (int i = 0; i < subdiv_points.size(); i++) {
            x_new[i] = x[i] - t * grad[i];
        }
        while (objective(x_new) > result - alpha * t * grad_size) {
            t = beta * t;
            for (int i = 0; i < subdiv_points.size(); i++) {
                x_new[i] = x[i] - t * grad[i];
            }
        }
        x = x_new;
        if (iter % 1000 == 0) {
            cout << "Starting iteration " << iter << endl;
            cout << "grad size squared:" << grad_size << endl;
            cout << "objective:" << result << endl;
        }
        iter++;
    }
    cout << "iteration count: " << iter << endl;
    cout << "grad size:" << sqrt(grad_size) << endl;
    cout << "objective:" << objective(x) << endl;
    return x;
}
void generateVisualization() {
    // Visualization
    //EdgeData<double> initialGuess(*mesh);
    //EdgeData<double> finalSolution(*mesh);
    //for (Edge e : mesh->edges())
    //{
    //    initialGuess[e] = x_init[eInd[e]];
    //    finalSolution[e] = x[eInd[e]];
    //cout << x[eInd[e]] << endl;
    //cout << x_init[eInd[e]] << endl;
    //}
    //psMesh->addEdgeScalarQuantity("Initial Guess", x_init);
    psMesh->addEdgeScalarQuantity("Final Solution", sol);
    psMesh->addVertexScalarQuantity("curvature",
            geometry->vertexGaussianCurvatures);
    cout << "No of subdivision points:" << subdiv_points.size();

    polyscope::registerPointCloud("Initial Subdivision", subdiv_points);
    polyscope::registerPointCloud("Final Subdivision", fin);
}
void generateSVGs() {
    int n = 0;
    for (Face f: mesh->faces()) {
        double ij, jk, ki, a_ij, a_jk, a_ki;
        Halfedge h = f.halfedge();
        ij = geometry->edgeLength(h.edge()) * 100;
        a_ij = (sol)[eInd[h.edge()]];
        h = h.next();
        jk = geometry->edgeLength(h.edge()) * 100;
        a_jk = (sol)[eInd[h.edge()]];
        h = h.next();
        ki = geometry->edgeLength(h.edge()) * 100;
        a_ki = (sol)[eInd[h.edge()]];
        CAT c = CAT(ij, jk, ki, a_ij, a_jk, a_ki);
        cout << ij << endl;
        cout << c.a_ij << endl;
        c.to_svg("current" + std::to_string(n) + ".svg");
        n++;
    }
}


int main(int argc, char **argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Optimization");
    args::Positional<std::string> inputFilename(parser, "mesh", "A mesh file.");

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    // Make sure a mesh name was given
    if (!inputFilename) {
        std::cerr << "Please specify a mesh file as argument" << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize polyscope
    cout << "Initialized" << endl;
    polyscope::init();

    // Load mesh
    std::tie(mesh, geometry) = loadMesh(args::get(inputFilename));

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
            polyscope::guessNiceNameFromPath(args::get(inputFilename)),
            geometry->inputVertexPositions, mesh->getFaceVertexList(),
            polyscopePermutations(*mesh));

    cout << "starting optimization" << endl;
    initializeQuantities();
    generateConstraints();
    subdivision();

    vector<Vector3> x = subdiv_points;
    cout << "OBJECTIVE" << objective(x) << endl;
    //fin = descent();
    generateVisualization();
    //generateSVGs();
    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
