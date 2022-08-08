#include "args/args.hxx"
#include <string>
#include <iostream>
#include <memory>


#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
using std::string;
using std::shared_ptr;

#include "Common.h"


#include "IntrinsicFlattening.h"
#include "EmbeddingOptimization.h"
#include "CircleWrapper.h"

#include "nasoq/nasoq_eigen.h"
typedef Eigen::Triplet<double> T;

shared_ptr<ManifoldSurfaceMesh> mesh;
shared_ptr<VertexPositionGeometry> geometry;
polyscope::SurfaceMesh *psMesh;

void planarMapping(int N) {
    for (int m = 0; m <= N; m++){
        IntrinsicFlattening flattener(mesh, geometry);
        SolutionData sol = flattener.solveFromPlane((double)m/(double)N);
        VertexData<double> vertexSum(*mesh);
        for (Vertex v: mesh->vertices()) {
            double accum = 0;
            for (Corner c: v.adjacentCorners()) {
                accum += sol.betas[c];
            }
            vertexSum[v] = accum;
        }
        /*
           psMesh->addVertexScalarQuantity("vertex sums", vertexSum);
        //psMesh->addHalfedgeScalarQuantity("solved", sol.betas);
        psMesh->addEdgeScalarQuantity("thetas", sol.thetas);
        //psMesh->addVertexScalarQuantity("infinite vertex", sol.infVertex);
        psMesh->addFaceScalarQuantity("faces", sol.fMask);
        */

        CircleWrapper patterns(mesh, sol, psMesh);
        patterns.solve("fin" + std::string(3 - std::to_string(m).length(), '0') +  std::to_string(m) );
    }
}

EmbeddingOptimization* E;
float t = 1e-4;
void embedding(int N) {
    IntrinsicFlattening flattener(mesh, geometry);
    CornerData<double> beta = flattener.solveIntrinsicOnly();
    E = (new EmbeddingOptimization(mesh, geometry, beta));
    auto [submesh, subgeometry] = E->solve(N);

}
void surfaceToPlane() {
    IntrinsicFlattening flattener(mesh, geometry);
    EdgeData<double> intersectionAngles = flattener.solveKSS();
    CircleWrapper patterns(mesh, intersectionAngles, psMesh);
    patterns.solveKSS();
}
void myCallback() {

    // Since options::openImGuiWindowForUserCallback == true by default, 
    // we can immediately start using ImGui commands to build a UI


    ImGui::InputFloat("param value", &t, 0.01f, 1.0f);  // set a float variable

    if (ImGui::Button("run subroutine")) {
        E->optimize(t);
    }
}
void nasoqTest() {
    vector<T> HList = {T(0,0,1), T(1,1,1)};
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> H(2,2);
    H.setFromTriplets(HList.begin(), HList.end());

    vector<T> AList;
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> A(1,2); 
    A.setFromTriplets(AList.begin(), AList.end());

    vector<T> CList = {T(1,1,-1)};
    Eigen::SparseMatrix<double,Eigen::ColMajor,int> C(2,2); 
    C.setFromTriplets(CList.begin(), CList.end());

    Eigen::VectorXd q(2);
    q[0] = 0; q[1] = 0;
    Eigen::VectorXd b(1);
    b[0] = 0;
    Eigen::VectorXd d(2);
    d[1] = -2;



    /// New settings if provided
    int iter;
    std::string nasoq_mode;
    double pert, eps, tol;
    nasoq::QPSettings *qs = NULL;

    /// output vectors
    Eigen::VectorXd x, y, z;


    /// call the wrapper.
    int ret = nasoq::quadprog(H,q,A,b,C,d,x,y,z,qs);
    for (int i = 0; i < 2; i++) {
        cout << "wah" << x[i] << endl;
    }
    
}
int main(int argc, char **argv) {
    // Configure the argument parser
    string inputMeshPath;
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
        inputMeshPath = "/home/elu/repos/catopt/meshes/spotwithhole.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/SmallDisk.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/BumpyTorusPatch.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/plane.obj";
        inputMeshPath = "/home/elu/repos/catopt/meshes/beanhole.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/nonconvex2.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/test.obj";
        //inputMeshPath = "/home/elu/repos/catopt/meshes/patch.obj";
    } else {
        inputMeshPath = args::get(inputFilename);
    }
    std::tie(mesh, geometry) = readManifoldSurfaceMesh(inputMeshPath);
    // polyscope sanity checks
    polyscope::init();
    /*
       psMesh = polyscope::registerSurfaceMesh(
       "waaah",
       geometry->inputVertexPositions, mesh->getFaceVertexList(),
       polyscopePermutations(*mesh));
       */

    mesh->compress();
    polyscope::state::userCallback = myCallback;
    //planarMapping(100);
    //embedding(5);

    //polyscope::show();
    nasoqTest();

    return EXIT_SUCCESS;
}
