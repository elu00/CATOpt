#pragma once
#include "Common.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

class CircleWrapper {
    public:
        CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, SolutionData sol, polyscope::SurfaceMesh *psMesh);
        CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> intersectionAngles, polyscope::SurfaceMesh *psMesh);
        void solve(std::string name = "fin");
        void solveKSS();
        void uvSVG(std::string filename, EdgeData<bool> eMask);
    private:
        // geometric data
        shared_ptr<ManifoldSurfaceMesh> mesh;
        CornerData<double> beta; 
        EdgeData<double> theta;
        VertexData<Eigen::Vector2d> uv;
        Eigen::Vector2d center;
        double invRadius;
        void circleInversion();
        Vector<double> circleSol;
        Vertex infVertex;
        EdgeData<bool> eMask; 
        EdgeData<bool> eBdry; 
        FaceData<bool> fMask;
        Eigen::VectorXd thetas;
        polyscope::SurfaceMesh *psMesh;

        void setOffsets();
};