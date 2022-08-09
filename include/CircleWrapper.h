#pragma once
#include "Common.h"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

class CircleWrapper {
    public:
        CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, CornerData<double> betas, polyscope::SurfaceMesh *psMesh);
        CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> intersectionAngles, polyscope::SurfaceMesh *psMesh);
        void solve(std::string name = "fin");
        void uvSVG(std::string filename);
    private:
        // geometric data
        shared_ptr<ManifoldSurfaceMesh> mesh;
        CornerData<double> beta; 
        EdgeData<double> theta;
        VertexData<Eigen::Vector2d> uv;
        Vector<double> circleSol;
        Eigen::VectorXd thetas;
        polyscope::SurfaceMesh *psMesh;

        void setOffsets();
};
