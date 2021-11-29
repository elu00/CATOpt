#pragma once
#include "Common.h"

class CircleWrapper {
    public:
        CircleWrapper(shared_ptr<ManifoldSurfaceMesh> mesh, SolutionData sol);
        void solve();
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

        void setOffsets();
};