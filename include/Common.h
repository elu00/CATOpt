#pragma once
#include <memory>
#include <cmath>
#include <iostream>
#include "glm/vec3.hpp"

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/numerical/linear_solvers.h"
using std::shared_ptr;
using std::cout;
using std::endl;
using namespace geometrycentral;
using namespace geometrycentral::surface;

struct SolutionData {

   // for Riemann mapping only:
    Vertex infVertex; // vertex at infinity
    // marks edges/faces near vertex at infinity
    EdgeData<bool> eMask;
    EdgeData<bool> eBdry; 
    FaceData<bool> fMask;
    CornerData<double> betas; // corner angles of circular arc triangulation
    EdgeData<double> thetas; // circumcircle intersection angles


};

struct BezierTriangle {
    Vector2 p200;
    Vector2 p020;
    Vector2 p002;
    Vector2 p110;
    Vector2 p011;
    Vector2 p101;
    double w200;
    double w020;
    double w002;
    double w110;
    double w011;
    double w101;
};
