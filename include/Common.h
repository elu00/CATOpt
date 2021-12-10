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
    Vertex infVertex;
    EdgeData<bool> eMask; 
    EdgeData<bool> eBdry; 
    FaceData<bool> fMask;
    CornerData<double> betas;
    EdgeData<double> thetas;
    
};