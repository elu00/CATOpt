#pragma once

#include <tuple>
#include <map>
#include <iostream>
#include <vector>
#include <queue>
#include <memory>
#include <tuple>
using std::cout;
using std::vector;
using std::shared_ptr;
using std::endl;
using std::tuple;


#include "Common.h"
#include "Common.h"
#include "EdgeLengthOptimization.h"
#include "IntrinsicFlattening.h"

std::pair<EdgeData<double>, CornerData<double>> PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, 
EdgeData<double> l, CornerData<double> beta, 
        VertexData<double> VertexCurvatures, EdgeData<double> EdgeCurvatures) ;
CornerData<Vector2> LayoutMesh(shared_ptr<ManifoldSurfaceMesh> mesh, 
        EdgeData<double> l, EdgeData<bool> S);
void CATToSVG(shared_ptr<ManifoldSurfaceMesh> mesh, CornerData<Vector2> p,
        CornerData<double> beta, std::string filename);



