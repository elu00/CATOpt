#pragma once

#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <tuple>
#include <vector>
using std::cout;
using std::endl;
using std::shared_ptr;
using std::tuple;
using std::vector;

#include "Common.h"
#include "EdgeLengthOptimization.h"
#include "IntrinsicFlattening.h"

std::pair<EdgeData<double>, CornerData<double>>
PrescribeCurvature(shared_ptr<ManifoldSurfaceMesh> mesh, EdgeData<double> l,
                   CornerData<double> beta, VertexData<double> VertexCurvatures,
                   EdgeData<double> EdgeCurvatures);
CornerData<Vector2> LayoutMesh(shared_ptr<ManifoldSurfaceMesh> mesh,
                               EdgeData<double> l, EdgeData<bool> S);
void CATToSVG(shared_ptr<ManifoldSurfaceMesh> mesh, CornerData<Vector2> p,
              CornerData<double> beta, std::string filename);
