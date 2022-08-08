#include "IntrinsicFlattening.h"



IntrinsicFlattening::IntrinsicFlattening(shared_ptr<ManifoldSurfaceMesh> mesh,shared_ptr<VertexPositionGeometry> geometry):
    mesh(mesh), geometry(geometry) {
    geometry->requireEdgeLengths();
    geometry->requireVertexGaussianCurvatures();
    geometry->requireVertexAngleSums();
    nVertices = mesh->nVertices();
    nEdges = mesh->nEdges();
    nCorners = mesh->nCorners();
    nFaces = mesh->nFaces();
    c_ = mesh->getCornerIndices();
    e_ = mesh->getEdgeIndices();
    v_ = mesh->getVertexIndices();
    
    f_ = mesh->getFaceIndices();
    /* 
    angleDefects = geometry->vertexGaussianCurvatures;
    */

}


EdgeData<double> IntrinsicFlattening::solveKSS() {
    /*
    // Model initialization
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });
    // the witness of our "coherent angle system"
    Variable::t a = M->variable("a", nCorners, Domain::inRange(0, 2 * PI));

    // dummy vectors for building matrices
    vector<int> rows;
    vector<int> cols;
    vector<double> values;
    vector<double> rhs(0);

    // Flat Boundary Constraint: Total geodesic curvature is 0 on "most" of the boundary
    size_t excl = 4;
    size_t count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        if (count >= excl) {
            for (Corner c : v.adjacentCorners()) {
                rows.emplace_back(count);
                cols.emplace_back(c_[c]);
                values.emplace_back(1.);
            }
            rhs.push_back(PI);
        } else {
            rhs.push_back(0);
        }
        count++;
    }
    auto bdryFlat = sMatrix(count, nCorners, rows, cols, values);
    auto bdryPI = new_array_ptr(rhs);
    //M->constraint("Flat boundary", Expr::mul(bdryFlat, a), Domain::equalsTo(bdryPI));
    rhs.clear();

    FaceData<bool> fMask(*mesh, true);
    // halfedges that are now exterior to the mesh after deleting the vertex at infinity
    EdgeData<bool> eBdry(*mesh, false);
    
    // mark all the boundary edges as the new boundary
    for (Edge e : mesh->edges()) {
        eBdry[e] = e.isBoundary();
    }
    vector<double> origAngles(nCorners+1,0);
    for (Corner c: mesh->corners()) {
        origAngles[c_[c] + 1] = geometry->cornerAngle(c);
    }

    // just setup the objective: ||a - original angles||_{L^2}
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    M->constraint(Expr::sub(Expr::vstack(t, a),new_array_ptr(origAngles)), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << "Optimization Done" << endl;
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;
    EdgeData<double> thetaSolve(*mesh,0);

    auto asize = a->getSize();
    auto asol = a->level();
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior() && e.halfedge().isInterior()
                        ? (*asol)[c_[e.halfedge().next().next().corner()]]
                        : 0;
        double a2 =
            e.halfedge().twin().isInterior() && e.halfedge().twin().isInterior()
                ? (*asol)[c_[e.halfedge().twin().next().next().corner()]]
                : 0;
        thetaSolve[e] = PI - a1 - a2;
    }
    */
    EdgeData<double> thetaSolve(*mesh,0);
    return thetaSolve;
}

/*
void IntrinsicFlattening::buildIntersectionAngleConstraints(Model::t& M, CornerData<double>& beta, Variable::t& a) {
    // =========== Equation [1] =================================================
   // For each interior edge, add a constraint which ensures that
   // the coherent angle system exhibits the same circumcircle
   // intersection angles as the input CAT, namely
   //
   //   α0 + α1 = (β0 + β1 - β2 - β3 - β4 - β5)/2 + π       [1]
   //
   // where α are the angles in the coherent angle system, and β are the
   // CAT corner angles, and corners are indexed as below:
   //
   //        *
   //       /0\
   //      /   \
   //     /2   3\
   //    *-------*
   //     \5   4/
   //      \   /
   //       \1/
   //        *
   //

   // temporary arrays for building sparse linear system for Equation [1]
   vector<int> aRows;
   vector<int> aCols;
   vector<double> aVals;
   vector<double> aRhs = vector<double>(nEdges,0.); // right hand side of Equation [1]

   for (Edge e : mesh->edges()) {
      size_t ind = e_[e];

      // only preserve intersection angles on interior edges
      if (!e.isBoundary()) {

         // right-hand side ((β0 + β1 - β2 - β3 - β4 - β5)/2 + π)
         Halfedge h23 = e.halfedge();
         Halfedge h30 = h23.next();
         Halfedge h02 = h30.next();
         Halfedge h45 = h23.twin();
         Halfedge h51 = h45.next();
         Halfedge h14 = h51.next();

         Corner c0 = h02.corner();
         Corner c1 = h14.corner();
         Corner c2 = h23.corner();
         Corner c3 = h30.corner();
         Corner c4 = h45.corner();
         Corner c5 = h51.corner();

         aRhs[ind] = PI;
         aRhs[ind] += beta[c0]/2.;
         aRhs[ind] += beta[c1]/2.;
         aRhs[ind] -= beta[c2]/2.;
         aRhs[ind] -= beta[c3]/2.;
         aRhs[ind] -= beta[c4]/2.;
         aRhs[ind] -= beta[c5]/2.;

         // left-hand side (α0 + α1)
         aRows.push_back(ind);
         aCols.push_back(c_[c0]);
         aVals.push_back(1.);

         aRows.push_back(ind);
         aCols.push_back(c_[c1]);
         aVals.push_back(1.);
      }
   }

   // add Equation [1] to the solver
   auto aLhsMatrix = sMatrix(nEdges, nCorners, aRows, aCols, aVals);
   auto aRhsVector = new_array_ptr(aRhs);
   M->constraint("Intersection Agreement", 
         (Expr::mul(aLhsMatrix, a)),
         Domain::equalsTo(aRhsVector));
    return;
}
void IntrinsicFlattening::buildFaceConstraints(Model::t& M, Variable::t& a){
 // =========== Equation [2] =================================================
   
   // For each triangle, add a constraint that says that the
   // interior angles from the coherent angle system sum to π:
   //
   //    α0 + α1 + α2 = π       [2]
   //
   //        *
   //       /0\
   //      /   \
   //     /     \
   //    /1     2\
   //   *---------*
   //

   // temporary arrays for building sparse linear system for Equation [2]
   vector<int> triRows;
   vector<int> triCols;
   vector<double> triValues;
   vector<double> triRhs = vector<double>(nFaces, PI); // set all values to π

   for (Face f : mesh->faces()) {
      for (Corner c : f.adjacentCorners()) {
         triRows.emplace_back(f_[f]);
         triCols.emplace_back(c_[c]);
         triValues.emplace_back(1.); 
      }
   }

   // add Equation [2] to the solver
   auto triLhsMatrix = sMatrix(nFaces, nCorners, triRows, triCols, triValues);
   auto triRhsVector = new_array_ptr(triRhs);
   M->constraint("triangle sum constraint",
         Expr::mul(triLhsMatrix, a),
         Domain::equalsTo(triRhsVector));
    return;
}
void IntrinsicFlattening::buildVertexConstraints(Model::t& M, Variable::t& a){
   // =========== Equation [3] =================================================

   // For each interior vertex, add a constraint that says the
   // angles around the vertex sum to 2π:
   //
   //     _______
   //    /\     /\
   //   /  \   /  \
   //  /   2\1/0   \
   //  ------*------
   //  \   3/…\n   /
   //   \  /   \  /
   //    \/_____\/
   //
   //   

   // temporary vectors for building sparse matrices
   vector<int> vRows;
   vector<int> vCols;
   vector<double> vValues;
   vector<double> vRhs;


   // Sum to 2pi around each vertex
   vRhs = vector<double>(nVertices, 0);
   for (Vertex v : mesh->vertices()) {
      if(!v.isBoundary()) {
         vRhs[v_[v]] = 2*PI;
         for (Corner c : v.adjacentCorners()) {
            vRows.emplace_back(v_[v]);
            vCols.emplace_back(c_[c]);
            vValues.emplace_back(1.); 
         }
      } 
   }
   // add Equation [3] to the solver
   auto vLhsMatrix = sMatrix(nVertices, nCorners, vRows, vCols, vValues);
   auto vRhsVector = new_array_ptr(vRhs);
   M->constraint("vertex sum constraint", Expr::mul(vLhsMatrix, a), Domain::equalsTo(vRhsVector));
   return;
}
void IntrinsicFlattening::buildDelaunayConstraints(Model::t& M, Variable::t& a){
   // =========== Equation [4] =================================================
   // For each interior edge, add a constraint which ensures that
   // the two opposite angles in the CAS satisfy the Delaunay condition , namely
   //
   //   α0 + α1 < π       [4]
   //
   // where α are the angles in the coherent angle system and corners are indexed as below:
   //
   //        *
   //       /0\
   //      /   \
   //     /2   3\
   //    *-------*
   //     \5   4/
   //      \   /
   //       \1/
   //        *
   //


    // temporary vectors for building sparse matrices
    vector<int> dRows;
    vector<int> dCols;
    vector<double> dValues;
    vector<double> dRhs(nEdges, 1);

    // local delaunay constraint
    for (Edge e : mesh->edges()) {
        if(!e.isBoundary()) {
            size_t eInd = e_[e];
            dRhs[eInd] = PI;
            dRows.emplace_back(eInd);
            dCols.emplace_back(c_[e.halfedge().next().next().corner()]);
            dValues.emplace_back(1.); 

            dRows.emplace_back(eInd);
            dCols.emplace_back(c_[e.halfedge().twin().next().next().corner()]);
            dValues.emplace_back(1.); 
        } 
    }
    auto dLhsMatrix = sMatrix(nEdges, nCorners, dRows, dCols, dValues);
    auto dRhsVector = new_array_ptr(dRhs);
    M->constraint("Delaunay", Expr::mul(dLhsMatrix, a), Domain::lessThan(dRhsVector));
}


void IntrinsicFlattening::buildBoundaryObjective(Model::t& M, Variable::t& a, Variable::t& t, size_t excl, double interpolationWeight){
    vector<int> oRows;
    vector<int> oCols;
    vector<double> oValues;
    vector<double> oRhs;
    // Flat Boundary Constraint: Total geodesic curvature is 0 on "most" of the boundary
    size_t count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        count++;
    }
    double bdryCount = (double) count;
    count = 0;
    for (Vertex v : mesh->boundaryLoop(0).adjacentVertices()) {
        if (count >= excl) {
            double origSum = 0;
            for (Corner c : v.adjacentCorners()) {
                oRows.emplace_back(count);
                oCols.emplace_back(c_[c]);
                oValues.emplace_back(1.);
                origSum += geometry->cornerAngle(c);
            }
            oRhs.push_back((interpolationWeight * (PI -  2 * PI/bdryCount) + (1-interpolationWeight) * origSum));
        } else {
            oRhs.push_back(0);
        }
        count++;
    }
    auto oLhsMatrix = sMatrix(count, nCorners, oRows, oCols, oValues);
    auto oRhsTargetVector = new_array_ptr(oRhs);
    //M->constraint("Flat boundary", Expr::mul(bdryFlat, a), Domain::equalsTo(bdryPI));

    // just setup the objective: L^2 of curvature differences
    M->constraint(Expr::vstack(t, Expr::sub(Expr::mul(oLhsMatrix, a),oRhsTargetVector)), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);

}
void IntrinsicFlattening::buildOffsetConstraints(Model::t& M, Variable::t& alpha, Variable::t& beta){
    //TODO: comment this
    vector<int> offsetRows;
    vector<int> offsetCols;
    vector<double> offsetValues;
    // build LHS
    for (Corner c: mesh->corners()) {
        offsetCols.push_back(e_[c.halfedge().edge()]);
        offsetRows.push_back(c_[c]);
        offsetValues.push_back(1.);
        offsetCols.push_back(e_[c.halfedge().twin().next().edge()]);
        offsetRows.push_back(c_[c]);
        offsetValues.push_back(1.);
    }
    // build RHS
    vector<double> offsetRhs(nCorners);
    for (Corner c: mesh->corners()) {
        offsetRhs[c_[c]]=geometry->cornerAngle(c);
    }
    auto offsetLhsMatrix = sMatrix(nCorners, nEdges, offsetRows, offsetCols, offsetValues);
    auto offsetRhsVector = new_array_ptr(offsetRhs);
    M->constraint("Alpha-Beta constraint", Expr::sub(beta,Expr::mul(offsetLhsMatrix, alpha)), Domain::equalsTo(offsetRhsVector));
}
*/
CornerData<double> IntrinsicFlattening::solveIntrinsicOnly() {
    /*
    // Initialize MOSEK solver
    Model::t M = new Model();
    auto _M = finally([&]() { M->dispose(); });

    // Angle offsets per edge representing the deviation from the 
    // Euclidean angles intrinsically
    Variable::t alpha = M->variable("alpha", nEdges, Domain::inRange(-2 * PI, 2 * PI));

    // The intrisic "CAT angle" at each corner
    Variable::t Beta = M->variable("beta", nCorners, Domain::inRange(0, 2*PI));

    // the betas are given by offsets by alpha
    buildOffsetConstraints(M, alpha, Beta);

    // Vertex sums of betas to 2 pi for each interior vertex
    buildVertexConstraints(M, Beta);

    // Dummy variable for the objective value
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    vector<double> originalAngles(nCorners);
    for (Corner c: mesh->corners()) {
        originalAngles[c_[c]]=geometry->cornerAngle(c);
    }
    auto originalAnglesVector = new_array_ptr(originalAngles);

    // Minimize ||alpha||^2
    M->constraint(Expr::vstack(t, Expr::sub(Beta, originalAnglesVector)), Domain::inQCone());
    M->objective(ObjectiveSense::Minimize, t);
    M->solve();

    cout << M->getProblemStatus() << endl;
    cout << "Optimization done for intrinsic only flattening" << endl;
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;

    CornerData<double> beta(*mesh);
    //auto alphaVal = alpha->level();
    //for (int i = 0; i < nEdges; i++) {
    //    cout << (*alphaVal)[i] << endl;
    //}
    auto xVal = Beta->level();
    for (Corner c: mesh->corners()) {
        //cout << "setting index " <<  c_[c] << " to " << (*xVal)[c_[c]] << endl;
        beta[c] = (*xVal)[c_[c]];
        //beta[c] = PI;
    }
    for (Vertex v: mesh->vertices()) {
        if (!v.isBoundary()) {
            double accum = 0.;
            for (Corner c: v.adjacentCorners()) {
                accum += beta[c];
            }
            cout << "angle sum: " << accum << endl;
        }
    }
    //DEBUG'
    //for (Corner c: mesh->corners()) {
    //    cout << beta[c] << endl;
    //}

    //for (int i = 0; i < nCorners; i++) {
    //    //beta[i] = (*xVal)[i];
    //    cout << beta[i] << "wah" << endl;
    //}
    */

    //DEBUG
    CornerData<double> beta(*mesh);
    return beta;
}

// Given a CAT in the plane, and an assignment of new boundary curvatures (defined inline),
// returns intersection angles for a conformally equivalent CAT, and new CAT corner angles (which are the same as input)
SolutionData IntrinsicFlattening::solveFromPlane(double interpolationWeight) {

    // get corner angles of input CAT mesh
    CornerData<double> beta(*mesh);
    for (Corner c : mesh->corners()) {
        beta[c] = geometry->cornerAngle(c);
    }

    /*
    // Initialize MOSEK solver
    Model::t M = new Model();
    auto _M = finally([&]()
                      { M->dispose(); });

    // The values a at corners are the "angles" of a coherent angle system
    // compatible with the prescribed intersection angles (including new boundary angles)
    Variable::t a = M->variable("a", nCorners, Domain::inRange(0, 2 * PI));

    // intersection angles from a = intersection angles from beta
    buildIntersectionAngleConstraints(M, beta, a);
    // CAS sums to pi within each triangle
    buildFaceConstraints(M, a);
    // Vertex sums to 2 pi for each interior vertex
    buildVertexConstraints(M, a);
    // Delaunay constraint
    // Since the input mesh is assumed to be Delaunay this shouldn't actually be necessary?
    buildDelaunayConstraints(M, a);

    // Dummy variable for the objective value
    Variable::t t = M->variable("t", 1, Domain::unbounded());
    buildBoundaryObjective(M, a, t, 2, interpolationWeight);
    M->solve();
    cout << M->getProblemStatus() << endl;
    cout << "Optimization Done" << endl;
    std::cout << "Optimal primal objective: " << M->primalObjValue() << endl;

    EdgeData<double> thetaSolve(*mesh, 0);
    auto asize = a->getSize();
    auto asol = a->level();
    for (Edge e : mesh->edges()) {
        double a1 = e.halfedge().isInterior() && e.halfedge().isInterior()
                        ? (*asol)[c_[e.halfedge().next().next().corner()]]
                        : 0;
        double a2 =
            e.halfedge().twin().isInterior() && e.halfedge().twin().isInterior()
                ? (*asol)[c_[e.halfedge().twin().next().next().corner()]]
                : 0;
        thetaSolve[e] = PI - a1 - a2;
    }

    */
    FaceData<bool> fMask(*mesh, true);
    EdgeData<bool> eMask(*mesh, true);
    EdgeData<bool> eBdry(*mesh, false);
    for (Edge e : mesh->edges())
    {
        eBdry[e] = e.isBoundary();
    }
    // DEBUG: TEMPORARY
    EdgeData<double> thetaSolve(*mesh, 0);
    return {mesh->vertex(0), eMask, eBdry, fMask, beta, thetaSolve};
}
