#include "CatOpt.h"
#include "Solver.h"
#include "CirclePatterns.h"

void CatOpt::circlePatterns() {
    targetAngles = CornerData<double>(*mesh);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        if (h.isInterior()) {
            double angle = geometry->cornerAngle(C);
            
            //cout << "Inital:" << angle << " Next:";
            //angle += sol[eInd[h.edge()]];
            //angle += sol[eInd[h.next().next().edge()]];

            targetAngles[C] = angle;
            //cout << angle << endl;
            //cout << "angle is" << angle;
        }
    }
    /*
    for (Vertex v: mesh->vertices()) {
        double accum = 0;
        for (Corner C : v.adjacentCorners()) {
            accum += targetAngles[C];
        }
        if (v.isBoundary()) {
            cout << "boundary" << accum << endl;
        } else {
            cout << "not boundary" << accum << endl;
        }
    }
    */
    CirclePatterns prob(mesh, 0, sol, eInd, vInd, fInd, targetAngles);
    cout << "starting parameterization" << endl;
    prob.parameterize();
    cout << "parameterization done" << endl;
    prob.dbgSVG("wog.svg");
}