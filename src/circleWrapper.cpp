#include "CatOpt.h"
#include "Solver.h"
#include "CirclePatterns.h"

void CatOpt::circlePatterns() {
    targetAngles = CornerData<double>(*mesh);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        if (h.isInterior()) {
            double angle = geometry->cornerAngle(C);
            cout << "Inital:" << angle << " Next:";
            angle += sol[eInd[h.edge()]];
            angle += sol[eInd[h.next().next().edge()]];
            targetAngles[C] = angle;
            cout << angle << endl;
            //cout << "angle is" << angle;
        }
    }
    CirclePatterns prob(mesh, 0, sol, eInd, vInd, fInd, targetAngles);
    cout << "starting parameterization" << endl;
    prob.parameterize();
    cout << "parameterization done" << endl;
    prob.dbgSVG("wog.svg");
}