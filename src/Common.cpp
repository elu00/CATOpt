#include "Common.h"

double Angle (Vector2 u, Vector2 v) {
    return atan2(cross(u,v), dot(u,v));
}
// solves for angle offsets alpha given Euclidean angles ti and beta values
tuple<double, double, double> bendAngles(double t1, double t2, double t3, double b1, double b2, double b3) {
    double aij = (b1+b2-b3-t1-t2+t3)/2;
    double ajk = (-b1+b2+b3+t1-t2-t3)/2;
    double aki = (b1-b2+b3-t1+t2-t3)/2;
    return {aij, ajk, aki};
}
// isometrically projects triangle i j k to the plane.
// for convenience, j is assumed to be on the x-axis, and i is set to 0
tuple<Vector2,Vector2,Vector2> projectToPlane(Vector3 i, Vector3 j, Vector3 k) {
    j -= i; k -= i;
    Vector3 U = j.unit();
    Vector3 V = cross(j,cross(j,k)).unit();
    return {{0.,0.}, {j.norm(),0}, {dot(U, k), dot(V,k)}};
}

BezierTriangle Coefficients (Vector3 I, Vector3 J, Vector3 K, double Bi, double Bj, double Bk) {
    auto [i,j,k] = projectToPlane(I,J,K);
    Vector2 eij = i-j;
    Vector2 ejk = j-k;
    Vector2 eki = k-i;
    Vector2 nij = eij.rotate90();
    Vector2 njk = ejk.rotate90();
    Vector2 nki = eki.rotate90();

    double thetai = Angle(-nij,nki);
    double thetaj = Angle(-njk,nij);
    double thetak = Angle(-nki,njk);

    auto [aij,ajk,aki] = bendAngles(thetai,thetaj,thetak,Bi,Bj,Bk);

    Vector2 mij = (i+j)/2;
    Vector2 mjk = (j+k)/2;
    Vector2 mki = (k+i)/2;

    Vector2 p200 = i;
    Vector2 p020 = j;
    Vector2 p002 = k;
    double w200 = 1, w020 = 1, w002 = 1;

    Vector2 p110 = mij + tan(aij)*nij/2;
    double w110 = cos(aij);

    Vector2 p011 = mjk + tan(ajk)*njk/2;
    double w011 = cos(ajk);

    Vector2 p110 = mki + tan(aki)*nki/2;
    double w110 = cos(aki);
    
    return { p200, p020, p002, p110, p011, p101, w200, w020, w002, w110, w011, w101 };
}
Vector2 RationalBezierTriangle(BezierTriangle T, double t1, double t2, double t3) {

    double B200 = t1 * t1;
    double B020 = t2 * t2;
    double B002 = t3 * t3;
    double B110 = 2 * t1 * t2;
    double B011 = 2 * t2 * t3;
    double B101 = 2 * t1 * t3;
    Vector2 y = 
        B200 * T.w200 * T.p200 +
        B020 * T.w020 * T.p020 +
        B002 * T.w002 * T.p002 +
        B110 * T.w110 * T.p110 +
        B011 * T.w011 * T.p011 +
        B101 * T.w101 * T.p101;
    double h = 
        B200 * T.w200 +
        B020 * T.w020 +
        B002 * T.w002 +
        B110 * T.w110 +
        B011 * T.w011 +
        B101 * T.w101;
    return y/h;
}
