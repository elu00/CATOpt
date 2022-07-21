#include "Common.h"

double Angle (Eigen::Vector2d u, Eigen::Vector2d v) {
    return atan2(u.x()*v.y() - u.y() * v.x(), u.x()*v.x() + u.y()*v.y());
}
BezierTriangle Coefficients () {
    return {};
}
Eigen::Vector2d RationalBezierTriangle(BezierTriangle T, double t1, double t2, double t3) {

    double B200 = t1 * t1;
    double B020 = t2 * t2;
    double B002 = t3 * t3;
    double B110 = 2 * t1 * t2;
    double B011 = 2 * t2 * t3;
    double B101 = 2 * t1 * t3;
    Eigen::Vector2d y = 
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
