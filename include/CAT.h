#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <tuple>

#include "glm/glm.hpp"
#include "glm/gtx/norm.hpp"
#include "glm/gtx/string_cast.hpp"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
const double PI_2 = 1.57079632679489661923;
const double PI_4 = 0.785398163397448309616;

using glm::acos;
using glm::dot;
using glm::isnan;
using glm::vec2;
using glm::vec3;
using std::abs;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// Triangle Stuff
extern "C" {
    void triangulate(char *, struct triangulateio *, struct triangulateio *,
                     struct triangulateio *);
    void trifree(void *memptr);
}

double angle(vec2 v1, vec2 v2) {
    return acos(dot(normalize(v1), normalize(v2)));
}
vec2 prod(vec2 a, vec2 b) {
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}
vec2 div(vec2 a, vec2 b) {
    return vec2(((a.x * b.x + a.y * b.y) / (b.x * b.x + b.y * b.y)), ((a.y * b.x - a.x * b.y) / (b.x * b.x + b.y * b.y)));
}
double sc(double c) {
    return sqrt(c * c - 1.);
}

struct S2Triangle {
    vec3 i, j, k;
    double area;
};
struct H2Triangle {
    vec3 i, j, k;
    double area;
};

struct E2Triangle {
    vec2 i, j, k;
};


struct MobiusInfo {
    vec2 a, b, c, d;
};

vec2 ForwardMobius(vec2 z, MobiusInfo m) {
    return div((prod(m.a, z) + m.b), (prod(m.c, z) + m.d));
}

vec2 InvMobius(vec2 z, MobiusInfo m) {
    return div((prod(m.d, z) - m.b), (-prod(m.c, z) + m.a));
}

vec2 ComplexDet(vec2 a, vec2 b, vec2 c, vec2 d, vec2 e, vec2 f, vec2 g, vec2 h, vec2 i) {
    return prod(a, (prod(e, i) - prod(f, h))) - prod(b, (prod(d, i) - prod(f, g))) + prod(c, (prod(d, h) - prod(e, g)));
}

vec2 ForwardSphericalProj(vec3 wh) {
    return vec2(wh.x / (1. - wh.z), wh.y / (1. - wh.z));
}
vec3 InverseSphericalProj(vec2 z) {
    //GLM PORT
    return vec3((2.f * z) / (1.f + dot(z, z)), (-1.f + dot(z, z)) / (1.f + dot(z, z)));
}

double SphericalArea(vec3 A, vec3 B, vec3 C) {
    return 2. * atan((dot(A, cross(B, C))) / (1. + dot(A, B) + dot(B, C) + dot(A, C)));
}
vec3 SphericalBary(vec3 p, vec3 i, vec3 j, vec3 k) {
    double area = SphericalArea(i,j,k);
    double u = (SphericalArea(j, k, p) / area);
    double v = (SphericalArea(k, i, p) / area);
    double w = (SphericalArea(i, j, p) / area);
    return vec3(u, v, w);
}
vec2 ForwardHyperbolicProj(vec3 wh) {
    return vec2(wh.x / (1. + wh.z), wh.y / (1. + wh.z));
}
vec3 InverseHyperbolicProj(vec2 z) {
    return vec3(2.f * z / (1.f - dot(z, z)), (1.f + dot(z, z)) / (1.f - dot(z, z)));
}
double HDist(vec3 i, vec3 j) {
    return ((i.x * j.x + i.y * j.y) - i.z * j.z);
}

double HyperbolicArea(vec3 A, vec3 B, vec3 C) {
    double cosha = -HDist(B, C);
    double coshb = -HDist(A, C);
    double coshc = -HDist(B, A);
    double alpha = acos((-cosha + coshb * coshc) / (sc(coshb) * sc(coshc)));
    double beta = acos((-coshb + cosha * coshc) / (sc(cosha) * sc(coshc)));
    double gamma = acos((-coshc + coshb * cosha) / (sc(coshb) * sc(cosha)));
    return PI - alpha - beta - gamma;
}
vec3 HyperbolicBary(vec3 p, H2Triangle T) {
    double u = HyperbolicArea(T.j, T.k, p) / T.area;
    double v = HyperbolicArea(T.i, T.k, p) / T.area;
    double w = HyperbolicArea(T.i, T.j, p) / T.area;
    if (isnan(u) || isnan(v) || u + v + w > 1.01) {
        return vec3(1., 1., 1.);
    }
    return vec3(u, v, w);
}

vec3 BaryCalc(vec2 p, E2Triangle T) {
    double u, v, w;
    vec2 v0 = T.j - T.i, v1 = T.k - T.i, v2 = p - T.i;
    double d00 = dot(v0, v0);
    double d01 = dot(v0, v1);
    double d11 = dot(v1, v1);
    double d20 = dot(v2, v0);
    double d21 = dot(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.f - v - w;
    /*
    if (u < 0. || v < 0. || w < 0.)
    {
        // DEBUG
        cout << "Something bad happened" << endl;
        cout << u << endl;
        cout << v << endl;
        cout << w << endl;
        cout << to_string(T.i) << endl;
        cout << to_string(T.j) << endl;
        cout << to_string(T.k) << endl;
        cout << to_string(p) << endl;
        return vec3(1., 1., 1.);
    }
    */
    return vec3(u, v, w);
} 
// given a point x and non-collinear points i, j, k, projects x to the plane spanned by the three and computes barycentric coordinates
vec3 Bary3D (vec3 x, vec3 i, vec3 j, vec3 k) {
    vec3 u = j - i;
    vec3 v = k - i;
    vec3 n = cross( u, v );
    double oneOver4ASquared= 1.0 / dot( n, n );
    vec3 w = x - i;
    return vec3(1.0 - dot( cross( w, v ), n ) * oneOver4ASquared -dot( cross( u, w ), n ) * oneOver4ASquared ,
                dot( cross( w, v ), n ) * oneOver4ASquared,dot( cross( u, w ), n ) * oneOver4ASquared);
}
MobiusInfo SolveTrans(vec2 z1, vec2 z2, vec2 z3, vec2 w1, vec2 w2, vec2 w3) {
    vec2 z1w1 = prod(z1, w1);
    vec2 z2w2 = prod(z2, w2);
    vec2 z3w3 = prod(z3, w3);
    vec2 a = ComplexDet(z1w1, w1, vec2(1., 0), z2w2, w2, vec2(1, 0), z3w3, w3, vec2(1, 0));
    vec2 b = ComplexDet(z1w1, z1, w1, z2w2, z2, w2, z3w3, z3, w3);
    vec2 c = ComplexDet(z1, w1, vec2(1., 0), z2, w2, vec2(1., 0), z3, w3, vec2(1, 0));
    vec2 d = ComplexDet(z1w1, z1, vec2(1., 0), z2w2, z2, vec2(1., 0), z3w3, z3, vec2(1, 0));
    return MobiusInfo{a, b, c, d};
}
class CAT {
public:
    vec2 i, j, k;
    double a_ij, a_jk, a_ki;
    string to_string() {
        // header
        std::stringstream ss;
        ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
           << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
           << "<svg width=\"" + std::to_string(dims) + "\" height=\"" + std::to_string(dims) + "\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
        /*
               ss  << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl
               << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\"" << endl
               << " \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">" << endl 
               << "<svg width=\"" << 100 
               << "\" height=\"" << 100 << "\">" << endl
               << "<title>Collection of Geometry</title>" << endl;
               */
        // points
        ss << "<circle cx=\"" << i.x << "\" cy=\"" << i.y << "\" r=\"2\"/>" << endl;
        ss << "<circle cx=\"" << j.x << "\" cy=\"" << j.y << "\" r=\"2\"/>" << endl;
        ss << "<circle cx=\"" << k.x << "\" cy=\"" << k.y << "\" r=\"2\"/>" << endl;
        // arcs
        ss << to_arc(i, j, a_ij) << endl
           << to_arc(j, k, a_jk) << endl
           << to_arc(k, i, a_ki) << endl;
        ss << to_subdiv(i, j, a_ij, 10)
           << endl
           << to_subdiv(j, k, a_jk, 10) << endl
           << to_subdiv(k, i, a_ki, 10) << endl;

        // footer
        ss << "</svg>";

        return ss.str();
    }
    //default initializer to make stuff play nice
    CAT() {
    }
    //overload for just sidelengths
    CAT(double ij, double jk, double ki, double a1, double a2, double a3, bool normalize = false) : a_ij(a1), a_jk(a2), a_ki(a3) {
        i = vec2(0, 0);
        j = vec2(ij, 0);
        double angle = acos((-jk * jk + ki * ki + ij * ij) / (2 * ki * ij));
        k = vec2(ki * cos(angle), ki * sin(angle));
        std::vector<double> temp = {ij, jk, ki};
        double max_len = *std::max_element(temp.begin(), temp.end());
        dims = (max_len)*3;
        if (normalize) {
            i = vec2(dims / 4, dims / 4);
            j += i;
            k += i;
        }
    }
    // overload for specifying coordinates
    CAT(double i1, double i2, double j1, double j2, double k1, double k2, double a1, double a2, double a3, bool normalize = false) : a_ij(a1), a_jk(a2), a_ki(a3) {
        i = vec2(i1, i2);
        j = vec2(j1, j2);
        k = vec2(k1, k2);

        std::vector<double> temp = {i1, i2, j1, j2, k1, k2};
        double max_coord = *std::max_element(temp.begin(), temp.end());
        double min_coord = *std::min_element(temp.begin(), temp.end());
        dims = (max_coord - min_coord) * 3;
        if (normalize) {
            j -= i;
            k -= i;
            i = vec2(dims / 4, dims / 4);
            j += i;
            k += i;
        }
    }
    void to_svg(string filename) {
        std::ofstream f(filename, std::ofstream::out);
        f << to_string();
        f.close();
        return;
    }
    std::tuple<vector<double>, vector<int>> triangulation(int n) {
        vector<double> pointList;
        vector<double> pl1 = to_subdiv_vec(i, j, a_ij, n),
                       pl2 = to_subdiv_vec(j, k, a_jk, n),
                       pl3 = to_subdiv_vec(k, i, a_ki, n);
        pointList.insert(pointList.end(), pl1.begin(), pl1.end());
        pointList.insert(pointList.end(), pl2.begin(), pl2.end());
        pointList.insert(pointList.end(), pl3.begin(), pl3.end());
        vector<int> pointMarkerList(n, 1);
        vector<int> edgeList = {0};
        for (int i = 1; i < 3 * n; i++){
            edgeList.push_back(i);
            edgeList.push_back(i);
        }
        edgeList.push_back(0);
        triangulateio t = {
            pointList.data(),
            NULL,                   //pointattributelist;
            pointMarkerList.data(), //int *pointmarkerlist;
            3 * n,                  // number of points
            0,                      //numberofpointattributes;
            NULL,                   //int *trianglelist;                     
            NULL,                   //double *triangleattributelist;          
            NULL,                   //double *trianglearealist;               
            NULL,                   //int *neighborlist;                 
            0,                      //int numberoftriangles;                    
            0,                      //int numberofcorners;                  
            0,                      //int numberoftriangleattributes;          

            edgeList.data(), //int *segmentlist;           
            NULL,            //int *segmentmarkerlist;     
            3 * n,           //int numberofsegments;       

            NULL, //double *holelist;              
            0,    //int numberofholes;             

            NULL, //double *regionlist;            
            0,    //int numberofregions;           

            NULL, //int *edgelist;                 
            NULL, //int *edgemarkerlist;           
            NULL, //double *normlist;              
            0,    //int numberofedges;             
        };
        triangulateio out = {};
        triangulate((char *)"pzYqQ", &t, &out, NULL);
        vector<double> finalPoints(out.pointlist, out.pointlist + 2 * out.numberofpoints);
        vector<int> triangles(out.trianglelist, out.trianglelist + 3 * out.numberoftriangles);
        //trifree(&out);
        return std::make_pair(finalPoints, triangles);
    }
    string triangulation_svg(int n) {
        // header
        std::stringstream ss;
        ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
           << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
           << "<svg width=\"" + std::to_string(dims) + "\" height=\"" + std::to_string(dims) + "\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
        vector<double> finalPoints;
        vector<int> triangles;
        std::tie(finalPoints, triangles) = triangulation(n);
        for (int i = 0; i < triangles.size() / 3; i++) {
            ss << "<path d=\"M" << finalPoints[2 * triangles[3 * i]] << " " << finalPoints[2 * triangles[3 * i] + 1] << " "
               << "L" << finalPoints[2 * triangles[3 * i + 1]] << " " << finalPoints[2 * triangles[3 * i + 1] + 1] << " "
               << "L" << finalPoints[2 * triangles[3 * i + 2]] << " " << finalPoints[2 * triangles[3 * i + 2] + 1] << " Z"
               << "\" fill=\"none\" stroke=\"green\" stroke-width=\"1\" />";
        }

        // footer
        ss << "</svg>";

        return ss.str();
    }
    vec3 planeToBary(vec2 x) {
        double alpha = angle(j - i, k - i) + a_ij + a_ki;
        double beta = angle(i - j, k - j) + a_ij + a_jk;
        double gamma = angle(i - k, j - k) + a_ki + a_jk;
        if (a_ij + a_jk + a_ki == 0.) {
            // Euclidean case
            vec2 ti, tj, tk;
            ti = vec2(0., 0.);
            double t_i = alpha;
            double t_j = beta;
            tj = vec2(cos(t_i) + sin(t_i) / tan(t_j), 0.);
            tk = vec2(cos(t_i), sin(t_i));
            MobiusInfo M = SolveTrans(ti, tj, tk, i, j, k);
            return BaryCalc(InvMobius(x, M), E2Triangle{ti, tj, tk}); 
        } else if (a_ij + a_jk + a_ki > 0.) {
            // spherical case
            vec3 ti, tj, tk;
            ti = vec3(0., 0., -1.);

            // law of cosines stuff to find side lenghts
            // a = cosine of length of jk, b = same for ik, c = same for ij
            double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
            double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
            double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));

            double tj1 = sqrt(1. - c * c);
            tj = vec3(tj1, 0., -c);
            double tk1 = (a - b * c) / tj1;
            double tk2 = sqrt(1. - tk1 * tk1 - b * b);
            tk = vec3(tk1, tk2, -b);
            MobiusInfo M = SolveTrans(ForwardSphericalProj(ti), ForwardSphericalProj(tj), ForwardSphericalProj(tk), i, j, k);

            return SphericalBary(InverseSphericalProj(InvMobius(x,M)), ti, tj, tk);
        } else {
            // hyperbolic case
            vec3 ti, tj, tk;
            // law of cosines stuff to find side lenghts
            // a = cosh of length of jk, b = cosh of length of ik, c = cosh of length of ij
            double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
            double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
            double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
            double tj1 = sqrt(c * c - 1.);

            double tk1 = (-a + b * c) / tj1;
            double tk2 = sqrt(-1. - tk1 * tk1 + b * b);

            ti = vec3(0., 0., 1.);
            tj = vec3(tj1, 0., c);
            tk = vec3(tk1, tk2, b);
            MobiusInfo M = SolveTrans(ForwardHyperbolicProj(ti), ForwardHyperbolicProj(tj), ForwardHyperbolicProj(tk), i, j, k);
            vec3 x_hyp = InverseHyperbolicProj(InvMobius(x,M));
            vec3 n = cross(tj - ti, tk - ti);
            return Bary3D(x_hyp * dot(ti, n)/dot(x_hyp, n) , ti, tj, tk);
        }
    }

private:
    int dims;
    string to_arc(vec2 a, vec2 b, double angle) {
        double radius = glm::distance(a, b) / abs(2 * sin(angle));
        string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
        // sweep flag is 1 if going outward, 0 if going inward
        string sweepFlag = angle < 0 ? "0" : "1";
        std::stringstream ss;
        ss << "<path d=\"M" << a.x << "," << a.y << " A" << radius << ","
           << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
           << b.x << "," << b.y << "\" fill=\"red\" stroke=\"green\" stroke-width=\"1\" />";
        return ss.str();
    }
    vector<double> to_subdiv_vec(vec2 a, vec2 b, double angle, int n) {
        vector<double> pts;
        if (angle == 0) {
            for (int i = 0; i < n; i++) {
                pts.push_back((i/n) * a.x + ((n-i)/n) * b.x);
                pts.push_back((i/n) * a.y + ((n-i)/n) * b.y);
            }
            return pts;
        } else {
            double radius = glm::distance(a, b) / abs(2 * sin(angle));
            vec2 c = center(a, b, angle);
            double theta = atan2(a.y - c.y, a.x - c.x);
            for (int i = 0; i < n; i++) {
                pts.push_back(c.x + cos(theta + 2 * i * angle / n) * radius);
                pts.push_back(c.y + sin(theta + 2 * i * angle / n) * radius);
            }
            return pts;
        }
    }
    string to_subdiv(vec2 a, vec2 b, double angle, int n) {
        double radius = glm::distance(a, b) / abs(2 * sin(angle));
        std::stringstream ss;
        ss << "<path d=\"M" << a.x << " " << a.y << " ";
        vec2 c = center(a, b, angle);
        double theta = atan2(a.y - c.y, a.x - c.x);
        for (int i = 1; i < n; i++)
        {
            ss << "L" << c.x + cos(theta + 2 * i * angle / n) * radius << " "
               << c.y + sin(theta + 2 * i * angle / n) * radius << " ";
        }
        ss << "\" fill=\"none\" stroke=\"green\" stroke-width=\"1\" />";
        //ss << "<circle cx=\"" << c.x << "\" cy=\"" << c.y << "\" r=\"2\"/>";
        return ss.str();
    }
    vec2 center(vec2 a, vec2 b, double angle) {
        // |b-a| = 2 sin(alpha) r, so radius = |b-a| / 2 sin(alpha)
        double radius = glm::distance(a, b) / abs(2 * sin(angle));
        vec2 diff = b - a;
        // counterclockwise rotation by 90 degrees
        vec2 offset = vec2(-diff.y, diff.x);
        if (angle <= 0) {
            offset *= -cos(angle) * radius / glm::length(diff);
        } else {
            offset *= cos(angle) * radius / glm::length(diff);
        }
        return (a + b) / 2.f + offset;
    }
};


vec2 baryToPlane(vec2 i, vec2 j, vec2 k, double a_ij, double a_jk, double a_ki, float i_b, float j_b, float k_b) {
    double alpha = angle(j - i, k - i) + a_ij + a_ki;
    double beta = angle(i - j, k - j) + a_ij + a_jk;
    double gamma = angle(i - k, j - k) + a_ki + a_jk;
    // DEBUG
    //cout << endl << "expected area ratios: " << endl << i_b << "," << j_b << "," << k_b << endl;
    if (a_ij + a_jk + a_ki == 0.) {
        // Euclidean case
        vec2 ti, tj, tk;
        ti = vec2(0., 0.);
        double t_i = alpha;
        double t_j = beta;
        tj = vec2(cos(t_i) + sin(t_i) / tan(t_j), 0.);
        tk = vec2(cos(t_i), sin(t_i));
        // DEBUG
        vec2 p = i_b * ti + j_b * tj + k_b * tk;
        //cout << "actual area ratios (E2): " << endl <<
        string temp = glm::to_string(BaryCalc(p, E2Triangle{ti, tj, tk})); // << endl;
        //
        MobiusInfo M = SolveTrans(ti, tj, tk, i, j, k);
        return ForwardMobius(i_b * ti + j_b * tj + k_b * tk, M);
    } else if (a_ij + a_jk + a_ki > 0.) {
        // spherical case
        vec3 ti, tj, tk;
        ti = vec3(0., 0., -1.);

        // law of cosines stuff to find side lenghts
        // a = cosine of length of jk, b = same for ik, c = same for ij
        double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));

        /*
        cout << "hmm" << endl;
        cout << alpha << endl;
        cout << beta << endl;
        cout << gamma << endl;
        cout << c << endl;
        */
        double tj1 = sqrt(1. - c * c);
        tj = vec3(tj1, 0., -c);
        double tk1 = (a - b * c) / tj1;
        double tk2 = sqrt(1. - tk1 * tk1 - b * b);
        tk = vec3(tk1, tk2, -b);
        MobiusInfo M = SolveTrans(ForwardSphericalProj(ti), ForwardSphericalProj(tj), ForwardSphericalProj(tk), i, j, k);
        vec3 coord = i_b * ti + j_b * tj + k_b * tk;
        coord = coord / l2Norm(coord);

        /*
        cout << "vertex stuff i:" << glm::to_string((ti)) << endl;
        cout << "vertex stuff j:" << glm::to_string((tj)) << endl;
        cout << "vertex stuff k:" << glm::to_string((tk)) << endl;
        cout << "coords" << glm::to_string((coord)) << endl;

        */
        //cout << "vertex stuff p:" << glm::to_string((coord)) << endl;
        //

        // DEBUG
        //cout << "actual area ratios (S2): " << endl << glm::to_string(SphericalBary(coord, S2Triangle{ti, tj, tk, SphericalArea(ti,tj,tk)})) << endl;
        //

        return ForwardMobius(ForwardSphericalProj(coord), M);
    } else {
        // hyperbolic case
        vec3 ti, tj, tk;
        // law of cosines stuff to find side lenghts
        // a = cosh of length of jk, b = cosh of length of ik, c = cosh of length of ij
        double a = (cos(alpha) + cos(gamma) * cos(beta)) / (sin(gamma) * sin(beta));
        double b = (cos(beta) + cos(alpha) * cos(gamma)) / (sin(alpha) * sin(gamma));
        double c = (cos(gamma) + cos(alpha) * cos(beta)) / (sin(alpha) * sin(beta));
        double tj1 = sqrt(c * c - 1.);

        double tk1 = (-a + b * c) / tj1;
        double tk2 = sqrt(-1. - tk1 * tk1 + b * b);

        ti = vec3(0., 0., 1.);
        tj = vec3(tj1, 0., c);
        tk = vec3(tk1, tk2, b);
        MobiusInfo M = SolveTrans(ForwardHyperbolicProj(ti), ForwardHyperbolicProj(tj), ForwardHyperbolicProj(tk), i, j, k);
        vec3 coord = i_b * ti + j_b * tj + k_b * tk;
        coord = coord / (float)sqrt(-coord.x * coord.x - coord.y * coord.y + coord.z * coord.z);
        // DEBUG
        //cout << "actual area ratios (H2): " << endl << glm::to_string(HyperbolicBary(coord, H2Triangle{ti, tj, tk, HyperbolicArea(ti,tj,tk)})) << endl;

        return ForwardMobius(ForwardHyperbolicProj(coord), M);
    }
}
double l2DistSquared(double ij, double jk, double ki, double a_ij, double a_jk, double a_ki, float i1, float j1, float k1, float i2, float j2, float k2) {
    vec2 i = vec2(0, 0);
    vec2 j = vec2(ij, 0);
    double angle = acos((-jk * jk + ki * ki + ij * ij) / (2 * ki * ij));
    vec2 k = vec2(ki * cos(angle), ki * sin(angle));
    /*
    cout << "i:" << endl << glm::to_string(i) << endl;
    cout << "j:" << endl << glm::to_string(j) << endl;
    cout << "k:" << endl << glm::to_string(k) << endl;
    */
    vec2 p1 = baryToPlane(i, j, k, a_ij, a_jk, a_ki, i1, j1, k1);
    vec2 p2 = baryToPlane(i, j, k, a_ij, a_jk, a_ki, i2, j2, k2);
    //cout << "p1:" << endl << glm::to_string(p1) << endl;
    //cout << "p2:" << endl << glm::to_string(p2) << endl;
    return distance2(p1, p2);
}
