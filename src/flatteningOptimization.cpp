#include "CatOpt.h"



void CatOpt::dbgSVG(string filename) {
    std::ofstream ss(filename, std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
       << "<svg width=\"2000\" height=\"2000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
    for (Vertex v : mesh->vertices()) {
        auto thing = flattened[vInd[v]];
        ss << "<circle cx=\"" << shift(thing.x) << "\" cy=\"" << shift(thing.y) << "\" r=\"1\"/>" << endl;
    }
    for (Edge e : mesh->edges()) {
        Vector2 i = flattened[vInd[e.halfedge().vertex()]];
        Vector2 j = flattened[vInd[e.halfedge().twin().vertex()]];
        double angle = alphas[eInd[e]];
        // FROM NORMALIZATION
        double radius = 500 * norm(i - j) / abs(2 * sin(angle));
        if (abs(angle) > 1e-7) {
            string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
            // sweep flag is 1 if going outward, 0 if going inward
            string sweepFlag = angle < 0 ? "0" : "1";
            ss << "<path d=\"M" << shift(i.x) << "," << shift(i.y) << " A" << radius << ","
            << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
            << shift(j.x) << "," << shift(j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
        } else {
            ss << "<line x1=\"" << shift(i.x) << "\" x2=\"" << shift(j.x) 
            << "\" y1=\"" << shift(i.y) <<"\" y2=\"" << shift(j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
        }
        
    }
    // footer
    ss << "</svg>";
}

double CatOpt::confObjective(vector<Vector2> &flattened, vector<double>& alphas) {
    double residual = 0.;
    bool first = true;
    for (Vertex v: mesh->vertices()) {
        double accum = 0.;
        double targetAccum = 0;
        //double vAngle = geometry->vertexAngleSums[v];
        for (Corner C : v.adjacentCorners()) {
            targetAccum += targetAngles[C];
            Halfedge h = C.halfedge();
            Halfedge ab = h;
            Halfedge bc = h.next();
            Halfedge ca = h.next().next();
            if (h.isInterior()) {
                const Vector2 &a = flattened[vInd[ab.vertex()]];
                const Vector2 &b = flattened[vInd[bc.vertex()]];
                const Vector2 &c = flattened[vInd[ca.vertex()]];
                auto u = b - a;
                auto v = c - a;
                double angle = orientedAngle(u,v);
                // alphas are offsets from canonical halfedge orientation
                if (h.edge().halfedge() == h) {
                    angle += alphas[eInd[h.edge()]];
                } else {
                    angle -= alphas[eInd[h.edge()]];
                }
                if (ca.edge().halfedge() == ca) {
                    angle += alphas[eInd[ca.edge()]];
                } else {
                    angle -= alphas[eInd[ca.edge()]];
                }

                accum += angle;
            }
            
        }
        if (first) {
            if (v.isBoundary()) {
                //cout << "Actual boundary sum: " << accum << "target: " << targetAccum << endl;
            } 
            first = false;
        }
    }
    //for (Corner C: mesh->corners()) cout << targetAngles[C] << endl;
    double result = 0.;
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        if (h.isInterior()) {
            const Vector2 &a = flattened[vInd[ab.vertex()]];
            const Vector2 &b = flattened[vInd[bc.vertex()]];
            const Vector2 &c = flattened[vInd[ca.vertex()]];
            auto u = b - a;
            auto v = c - a;
            double angle = orientedAngle(u,v);
            // alphas are offsets from canonical halfedge orientation
            if (h.edge().halfedge() == h) {
                angle += alphas[eInd[h.edge()]];
            } else {
                angle -= alphas[eInd[h.edge()]];
            }
            if (ca.edge().halfedge() == ca) {
                angle += alphas[eInd[ca.edge()]];
            } else {
                angle -= alphas[eInd[ca.edge()]];
            }
            //cout << "target:" << targetAngles[C] << "   " << "actual:" << angle << endl;
            // check: sign issue?
            if (ab.vertex().isBoundary()) {
                //cout << "BOUNDARY " << angle << " target " << targetAngles[C] << endl;
            }
            result += sqr(angle - targetAngles[C]);
        }
    }
    return result/2;
}
void CatOpt::confGradient(vector<Vector2> &x, vector<double> &alphas, vector<Vector2> &flattened_grad, vector<double>& alphas_grad) {
    for (int i = 0; i < nVertices; i++) {
        flattened_grad[i] = Vector2::zero();
    }
    fill(alphas_grad.begin(), alphas_grad.end(), 0);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        Halfedge ab = h;
        Halfedge bc = h.next();
        Halfedge ca = h.next().next();
        if (h.isInterior()) {
            const Vector2 &a = flattened[vInd[C.vertex()]];
            const Vector2 &b = flattened[vInd[bc.vertex()]];
            const Vector2 &c = flattened[vInd[ca.vertex()]];
            Vector2 u = b - a;
            Vector2 v = c - a;
            // this rotate is ccw
            Vector2 bGrad = -u.rotate90()/u.norm2();
            Vector2 cGrad = v.rotate90()/v.norm2();


            double angle = orientedAngle(u,v);
            //cout << "angle is" << angle << "\n";
            // alphas are offsets from canonical halfedge orientation
            if (h.edge().halfedge() == h) {
                angle += alphas[eInd[h.edge()]];
            } else {
                angle -= alphas[eInd[h.edge()]];
            }
            // fix second conditional in the objective
            if (ca.edge().halfedge() == ca) {
                angle += alphas[eInd[ca.edge()]];
            } else {
                angle -= alphas[eInd[ca.edge()]];
            }
            double diff = angle - targetAngles[C];
            flattened_grad[vInd[bc.vertex()]] += diff * bGrad;
            flattened_grad[vInd[ca.vertex()]] += diff * cGrad;
            // TODO double check this
            flattened_grad[vInd[C.vertex()]] -= diff*(bGrad + cGrad);
            if (h.edge().halfedge() == h) {
                alphas_grad[eInd[h.edge()]] += diff;
            } else {
                alphas_grad[eInd[h.edge()]] -= diff;
            }
            if (ca.edge().halfedge() == ca) {
                alphas_grad[eInd[ca.edge()]] += diff;
            } else {
                alphas_grad[eInd[ca.edge()]] -= diff;
            }
        }
    }
}
void CatOpt::confStep(int n) {
    //for (auto& v : subdiv_points) v *= 1.1;
    cout << "Starting planar optimization" << endl;
    vector<Vector2> flattened_new = flattened;
    vector<Vector2> flattened_grad(nVertices);
    vector<double> alphas_new = alphas;
    vector<double> alphas_grad(nEdges);
    double grad_size = 1.;

    for (int m = 0; m < n; m++){
    //while (sqrt(grad_size) > ep) {
        double result = confObjective(flattened, alphas);
        confGradient(flattened, alphas, flattened_grad, alphas_grad);
        grad_size = 0; // change this
        for (Vector2 a: flattened_grad) grad_size += a.norm2();
        for (double a: alphas_grad) grad_size += a*a;
        
        //double t = .00001;
        double t = 0.01;
        int steps = 0;
        while (confObjective(flattened_new, alphas_new) > result - alpha * t * grad_size) {
            steps++;
            t = beta * t;
            for (int i = 0; i < nVertices; i++) {
                flattened_new[i] = flattened[i] - t * flattened_grad[i];
                //cout << "flat" << flattened_grad[i].x << " " << flattened_grad[i].y << endl;
            }
            for (int i = 0; i < nEdges; i++) {
                alphas_new[i] = alphas[i] - t * alphas_grad[i];
                //cout << "a" << alphas_grad[i] << endl;
            }
            /*
            int i = 6;
            flattened_new[i] = flattened[i] + Vector2{t,0};
            //alphas_new[i] = alphas[i] + t;
            double new_result = confObjective(flattened_new, alphas_new);
            cout << "Objective change:" << new_result - result << endl;
            cout <<  "Estimated change:" << t *(flattened_grad[i].x) << endl;
            */
        }
        
        flattened = flattened_new;
        alphas = alphas_new;
        if (iter % 1000 == 0) {
            cout << "Starting iteration " << iter << endl;
            cout << "t: " << t << endl;
            cout << "grad" << grad_size << endl;
            cout << "steps: " << steps << endl;
            //cout << "grad size squared:" << grad_size << endl;
            cout << "objective:" << result << endl;
            dbgSVG("step" + std::to_string(iter + 1) + ".svg");
            //std::ofstream s ("chug" + std::to_string(iter), std::ofstream::out);
        }
        /*
        if (iter == 9999) {
            cout << "final";
            dbgOutput("final");
        }
        */
        iter++;
    }
    //cout << "iteration count: " << iter << endl;
    //cout << "grad size:" << sqrt(grad_size) << endl;
    //cout << "objective:" << objective(x1, x2, x3) << endl;
    //vector<Vector3> pos = subdiv_points;
}

void CatOpt::loadModel(const std::string& inputPath, bff::Model& model,
			   std::vector<bool>& surfaceIsClosed) {
	std::string error;
	if (bff::MeshIO::read(inputPath, model, error)) {
        bff::Mesh& mesh = model[0];
        int nBoundaries = (int)mesh.boundaries.size();
        if (nBoundaries >= 1) {
            // mesh has boundaries
            int eulerPlusBoundaries = mesh.eulerCharacteristic() + nBoundaries;
            if (eulerPlusBoundaries == 2) {
                // fill holes if mesh has more than 1 boundary
                if (nBoundaries > 1) {
                    if (bff::HoleFiller::fill(mesh)) {
                        // all holes were filled
                        surfaceIsClosed[0] = true;
                    }
                }
            } else {
                // mesh probably has holes and handles
                bff::HoleFiller::fill(mesh, true);
                bff::Generators::compute(mesh);
            }

        } else if (nBoundaries == 0) {
            if (mesh.eulerCharacteristic() == 2) {
                // mesh is closed
                surfaceIsClosed[0] = true;
            } else {
                // mesh has handles
                bff::Generators::compute(mesh);
            }
        }

	} else {
		std::cerr << "Unable to load file: " << inputPath << ". " << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

void CatOpt::flatten(bff::Model& model, const std::vector<bool>& surfaceIsClosed,
			 int nCones, bool flattenToDisk, bool mapToSphere) {
    bff::Mesh& mesh = model[0];
    bff::BFF bff(mesh);

    if (nCones > 0) {
        std::vector<bff::VertexIter> cones;
        bff::DenseMatrix coneAngles(bff.data->iN);
        int S = std::min(nCones, (int)mesh.vertices.size() - bff.data->bN);

        if (bff::ConePlacement::findConesAndPrescribeAngles(S, cones, coneAngles, bff.data, mesh)
            == bff::ConePlacement::ErrorCode::ok) {
            if (!surfaceIsClosed[0] || cones.size() > 0) {
                bff::Cutter::cut(cones, mesh);
                bff.flattenWithCones(coneAngles, true);
            }
        }
    } else {
        if (surfaceIsClosed[0]) {
                std::cerr << "Surface is closed. Either specify nCones or mapToSphere." << std::endl;
                exit(EXIT_FAILURE);
        } else {
                bff::DenseMatrix u(bff.data->bN);
                bff.flatten(u, true);
        }
    }
}


void CatOpt::conformalFlatten() {
	// parse command line options
	std::string inputPath = inputMeshPath;
	int nCones = 0;
	bool flattenToDisk = false;
	bool mapToSphere = false;
	bool normalizeUVs = false;
	// load model
	std::vector<bool> surfaceIsClosed(1);
	loadModel(inputPath, model, surfaceIsClosed);

	// set nCones to 8 for closed surfaces`
    if (surfaceIsClosed[0] && !mapToSphere && nCones < 3) {
        std::cout << "Setting nCones to 8." << std::endl;
        nCones = 8;
    }

    // flatten
    flatten(model, surfaceIsClosed, nCones, flattenToDisk, mapToSphere);
    for (bff::VertexCIter v = model[0].vertices.begin(); v != model[0].vertices.end(); v++) {
        flattened.push_back({v->wedge()->uv.x, v->wedge()->uv.y});
        //cout << v->wedge()->uv.x << endl;
        //cout << v->wedge()->uv.y << endl;
    }
    // DEBUG
    alphas = vector<double>(nEdges, 0);
    targetAngles = CornerData<double>(*mesh);
    for (Corner C: mesh->corners()) {
        Halfedge h = C.halfedge();
        if (h.isInterior()) {
            double angle = geometry->cornerAngle(C);
            angle += sol[eInd[h.edge()]];
            angle += sol[eInd[h.next().next().edge()]];
            targetAngles[C] = angle;
            //cout << "angle is" << angle;
        }
    }
    /*
    // DEBUG
    for (Vertex v: mesh->vertices()) {
        double accum = geometry->vertexAngleSums[v];
        for (Edge e: v.adjacentEdges()) {
            accum += (e.isBoundary() ? 1 : 2)*sol[eInd[e]];
        }
        if (v.isBoundary()) {
            cout << "boundary from solve" << accum << endl;
        } else {
            cout << "not boundary from solve" << accum << endl;
        }
    }
    
    // DEBUG
    for (Vertex v: mesh->vertices()) {
        double accum = 0;
        //double vAngle = geometry->vertexAngleSums[v];
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
    dbgSVG("step0.svg");
    //dbgOutput("chug0");
    confStep(50000);  

}
// just for validating the SVG formula I'm using
void CatOpt::testSVG() {
    std::ofstream ss("test.svg", std::ofstream::out);
    ss << "<?xml version=\"1.0\" standalone=\"no\" ?>" << endl
       << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
       << "<svg width=\"2000\" height=\"2000\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" >" << endl;
    ss << "<circle cx=\"" << 500 << "\" cy=\"" << (500) << "\" r=\"1\"/>" << endl;
    ss << "<circle cx=\"" << 500 << "\" cy=\"" << (1000) << "\" r=\"1\"/>" << endl;

    Vector2 i = {500, 500};
    Vector2 j = {1000, 500};
    double angle = 3*PI/4;
    double radius = norm(i - j) / abs(2 * sin(angle));
    string largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
    // sweep flag is 1 if going outward, 0 if going inward
    string sweepFlag = angle < 0 ? "0" : "1";
    ss << "<path d=\"M" << (i.x) << "," << (i.y) << " A" << radius << ","
       << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
       << (j.x) << "," << (j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
    ss << "<line x1=\"" << (i.x) << "\" x2=\"" << (j.x)
       << "\" y1=\"" << (i.y) << "\" y2=\"" << (j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
    
    double wog = 1.3;
    j = {500 + 700*cos(wog), 500 + 800*sin(wog)};
    angle = -(PI - wog - angle);
    radius = norm(i - j) / abs(2 * sin(angle));
    largeArcFlag = std::abs(2 * angle) <= 3.14159265358979323846264 ? "0" : "1";
    // sweep flag is 1 if going outward, 0 if going inward
    sweepFlag = angle < 0 ? "0" : "1";
    ss << "<path d=\"M" << (i.x) << "," << (i.y) << " A" << radius << ","
       << radius << " 0 " << largeArcFlag << " " << sweepFlag << " "
       << (j.x) << "," << (j.y) << "\" fill=\"none\" stroke=\"green\" stroke-width=\"2\" />" << endl;
    ss << "<line x1=\"" << (i.x) << "\" x2=\"" << (j.x)
       << "\" y1=\"" << (i.y) << "\" y2=\"" << (j.y) << "\" stroke=\"blue\" stroke-width=\"2\"/>" << endl;
    ss << "</svg>";
}