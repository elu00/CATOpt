#include "Solver.h"
#include <deque>
#define beta 0.9
#define EPSILON 1e-9
#define MAX_ITERS 20000

SolverF::SolverF():
n(0)
{
    obj.reserve(MAX_ITERS);
}

void SolverF::gradientDescent()
{
    int k = 1;
    double f = 0.0, tp = 1.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    obj.clear(); obj.push_back(f);
    Eigen::VectorXd xp = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(n);

    while (true) {
        // compute momentum term
        v = x;
        if (k > 1) v += (k-2)*(x - xp)/(k+1);
        
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, v);
        
        // compute step size
        double t = tp;
        double fp = 0.0;
        Eigen::VectorXd xn = v - t*g;
        Eigen::VectorXd xnv = xn - v;
        handle->computeEnergy(fp, v);
        handle->computeEnergy(f, xn);
        while (f > fp + g.dot(xnv) + xnv.dot(xnv)/(2*t)) {
            t = beta*t;
            xn = v - t*g;
            xnv = xn - v;
            handle->computeEnergy(f, xn);
        }
    
        // update
        tp = t;
        xp = x;
        x  = xn;
        obj.push_back(f);
        k++;

        // check termination condition
        if (fabs(f - fp) < EPSILON || k > MAX_ITERS) break;
    }

    std::cout << "Iterations: " << k << std::endl;
}

void solvePositiveDefinite(Eigen::VectorXd& x,
                           const Eigen::VectorXd& b,
                           const Eigen::SparseMatrix<double>& A)
{
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(A);
    x = solver.solve(b);
}

void SolverF::newton()
{
    int k = 0;
    double f = 0.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    obj.clear(); obj.push_back(f);

    const double alpha = 0.5;
    while (true) {
        // compute update direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, x);

        Eigen::SparseMatrix<double> H(n, n);
        handle->computeHessian(H, x);

        Eigen::VectorXd p;
        solvePositiveDefinite(p, g, H);

        // compute step size
        double t = 1.0;
        double fp = f;
        handle->computeEnergy(f, x - t*p);
        while (f > fp - alpha*t*g.dot(p)) {
            t = beta*t;
            handle->computeEnergy(f, x - t*p);
        }

        // update
        x -= t*p;
        obj.push_back(f);
        k++;

        // check termination condition
        if (fabs(f - fp) < EPSILON || k > MAX_ITERS) break;
    }
    
    std::cout << "Iterations: " << k << std::endl;
}

bool solveTrustRegionQP(Eigen::VectorXd& x,
                        const Eigen::VectorXd& c,
                        const Eigen::SparseMatrix<double>& Q,
                        const double r)
{
    return true;
}

void SolverF::trustRegion()
{
    int k = 0;
    double f = 0.0;
    double r = 10.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    obj.clear(); obj.push_back(f);
    
    while (true) {
        // compute update length and direction
        Eigen::VectorXd g(n);
        handle->computeGradient(g, x);
        
        Eigen::SparseMatrix<double> H(n, n);
        handle->computeHessian(H, x);
        
        Eigen::VectorXd p;
        if (!solveTrustRegionQP(p, g, H, r)) {
            std::cout << "Unable to solve QP" << std::endl;
            break;
        }
        
        // compute rho
        double fp = f;
        handle->computeEnergy(f, x + p);
        double predicted = -p.dot(H*p) - g.dot(p);
        double actual = fp - f;
        double rho = actual/predicted;
        
        // update
        if (rho > 0) x += p;
        if (rho < 0.25) r *= 0.25;
        else if (rho > 0.75 && fabs(p.norm() - r) <= EPSILON) r = std::min(1000.0, 2*r);
        obj.push_back(f);
        k++;
        
        // check termination condition
        if (fabs(f - fp) < EPSILON || k > MAX_ITERS) break;
    }
    
    std::cout << "Iterations: " << k << std::endl;
}

void SolverF::lbfgs(int m)
{
    int k = 0;
    double f = 0.0;
    x = Eigen::VectorXd::Zero(n);
    handle->computeEnergy(f, x);
    obj.clear(); obj.push_back(f);
    Eigen::VectorXd g(n);
    handle->computeGradient(g, x);
    std::deque<Eigen::VectorXd> s;
    std::deque<Eigen::VectorXd> y;
    
    const double alpha = 1e-4;
    while (true) {
        // compute update direction
        int l = std::min(k, m);
        Eigen::VectorXd q = -g;
        
        Eigen::VectorXd a(l);
        for (int i = l-1; i >= 0; i--) {
            a(i) = s[i].dot(q) / y[i].dot(s[i]);
            q -= a(i)*y[i];
        }
        
        Eigen::VectorXd p = q;
        if (l > 0) p *= y[l-1].dot(s[l-1]) / y[l-1].dot(y[l-1]);
        
        for (int i = 0; i < l; i++) {
            double b = y[i].dot(p) / y[i].dot(s[i]);
            p += (a(i) - b)*s[i];
        }
        
        // compute step size
        double t = 1.0;
        double fp = f;
        handle->computeEnergy(f, x + t*p);
        while (f > fp + alpha*t*g.dot(p)) {
            t = beta*t;
            handle->computeEnergy(f, x + t*p);
        }
        
        // update
        Eigen::VectorXd xp = x;
        Eigen::VectorXd gp = g;
        x += t*p;
        handle->computeGradient(g, x);
        obj.push_back(f);
        k++;
        
        // update history
        if (k > m) {
            s.pop_front();
            y.pop_front();
        }
        s.push_back(x - xp);
        y.push_back(g - gp);
        
        // check termination condition
        if (fabs(f - fp) < EPSILON || k > MAX_ITERS) break;
    }
    
    std::cout << "Iterations: " << k << std::endl;
}
