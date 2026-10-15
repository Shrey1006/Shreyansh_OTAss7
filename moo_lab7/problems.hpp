// problems.hpp - problems P1..P14 (continuous ones) as objective + constraint-violation functions.
// All objectives are MINIMISED (maximisations are negated). cv is a normalised violation (0 = feasible).
// fmin/fmax are fixed per problem and are used to normalise objectives for ALL algorithms and metrics.
#pragma once
#include "moo_utils.hpp"
namespace moo {

struct Problem {
    std::string name;
    int n = 0, M = 2;
    Vec lo, hi;
    std::function<void(const Vec&, Vec&, double&)> fn;
    Vec fmin, fmax;          // normalisation bounds (raw space)
    Vec prefRef;             // decision-maker reference point (normalised space), used by R-NSGA-II
    bool integer = false;    // round decision variables (P2 integer variant)
    std::vector<Vec> trueFront;  // raw-space samples of the true front (may be empty)
    std::function<Vec(const Vec&, int)> grad;  // gradient of objective k (P13 only)
    mutable long nev = 0;    // evaluation counter

    Ind make(Vec x) const {
        for (int i = 0; i < n; i++) {
            x[i] = std::min(hi[i], std::max(lo[i], x[i]));
            if (integer) x[i] = std::round(x[i]);
        }
        Ind a; a.x = x; fn(a.x, a.f, a.cv); nev++;
        return a;
    }
    Vec norm(const Vec& f) const {
        Vec v(M);
        for (int i = 0; i < M; i++) v[i] = (f[i] - fmin[i]) / (fmax[i] - fmin[i]);
        return v;
    }
};

inline void fillTrue(Problem& p, const std::vector<Vec>& xs) {
    for (auto& x : xs) { Vec f; double cv; p.fn(x, f, cv); p.trueFront.push_back(f); }
}
inline Problem P1() {
    Problem p; p.name = "P1"; p.n = 1; p.M = 2; p.lo = {1}; p.hi = {10};
    p.fn = [](const Vec& x, Vec& f, double& cv) { f = {(x[0]-3)*(x[0]-3), (x[0]-7)*(x[0]-7)}; cv = 0; };
    p.fmin = {0, 0}; p.fmax = {16, 16};
    std::vector<Vec> xs; for (int k = 0; k < 200; k++) xs.push_back({3.0 + 4.0 * k / 199});
    fillTrue(p, xs); return p;
}
inline Problem P2(bool integer = false) {
    Problem p; p.name = integer ? "P2int" : "P2"; p.n = 1; p.M = 2; p.lo = {1}; p.hi = {10}; p.integer = integer;
    p.fn = [](const Vec& x, Vec& f, double& cv) { f = {100.0 / x[0], x[0]*x[0]}; cv = 0; };
    p.fmin = {10, 1}; p.fmax = {100, 100};
    std::vector<Vec> xs;
    if (integer) for (int k = 1; k <= 10; k++) xs.push_back({(double)k});
    else for (int k = 0; k < 200; k++) xs.push_back({1.0 + 9.0 * k / 199});
    fillTrue(p, xs); return p;
}
inline Problem P3() {  // ZDT2, n = 10
    Problem p; p.name = "P3"; p.n = 10; p.M = 2; p.lo = Vec(10, 0.0); p.hi = Vec(10, 1.0);
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        double s = 0; for (int i = 1; i < 10; i++) s += x[i];
        double g = 1 + s; f = {x[0], g * (1 - (x[0]/g)*(x[0]/g))}; cv = 0; };
    p.fmin = {0, 0}; p.fmax = {1, 1};
    for (int k = 0; k < 200; k++) { double t = k / 199.0; p.trueFront.push_back({t, 1 - t*t}); }
    return p;
}
inline Problem P4() {  // portfolio, repair by normalisation, f = (risk, -return)
    Problem p; p.name = "P4"; p.n = 4; p.M = 2; p.lo = Vec(4, 0.0); p.hi = Vec(4, 1.0);
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        static const double mu[4] = {0.08, 0.12, 0.15, 0.10};
        static const double S[4][4] = {{0.010,0.002,0.001,0.003},{0.002,0.040,0.006,0.004},{0.001,0.006,0.090,0.005},{0.003,0.004,0.005,0.020}};
        double s = 0; for (int i = 0; i < 4; i++) s += x[i];
        double w[4]; for (int i = 0; i < 4; i++) w[i] = s > 1e-12 ? x[i] / s : 0.25;
        double risk = 0, ret = 0;
        for (int i = 0; i < 4; i++) { ret += mu[i]*w[i]; for (int j = 0; j < 4; j++) risk += w[i]*S[i][j]*w[j]; }
        f = {risk, -ret}; cv = 0; };
    p.fmin = {0.007, -0.15}; p.fmax = {0.09, -0.08};
    return p;
}
inline Problem P5() {  // cantilever beam
    Problem p; p.name = "P5"; p.n = 2; p.M = 2; p.lo = {0.01, 0.02}; p.hi = {0.10, 0.20};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        double b = x[0], h = x[1];
        f = {7850.0 * 1.0 * b * h, 4.0 * 1000.0 / (200e9 * b * h*h*h)};
        cv = std::max(0.0, (2.4e-5 - b*h*h) / 2.4e-5); };
    p.fmin = {3.8, 2.5e-5}; p.fmax = {157.0, 0.0175};
    return p;
}
inline Problem P6() {  // Binh-Korn
    Problem p; p.name = "P6"; p.n = 2; p.M = 2; p.lo = {0, 0}; p.hi = {5, 3};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        f = {4*x[0]*x[0] + 4*x[1]*x[1], (x[0]-5)*(x[0]-5) + (x[1]-5)*(x[1]-5)};
        double c1 = std::max(0.0, ((x[0]-5)*(x[0]-5) + x[1]*x[1] - 25) / 25.0);
        double c2 = std::max(0.0, (7.7 - ((x[0]-8)*(x[0]-8) + (x[1]+3)*(x[1]+3))) / 7.7);
        cv = c1 + c2; };
    p.fmin = {0, 4}; p.fmax = {136, 50};
    std::vector<Vec> xs;
    for (int k = 0; k <= 100; k++) xs.push_back({3.0 * k / 100, 3.0 * k / 100});
    for (int k = 1; k <= 100; k++) xs.push_back({3.0 + 2.0 * k / 100, 3.0});
    fillTrue(p, xs); return p;
}
inline Problem P7() {  // water allocation, f = (-0.8x1, -1.2x2)
    Problem p; p.name = "P7"; p.n = 2; p.M = 2; p.lo = {10, 20}; p.hi = {90, 90};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        f = {-0.8 * x[0], -1.2 * x[1]}; cv = std::max(0.0, (x[0] + x[1] - 90) / 90.0); };
    p.fmin = {-56, -96}; p.fmax = {-8, -24};
    std::vector<Vec> xs; for (int k = 0; k <= 200; k++) { double a = 10 + 60.0 * k / 200; xs.push_back({a, 90 - a}); }
    fillTrue(p, xs); return p;
}
inline Problem P10(int M) {  // DTLZ2
    Problem p; p.name = "P10_M" + std::to_string(M); p.M = M; p.n = M + 9; p.lo = Vec(p.n, 0.0); p.hi = Vec(p.n, 1.0);
    p.fn = [M](const Vec& x, Vec& f, double& cv) {
        int n = (int)x.size(); double g = 0;
        for (int i = M - 1; i < n; i++) g += (x[i]-0.5)*(x[i]-0.5);
        f.assign(M, 0.0);
        const double h = M_PI / 2;
        for (int m = 1; m <= M; m++) {
            double v = 1 + g;
            for (int i = 1; i <= M - m; i++) v *= std::cos(x[i-1] * h);
            if (m > 1) v *= std::sin(x[M - m] * h);
            f[m-1] = v;
        }
        cv = 0; };
    p.fmin = Vec(M, 0.0); p.fmax = Vec(M, 1.0);
    int pdiv = (M == 3) ? 12 : (M == 5 ? 6 : 8);
    for (auto w : dasDennis(M, pdiv)) { double s = 0; for (double v : w) s += v*v; s = std::sqrt(s); for (double& v : w) v /= s; p.trueFront.push_back(w); }
    return p;
}
inline Problem P11() {  // production planning as a 2-objective problem (min -profit, min -units)
    Problem p; p.name = "P11"; p.n = 2; p.M = 2; p.lo = {0, 0}; p.hi = {400, 400};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        f = {-(4*x[0] + 6*x[1]), -(x[0] + x[1])};
        cv = std::max(0.0, (2*x[0] + 3*x[1] - 1200) / 1200.0) + std::max(0.0, (4*x[0] + 2*x[1] - 1600) / 1600.0); };
    p.fmin = {-2400, -500}; p.fmax = {0, 0};
    return p;
}
inline Problem P13() {  // toy multi-task problem, f_k = ||x - c_k||^2, with gradients
    Problem p; p.name = "P13"; p.n = 2; p.M = 2; p.lo = {-10, -10}; p.hi = {10, 10};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        f = {(x[0]-3)*(x[0]-3) + (x[1]-3)*(x[1]-3), (x[0]+3)*(x[0]+3) + (x[1]+3)*(x[1]+3)}; cv = 0; };
    p.grad = [](const Vec& x, int k) { double c = k == 0 ? 3.0 : -3.0; return Vec{2*(x[0]-c), 2*(x[1]-c)}; };
    p.fmin = {0, 0}; p.fmax = {72, 72};
    std::vector<Vec> xs; for (int k = 0; k <= 200; k++) { double t = -3 + 6.0 * k / 200; xs.push_back({t, t}); }
    fillTrue(p, xs); return p;
}
inline Problem P14() {  // laptop design, f = (weight, -runtime)
    Problem p; p.name = "P14"; p.n = 2; p.M = 2; p.lo = {13, 40}; p.hi = {17, 90};
    p.fn = [](const Vec& x, Vec& f, double& cv) {
        f = {0.4 + 0.05*x[0] + 0.012*x[1], -(0.2*x[1] / (0.5 + 0.03*x[0]))}; cv = 0; };
    p.fmin = {1.53, -20.2}; p.fmax = {2.13, -9.0};
    p.prefRef = {(1.8 - 1.53) / 0.6, (-12.0 + 20.2) / 11.2};   // reference point (1.8 kg, 12 h)
    std::vector<Vec> xs; for (int k = 0; k <= 200; k++) xs.push_back({13.0, 40 + 50.0 * k / 200});
    fillTrue(p, xs); return p;
}
inline Problem makeProblem(const std::string& key) {
    if (key == "P1") return P1();
    if (key == "P2") return P2(false);
    if (key == "P2int") return P2(true);
    if (key == "P3") return P3();
    if (key == "P4") return P4();
    if (key == "P5") { Problem p = P5(); p.prefRef = {0.1, 0.1}; return p; }
    if (key == "P6") return P6();
    if (key == "P7") return P7();
    if (key == "P10_3") return P10(3);
    if (key == "P10_5") return P10(5);
    if (key == "P11") return P11();
    if (key == "P13") return P13();
    if (key == "P14") return P14();
    throw std::runtime_error("unknown problem " + key);
}
}  // namespace moo
