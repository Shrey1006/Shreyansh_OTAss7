// grad.hpp - (36) MGDA and (37) GradNorm on the differentiable toy problem P13(a)
#pragma once
#include "classical.hpp"
namespace moo {
// min-norm point of the convex hull of the gradients: Frank-Wolfe on the simplex (closed form for M = 2)
inline Vec mgdaWeights(const std::vector<Vec>& G) {
    int M = (int)G.size();
    if (M == 2) { Vec d(G[0].size()); double nn = 0, nd = 0; for (size_t i = 0; i < d.size(); i++) { d[i] = G[0][i] - G[1][i]; nn += d[i]*d[i]; nd += (G[1][i]) * d[i]; }
        double a = nn < 1e-18 ? 0.5 : std::min(1.0, std::max(0.0, -nd / nn)); return {a, 1 - a}; }
    Vec a(M, 1.0 / M);
    for (int it = 0; it < 200; it++) {
        Vec d(G[0].size(), 0.0); for (int k = 0; k < M; k++) for (size_t i = 0; i < d.size(); i++) d[i] += a[k] * G[k][i];
        int best = 0; double bv = 1e300; for (int k = 0; k < M; k++) { double v = 0; for (size_t i = 0; i < d.size(); i++) v += G[k][i] * d[i]; if (v < bv) { bv = v; best = k; } }
        double gam = 2.0 / (it + 2); for (int k = 0; k < M; k++) a[k] *= (1 - gam); a[best] += gam; }
    return a;
}
struct GDTrace { std::vector<Vec> x, loss, w; };
// MGDA descent: x <- x - eta * sum alpha_k grad f_k, stops when the common direction vanishes (Pareto-stationary)
inline GDTrace run_mgda(const Problem& P, Vec x, int iters, double eta = 0.05) {
    GDTrace tr;
    for (int t = 0; t < iters; t++) {
        std::vector<Vec> G; for (int k = 0; k < P.M; k++) G.push_back(P.grad(x, k)); Vec a = mgdaWeights(G);
        Vec f; double cv; P.fn(x, f, cv); tr.x.push_back(x); tr.loss.push_back(f); tr.w.push_back(a);
        for (int i = 0; i < P.n; i++) { double d = 0; for (int k = 0; k < P.M; k++) d += a[k] * G[k][i]; x[i] -= eta * d; }
    }
    return tr;
}
// fixed equal weights baseline
inline GDTrace run_fixed(const Problem& P, Vec x, int iters, double eta = 0.05) {
    GDTrace tr;
    for (int t = 0; t < iters; t++) { Vec f; double cv; P.fn(x, f, cv); tr.x.push_back(x); tr.loss.push_back(f); tr.w.push_back(Vec(P.M, 1.0 / P.M));
        for (int i = 0; i < P.n; i++) { double d = 0; for (int k = 0; k < P.M; k++) d += P.grad(x, k)[i] / P.M; x[i] -= eta * d; } }
    return tr;
}
// GradNorm (Chen et al. 2018) : adapt task weights w_k so that G_k = w_k ||grad f_k|| matches mean(G) * r_k^alpha,
// r_k = relative inverse training rate (L_k(t)/L_k(0)) / mean.  Weight update by gradient step on L_grad = sum |G_k - target_k|.
inline GDTrace run_gradnorm(const Problem& P, Vec x, int iters, double eta = 0.05, double alpha = 1.5, double lrw = 0.05) {
    GDTrace tr; int M = P.M; Vec w(M, 1.0), L0;
    for (int t = 0; t < iters; t++) {
        Vec f; double cv; P.fn(x, f, cv); if (t == 0) L0 = f;
        std::vector<Vec> G; Vec gn(M), Gk(M); double Gm = 0;
        for (int k = 0; k < M; k++) { G.push_back(P.grad(x, k)); double s = 0; for (double v : G[k]) s += v*v; gn[k] = std::sqrt(s); Gk[k] = w[k] * gn[k]; Gm += Gk[k] / M; }
        Vec r(M); double rm = 0; for (int k = 0; k < M; k++) { r[k] = f[k] / std::max(L0[k], 1e-12); rm += r[k] / M; }
        for (int k = 0; k < M; k++) { double target = Gm * std::pow(r[k] / std::max(rm, 1e-12), alpha); double sgn = Gk[k] > target ? 1.0 : -1.0; w[k] = std::max(1e-3, w[k] - lrw * sgn * gn[k]); }
        double ws = 0; for (double v : w) ws += v; for (double& v : w) v *= M / ws;   // renormalise sum w = M
        tr.x.push_back(x); tr.loss.push_back(f); tr.w.push_back(w);
        for (int i = 0; i < P.n; i++) { double d = 0; for (int k = 0; k < M; k++) d += w[k] * G[k][i] / M; x[i] -= eta * d; }
    }
    return tr;
}
}  // namespace moo
