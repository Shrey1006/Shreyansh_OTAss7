// aco.hpp - discrete problem P8 (bi-objective TSP) and the permutation-based algorithms:
// 20 MOEA/D-ACO, 25 MOACO, plus permutation versions of NSGA-II (14) and MOMA (30) for Task D3.
#pragma once
#include "evo.hpp"
namespace moo {

struct TSP {
    int n = 0; std::vector<Vec> xy, d, t;
    Vec fmin, fmax;   // normalisation (set from a reference front)
    static TSP make(int n = 15, unsigned seed = 7) {   // coordinates in [0,100]^2, u_ij ~ U(0.5,2) symmetric, mt19937 seed 7
        TSP s; s.n = n; std::mt19937 g(seed); std::uniform_real_distribution<double> U(0, 100), V(0.5, 2.0);
        for (int i = 0; i < n; i++) { double x = U(g), y = U(g); s.xy.push_back({x, y}); }
        s.d.assign(n, Vec(n, 0.0)); s.t.assign(n, Vec(n, 0.0));
        for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
            double dd = std::hypot(s.xy[i][0] - s.xy[j][0], s.xy[i][1] - s.xy[j][1]), u = V(g);
            s.d[i][j] = s.d[j][i] = dd; s.t[i][j] = s.t[j][i] = dd * u; }
        return s; }
    static TSP fromMatrices(std::vector<Vec> d, std::vector<Vec> t) { TSP s; s.n = (int)d.size(); s.d = d; s.t = t; return s; }
    Vec eval(const std::vector<int>& p) const { double a = 0, b = 0; for (int i = 0; i < n; i++) { int u = p[i], v = p[(i + 1) % n]; a += d[u][v]; b += t[u][v]; } return {a, b}; }
    Vec normf(const Vec& f) const { Vec v(2); for (int i = 0; i < 2; i++) v[i] = (f[i] - fmin[i]) / std::max(fmax[i] - fmin[i], 1e-12); return v; }
};
struct Tour { std::vector<int> p; Vec f; };

inline std::vector<int> nearestNeighbour(const TSP& s, const std::vector<Vec>& cost, int start = 0) {
    std::vector<int> p{start}; std::vector<bool> used(s.n, false); used[start] = true;
    for (int k = 1; k < s.n; k++) { int c = p.back(), b = -1; for (int j = 0; j < s.n; j++) if (!used[j] && (b < 0 || cost[c][j] < cost[c][b])) b = j; p.push_back(b); used[b] = true; }
    return p;
}
inline Tour mkTour(const TSP& s, std::vector<int> p, long& nev) { Tour t; t.p = p; t.f = s.eval(p); nev++; return t; }

// probabilistic construction of one tour; attract(i,j) = desirability of edge i->j
inline std::vector<int> constructTour(const TSP& s, const std::function<double(int, int)>& attract, bool verbose = false) {
    std::vector<int> p{irand(s.n)}; std::vector<bool> used(s.n, false); used[p[0]] = true;
    if (verbose) std::cout << "    start city " << p[0] << "\n";
    for (int k = 1; k < s.n; k++) {
        int c = p.back(); std::vector<double> w(s.n, 0.0); double tot = 0;
        for (int j = 0; j < s.n; j++) if (!used[j]) { w[j] = attract(c, j); tot += w[j]; }
        double r = urand() * tot, acc = 0; int pick = -1;
        for (int j = 0; j < s.n; j++) if (!used[j]) { acc += w[j]; pick = j; if (acc >= r) break; }
        if (verbose) { std::cout << "    from " << c << ": probs"; for (int j = 0; j < s.n; j++) if (!used[j]) std::cout << " [" << j << "]=" << w[j] / tot; std::cout << " r=" << r / tot << " -> " << pick << "\n"; }
        p.push_back(pick); used[pick] = true;
    }
    return p;
}
inline void addTourArchive(std::vector<Tour>& A, const Tour& c, int maxsz = 100) {
    for (auto& a : A) if (dominates(a.f, c.f) || a.f == c.f) return;
    A.erase(std::remove_if(A.begin(), A.end(), [&](const Tour& a) { return dominates(c.f, a.f); }), A.end());
    A.push_back(c);
    while ((int)A.size() > maxsz) { std::vector<Vec> F; for (auto& a : A) F.push_back(a.f); auto d = crowdingDistance(F);
        A.erase(A.begin() + (std::min_element(d.begin(), d.end()) - d.begin())); }
}
struct ACOParams { double alpha = 1.0, beta = 2.0, rho = 0.1; int nbr = 3; int N = 20; };

// (25) MOACO: two pheromone matrices, each ant has its own lambda, edge desirability (tau1^l tau2^(1-l))^alpha (eta1^l eta2^(1-l))^beta
inline std::vector<Tour> run_moaco(const TSP& s, long budget, ACOParams pr = ACOParams()) {
    int n = s.n, m = 20; long nev = 0; std::vector<Vec> t1(n, Vec(n, 1.0)), t2(n, Vec(n, 1.0)), e1(n, Vec(n, 0.0)), e2(n, Vec(n, 0.0));
    double md = 0, mt = 0; for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { md += s.d[i][j]; mt += s.t[i][j]; } md /= n * (n - 1); mt /= n * (n - 1);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { e1[i][j] = md / s.d[i][j]; e2[i][j] = mt / s.t[i][j]; }
    Vec nn = {s.eval(nearestNeighbour(s, s.d))[0], s.eval(nearestNeighbour(s, s.t))[1]};
    std::vector<Tour> arch;
    while (nev + m <= budget) {
        std::vector<Tour> iter;
        for (int k = 0; k < m; k++) { double l = (double)k / (m - 1);
            auto att = [&](int i, int j) { return std::pow(std::pow(t1[i][j], l) * std::pow(t2[i][j], 1 - l), pr.alpha) * std::pow(std::pow(e1[i][j], l) * std::pow(e2[i][j], 1 - l), pr.beta); };
            iter.push_back(mkTour(s, constructTour(s, att), nev)); }
        for (auto& c : iter) addTourArchive(arch, c);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) { t1[i][j] *= (1 - pr.rho); t2[i][j] *= (1 - pr.rho); }
        for (auto& a : arch) for (int i = 0; i < n; i++) { int u = a.p[i], v = a.p[(i + 1) % n];
            double d1 = 0.5 * nn[0] / a.f[0], d2 = 0.5 * nn[1] / a.f[1]; t1[u][v] += d1; t1[v][u] += d1; t2[u][v] += d2; t2[v][u] += d2; }
    }
    return arch;
}
// (20) MOEA/D-ACO: N scalar subproblems (Tchebycheff, lambda_i); each has its own pheromone matrix and one ant per iteration.
// A new tour updates the current solution of every neighbouring subproblem it improves; pheromone of subproblem i is
// reinforced by the current solutions of its nbr neighbours.  verbose prints every step for the hand trace.
inline std::vector<Tour> run_moeadaco(const TSP& s, long budget, ACOParams pr = ACOParams(), bool verbose = false) {
    int n = s.n, N = pr.N; long nev = 0; std::vector<Vec> W(N);
    for (int i = 0; i < N; i++) { double a = N == 1 ? 0.5 : (double)i / (N - 1); W[i] = {std::max(a, 1e-3), std::max(1 - a, 1e-3)}; }
    std::vector<std::vector<Vec>> tau(N, std::vector<Vec>(n, Vec(n, 1.0)));
    std::vector<Vec> e1(n, Vec(n, 0.0)), e2(n, Vec(n, 0.0)); double md = 0, mt = 0;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { md += s.d[i][j]; mt += s.t[i][j]; } md /= n * (n - 1); mt /= n * (n - 1);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { e1[i][j] = md / s.d[i][j]; e2[i][j] = mt / s.t[i][j]; }
    std::vector<std::vector<int>> B(N);
    for (int i = 0; i < N; i++) { std::vector<int> idx(N); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return std::fabs(W[a][0] - W[i][0]) < std::fabs(W[b][0] - W[i][0]); }); B[i].assign(idx.begin(), idx.begin() + std::min(pr.nbr, N)); }
    Vec z = {1e18, 1e18};
    auto g = [&](const Vec& f, const Vec& w) { Vec v = s.normf(f); return std::max(w[0] * v[0], w[1] * v[1]); };   // ideal = 0 in normalised space
    std::vector<Tour> cur(N), EP;
    for (int i = 0; i < N; i++) { std::vector<int> p(n); std::iota(p.begin(), p.end(), 0); std::shuffle(p.begin(), p.end(), rng); cur[i] = mkTour(s, p, nev); addTourArchive(EP, cur[i]); }
    int it = 0;
    while (nev + N <= budget) {
        if (verbose) std::cout << " iteration " << ++it << "\n";
        for (int i = 0; i < N; i++) {
            auto att = [&](int a, int b) { return std::pow(tau[i][a][b], pr.alpha) * std::pow(std::pow(e1[a][b], W[i][0]) * std::pow(e2[a][b], W[i][1]), pr.beta); };
            if (verbose) std::cout << "  subproblem " << i << " lambda=(" << W[i][0] << "," << W[i][1] << ")\n";
            Tour c = mkTour(s, constructTour(s, att, verbose), nev); addTourArchive(EP, c);
            if (verbose) std::cout << "   tour f=(" << c.f[0] << "," << c.f[1] << ")\n";
            for (int j : B[i]) if (g(c.f, W[j]) < g(cur[j].f, W[j])) { cur[j] = c; if (verbose) std::cout << "   replaces current solution of subproblem " << j << "\n"; }
        }
        for (int i = 0; i < N; i++) {
            for (int a = 0; a < n; a++) for (int b = 0; b < n; b++) tau[i][a][b] *= (1 - pr.rho);
            for (int j : B[i]) { double dep = 0.5 / (1 + 5 * g(cur[j].f, W[i]));
                for (int k = 0; k < n; k++) { int u = cur[j].p[k], v = cur[j].p[(k + 1) % n]; tau[i][u][v] += dep; tau[i][v][u] += dep; } }
            if (verbose) { std::cout << "  tau[" << i << "] row0:"; for (int b = 0; b < n; b++) std::cout << " " << tau[i][0][b]; std::cout << "\n"; }
        }
    }
    return EP;
}

// ---- permutation NSGA-II (14p) and MOMA (30p) ----
inline std::vector<int> oxCross(const std::vector<int>& a, const std::vector<int>& b) {
    int n = (int)a.size(), i = irand(n), j = irand(n); if (i > j) std::swap(i, j);
    std::vector<int> c(n, -1); std::vector<bool> used(n, false);
    for (int k = i; k <= j; k++) { c[k] = a[k]; used[a[k]] = true; }
    int pos = (j + 1) % n; for (int k = 0; k < n; k++) { int v = b[(j + 1 + k) % n]; if (!used[v]) { c[pos] = v; used[v] = true; pos = (pos + 1) % n; } }
    return c;
}
inline void permMutate(std::vector<int>& p) {
    int n = (int)p.size(), i = irand(n), j = irand(n); if (i > j) std::swap(i, j);
    if (urand() < 0.5) std::reverse(p.begin() + i, p.begin() + j + 1); else std::swap(p[i], p[j]);
}
inline std::vector<Tour> run_nsga2_perm(const TSP& s, long budget, bool memetic = false) {
    int N = 100, n = s.n; long nev = 0; std::vector<Tour> pop;
    for (int i = 0; i < N; i++) { std::vector<int> p(n); std::iota(p.begin(), p.end(), 0); std::shuffle(p.begin(), p.end(), rng); pop.push_back(mkTour(s, p, nev)); }
    std::vector<int> rk(N); std::vector<double> cr(N);
    auto assign = [&](std::vector<Tour>& P) { std::vector<Vec> F; for (auto& t : P) F.push_back(t.f); auto fr = fastNonDominatedSortF(F); rk.assign(P.size(), 0); cr.assign(P.size(), 0);
        for (size_t r = 0; r < fr.size(); r++) { std::vector<Vec> G; for (int i : fr[r]) G.push_back(F[i]); auto d = crowdingDistance(G); for (size_t k = 0; k < fr[r].size(); k++) { rk[fr[r][k]] = (int)r; cr[fr[r][k]] = d[k]; } } };
    assign(pop);
    const int nLS = 5, lsSteps = 8;
    while (nev + N + (memetic ? nLS * lsSteps : 0) <= budget) {
        std::vector<Tour> off;
        auto sel = [&] { int a = irand(N), b = irand(N); if (rk[a] != rk[b]) return rk[a] < rk[b] ? a : b; return cr[a] >= cr[b] ? a : b; };
        while ((int)off.size() < N) { auto c = urand() < 0.9 ? oxCross(pop[sel()].p, pop[sel()].p) : pop[sel()].p; if (urand() < 0.3) permMutate(c); off.push_back(mkTour(s, c, nev)); }
        if (memetic) for (int l = 0; l < nLS; l++) {   // 2-opt style local search, random weighted objective
            Tour cu = off[irand(N)]; double w = urand(); auto sc = [&](const Vec& f) { Vec v = s.normf(f); return w * v[0] + (1 - w) * v[1]; }; double cs = sc(cu.f);
            for (int t = 0; t < lsSteps; t++) { auto p = cu.p; int i = irand(n), j = irand(n); if (i > j) std::swap(i, j); std::reverse(p.begin() + i, p.begin() + j + 1);
                Tour c = mkTour(s, p, nev); if (sc(c.f) < cs) { cu = c; cs = sc(c.f); } }
            off.push_back(cu); }
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); std::vector<Vec> F; for (auto& t : R) F.push_back(t.f);
        auto fr = fastNonDominatedSortF(F); std::vector<Tour> nxt;
        for (auto& f : fr) { if ((int)(nxt.size() + f.size()) <= N) { for (int i : f) nxt.push_back(R[i]); continue; }
            std::vector<Vec> G; for (int i : f) G.push_back(F[i]); auto d = crowdingDistance(G); std::vector<int> idx(f.size()); std::iota(idx.begin(), idx.end(), 0);
            std::sort(idx.begin(), idx.end(), [&](int a, int b) { return d[a] > d[b]; }); for (size_t k = 0; (int)nxt.size() < N; k++) nxt.push_back(R[f[idx[k]]]); break; }
        pop = nxt; assign(pop);
    }
    std::vector<Tour> out; std::vector<Vec> F; for (auto& t : pop) F.push_back(t.f); auto fr = fastNonDominatedSortF(F); for (int i : fr[0]) out.push_back(pop[i]); return out;
}
// reference front for P8 (union of long runs); also sets the normalisation bounds of the TSP
inline std::vector<Vec> tspReferenceFront(TSP& s) {
    s.fmin = {0, 0}; s.fmax = {s.eval(nearestNeighbour(s, s.d))[0] * 2, s.eval(nearestNeighbour(s, s.t))[1] * 2};   // provisional scale for the local search
    std::vector<Vec> all; for (int r = 0; r < 3; r++) { seed(1000 + r); for (auto& t : run_nsga2_perm(s, 100000, r == 2)) all.push_back(t.f); }
    auto nd = nondominatedSubset(all); s.fmin = {1e18, 1e18}; s.fmax = {-1e18, -1e18};
    for (auto& f : nd) for (int i = 0; i < 2; i++) { s.fmin[i] = std::min(s.fmin[i], f[i]); s.fmax[i] = std::max(s.fmax[i], f[i]); }
    return nd;
}
}  // namespace moo
