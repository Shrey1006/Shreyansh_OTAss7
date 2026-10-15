// classical.hpp - Algorithms 1-9 (Part 1 + 5E): Weighted Sum, eps-Constraint, Benson, Weight Metric,
// Value Function, Goal Programming, Lexicographic, TOPSIS, VIKOR.
// The scalar sub-problem is solved by our own derivative-free solver (1+1)-ES with 1/5th rule + restarts
// (allowed for C++: "write your own simple solver"). Sweeps, scalarisations, deviation variables and the
// priority logic are our own code.
#pragma once
#include "problems.hpp"
namespace moo {

inline Ind randInd(const Problem& P) {
    Vec x(P.n); for (int i = 0; i < P.n; i++) x[i] = urand(P.lo[i], P.hi[i]);
    return P.make(x);
}
inline double penalty(const Ind& a) { return g_deb ? 1e3 * a.cv : 0.0; }

inline Ind nelderMead(const Problem& P, const std::function<double(const Ind&)>& obj, const Ind& start, long evals);
// minimise obj(Ind) over the box. (1+1)-ES, success-rule step adaptation, 3 restarts.
inline Ind solveScalar(const Problem& P, const std::function<double(const Ind&)>& obj, long evals, const Vec* start = nullptr) {
    int restarts = 3; long per = std::max(10L, (evals * 3 / 4) / restarts);
    Ind best; double bestv = 1e300;
    for (int r = 0; r < restarts; r++) {
        Ind cur = (r == 0 && start) ? P.make(*start) : randInd(P);
        double cv = obj(cur) + penalty(cur), sigma = 0.2; int succ = 0, tries = 0;
        for (long t = 1; t < per; t++) {
            Vec x = cur.x; int forced = irand(P.n);
            for (int i = 0; i < P.n; i++)
                if (i == forced || urand() < std::min(1.0, 2.0 / P.n)) x[i] += sigma * (P.hi[i] - P.lo[i]) * nrand();
            Ind c = P.make(x); double v = obj(c) + penalty(c);
            tries++;
            if (v <= cv) { cur = c; cv = v; succ++; }
            if (tries == 10) { sigma *= (succ > 2) ? 1.5 : 0.65; sigma = std::min(0.5, std::max(1e-9, sigma)); succ = tries = 0; }
        }
        if (cv < bestv) { bestv = cv; best = cur; }
    }
    if (!P.integer && evals >= 200) { Ind pol = nelderMead(P, obj, best, evals / 4); if (obj(pol) + penalty(pol) < bestv) best = pol; }
    return best;
}


// Nelder-Mead simplex polish (own implementation) with restarts; handles narrow feasible edges better than the ES alone.
inline Ind nelderMead(const Problem& P, const std::function<double(const Ind&)>& obj, const Ind& start, long evals) {
    int n = P.n; long used = 0; Ind best = start; double bestv = obj(start) + penalty(start);
    double step = 0.05;
    while (used < evals) {
        std::vector<Ind> S(n + 1, best); std::vector<double> fv(n + 1, bestv);
        for (int i = 0; i < n; i++) { Vec x = best.x; double h = step * (P.hi[i] - P.lo[i]); x[i] += (x[i] + h <= P.hi[i]) ? h : -h; S[i + 1] = P.make(x); fv[i + 1] = obj(S[i + 1]) + penalty(S[i + 1]); used++; }
        auto ev = [&](const Vec& x, double& v) { Ind a = P.make(x); v = obj(a) + penalty(a); used++; return a; };
        for (int it = 0; it < 400 && used < evals; it++) {
            std::vector<int> o(n + 1); std::iota(o.begin(), o.end(), 0); std::sort(o.begin(), o.end(), [&](int a, int b) { return fv[a] < fv[b]; });
            std::vector<Ind> S2; std::vector<double> f2; for (int k : o) { S2.push_back(S[k]); f2.push_back(fv[k]); } S = S2; fv = f2;
            Vec c(n, 0.0); for (int k = 0; k < n; k++) for (int d = 0; d < n; d++) c[d] += S[k].x[d] / n;
            auto mix = [&](double t) { Vec x(n); for (int d = 0; d < n; d++) x[d] = c[d] + t * (S[n].x[d] - c[d]); return x; };
            double fr; Ind r = ev(mix(-1.0), fr);
            if (fr < fv[0]) { double fe; Ind e = ev(mix(-2.0), fe); if (fe < fr) { S[n] = e; fv[n] = fe; } else { S[n] = r; fv[n] = fr; } }
            else if (fr < fv[n - 1]) { S[n] = r; fv[n] = fr; }
            else { double fc; Ind cc = ev(mix(fr < fv[n] ? -0.5 : 0.5), fc);
                if (fc < std::min(fr, fv[n])) { S[n] = cc; fv[n] = fc; }
                else for (int k = 1; k <= n; k++) { Vec x(n); for (int d = 0; d < n; d++) x[d] = S[0].x[d] + 0.5 * (S[k].x[d] - S[0].x[d]); S[k] = ev(x, fv[k]); } }
        }
        int bi = (int)(std::min_element(fv.begin(), fv.end()) - fv.begin()); if (fv[bi] < bestv) { bestv = fv[bi]; best = S[bi]; }
        step *= 0.4;
    }
    return best;
}

inline std::vector<Vec> weightSet(int M, int K) {   // about K weight vectors on the simplex
    if (M == 2) { std::vector<Vec> w; for (int k = 0; k < K; k++) { double a = (double)k / (K - 1); w.push_back({a, 1 - a}); } return w; }
    int p = 1; while (true) { auto d = dasDennis(M, p); if ((int)d.size() >= K) return d; p++; }
}

// (1) Weighted Sum: min sum_i w_i * fn_i, vary w.
inline std::vector<Ind> run_weighted_sum(const Problem& P, long budget, int K = 20) {
    P.nev = 0; std::vector<Ind> out;
    for (auto& w : weightSet(P.M, K))
        out.push_back(solveScalar(P, [&](const Ind& a) { Vec v = P.norm(a.f); double s = 0; for (int i = 0; i < P.M; i++) s += w[i] * v[i]; return s; }, budget / K));
    return out;
}
// (2) eps-Constraint: min fn_0 s.t. fn_j <= eps_j (j>=1), eps swept over [0,1] in normalised space.
inline std::vector<Ind> run_eps_constraint(const Problem& P, long budget, int K = 20) {
    P.nev = 0; std::vector<Ind> out;
    for (int k = 0; k < K; k++) {
        Vec eps(P.M, 1.0);
        if (P.M == 2) eps[1] = (double)k / (K - 1); else for (int j = 1; j < P.M; j++) eps[j] = urand(0.3, 1.0);
        out.push_back(solveScalar(P, [&](const Ind& a) { Vec v = P.norm(a.f); double s = v[0];
            for (int j = 1; j < P.M; j++) s += 1e3 * std::max(0.0, v[j] - eps[j]); return s; }, budget / K));
    }
    return out;
}
// (3) Benson: from a feasible (usually dominated) point z, maximise the total improvement sum(z_i - f_i), f <= z.
inline std::vector<Ind> run_benson(const Problem& P, long budget, int K = 20) {
    P.nev = 0; std::vector<Ind> out;
    for (int k = 0; k < K; k++) {
        Ind z = randInd(P); int tries = 0; while (z.cv > 0 && tries++ < 200) z = randInd(P);
        Vec zn = P.norm(z.f);
        out.push_back(solveScalar(P, [&](const Ind& a) { Vec v = P.norm(a.f); double s = 0;
            for (int i = 0; i < P.M; i++) s += v[i] + 1e3 * std::max(0.0, v[i] - zn[i]); return s; }, std::max(30L, budget / K - 1), &z.x));
    }
    return out;
}
// (4) Weight metric: ( sum w_i |fn_i - z_i|^p )^(1/p), p = 1, 2 or infinity (p <= 0 means Tchebycheff). z = ideal point.
inline std::vector<Ind> run_weight_metric(const Problem& P, long budget, double p = 0, int K = 20) {
    P.nev = 0; std::vector<Ind> out;
    for (auto& w : weightSet(P.M, K))
        out.push_back(solveScalar(P, [&](const Ind& a) { Vec v = P.norm(a.f); double s = 0;
            for (int i = 0; i < P.M; i++) { double d = std::fabs(v[i]) ;
                if (p <= 0) s = std::max(s, w[i] * d); else s += w[i] * std::pow(d, p); }
            return p <= 0 ? s : std::pow(s, 1.0 / p); }, budget / K));
    return out;
}
// (5) Value function: DM utility U = prod (1.1 - fn_i)^{w_i} (Cobb-Douglas); maximise for several DM preference sets.
inline std::vector<Ind> run_value_function(const Problem& P, long budget, int K = 11) {
    P.nev = 0; std::vector<Ind> out;
    for (auto& w : weightSet(P.M, K))
        out.push_back(solveScalar(P, [&](const Ind& a) { Vec v = P.norm(a.f); double s = 0;
            for (int i = 0; i < P.M; i++) s -= w[i] * std::log(std::max(1.1 - v[i], 1e-9)); return s; }, budget / K));
    return out;
}
// (6) Goal programming: targets t_i (raw, minimisation form), one-sided deviation d+_i = max(0, f_i - t_i).
// min sum a_i d+_i.  a_i = 1 (unit weights) or 1/|t_i| (normalised weights).
inline Ind goal_programming(const Problem& P, const Vec& t, const Vec& a, long evals, Vec* dev = nullptr) {
    Ind r = solveScalar(P, [&](const Ind& s) { double z = 0; for (int i = 0; i < P.M; i++) z += a[i] * std::max(0.0, s.f[i] - t[i]); return z + 1e4 * s.cv; }, evals);
    if (dev) { dev->assign(P.M, 0.0); for (int i = 0; i < P.M; i++) (*dev)[i] = std::max(0.0, r.f[i] - t[i]); }
    return r;
}
// (7) Lexicographic: optimise priority 1, then priority 2 subject to f1 <= f1* + tol, ...
inline std::vector<Ind> lexicographic(const Problem& P, const std::vector<int>& order, double tolRel, long evals) {
    std::vector<Ind> sol; Vec bound(P.M, 1e300); long per = evals / (long)order.size();
    Ind cur;
    for (size_t s = 0; s < order.size(); s++) {
        int k = order[s];
        cur = solveScalar(P, [&](const Ind& a) { double z = a.f[k];
            for (size_t q = 0; q < s; q++) { int j = order[q]; z += 1e4 * std::max(0.0, a.f[j] - bound[j]); } return z + 1e4 * a.cv; }, per);
        bound[k] = cur.f[k] + tolRel * (P.fmax[k] - P.fmin[k]);
        sol.push_back(cur);
    }
    return sol;
}

// ---- MCDM ranking (not optimisers): decision matrix X[alt][crit], weights w, benefit[crit] ----
struct Ranking { std::vector<int> order; std::vector<double> score; };
inline Ranking topsis(const std::vector<Vec>& X, const Vec& w, const std::vector<bool>& benefit, bool verbose = false) {
    int m = (int)X.size(), n = (int)w.size();
    std::vector<Vec> V(m, Vec(n));
    for (int j = 0; j < n; j++) { double s = 0; for (int i = 0; i < m; i++) s += X[i][j]*X[i][j]; s = std::sqrt(s);
        if (verbose) std::cout << "  col " << j << " norm = " << s << "\n";
        for (int i = 0; i < m; i++) V[i][j] = w[j] * X[i][j] / s; }
    Vec best(n), worst(n);
    for (int j = 0; j < n; j++) { double mx = -1e300, mn = 1e300; for (int i = 0; i < m; i++) { mx = std::max(mx, V[i][j]); mn = std::min(mn, V[i][j]); }
        best[j] = benefit[j] ? mx : mn; worst[j] = benefit[j] ? mn : mx; }
    Ranking r; r.score.assign(m, 0.0);
    for (int i = 0; i < m; i++) { double dp = 0, dm = 0;
        for (int j = 0; j < n; j++) { dp += (V[i][j]-best[j])*(V[i][j]-best[j]); dm += (V[i][j]-worst[j])*(V[i][j]-worst[j]); }
        dp = std::sqrt(dp); dm = std::sqrt(dm); r.score[i] = dm / (dp + dm);
        if (verbose) { std::cout << "  alt " << i+1 << " weighted-normalised:"; for (double v : V[i]) std::cout << " " << v;
            std::cout << " | D+=" << dp << " D-=" << dm << " C=" << r.score[i] << "\n"; } }
    r.order.resize(m); std::iota(r.order.begin(), r.order.end(), 0);
    std::sort(r.order.begin(), r.order.end(), [&](int a, int b) { return r.score[a] > r.score[b]; });
    return r;
}
// VIKOR: S (group utility), R (individual regret), Q = v(S-S*)/(S^- - S*) + (1-v)(R-R*)/(R^- - R*); lower Q is better.
inline Ranking vikor(const std::vector<Vec>& X, const Vec& w, const std::vector<bool>& benefit, double v, bool verbose = false) {
    int m = (int)X.size(), n = (int)w.size();
    Vec fs(n), fm(n);
    for (int j = 0; j < n; j++) { double mx = -1e300, mn = 1e300; for (int i = 0; i < m; i++) { mx = std::max(mx, X[i][j]); mn = std::min(mn, X[i][j]); }
        fs[j] = benefit[j] ? mx : mn; fm[j] = benefit[j] ? mn : mx; }
    Vec S(m, 0.0), R(m, 0.0), Q(m, 0.0);
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) {
        double d = std::fabs(fs[j] - fm[j]) < 1e-14 ? 0 : w[j] * (fs[j] - X[i][j]) / (fs[j] - fm[j]);
        S[i] += d; R[i] = std::max(R[i], d); }
    double Ss = *std::min_element(S.begin(), S.end()), Sm = *std::max_element(S.begin(), S.end());
    double Rs = *std::min_element(R.begin(), R.end()), Rm = *std::max_element(R.begin(), R.end());
    for (int i = 0; i < m; i++) { Q[i] = v * (S[i]-Ss) / std::max(Sm-Ss, 1e-14) + (1-v) * (R[i]-Rs) / std::max(Rm-Rs, 1e-14);
        if (verbose) std::cout << "  alt " << i+1 << " S=" << S[i] << " R=" << R[i] << " Q=" << Q[i] << "\n"; }
    Ranking r; r.score = Q; r.order.resize(m); std::iota(r.order.begin(), r.order.end(), 0);
    std::sort(r.order.begin(), r.order.end(), [&](int a, int b) { return Q[a] < Q[b]; });
    return r;
}
}  // namespace moo
