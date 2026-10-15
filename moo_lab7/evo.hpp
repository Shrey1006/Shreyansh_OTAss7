// evo.hpp - evolutionary / decomposition / indicator algorithms:
// 10 VEGA, 11 MOGA, 12 NSGA, 13 NPGA, 14 NSGA-II, 15 NSGA-III, 16 SPEA2, 17 PESA-II, 18 MOEA/D, 19 MOEA/D-DE,
// 21 Constrained NSGA-II (Deb rules), 30 MOMA, 31 RVEA, 32 IBEA, 33 SMS-EMOA, 34 HypE, 35 R-NSGA-II.
#pragma once
#include "classical.hpp"
namespace moo {

inline int popSize(const Problem& P) { return P.M <= 3 ? 100 : 126; }
// objective vector used for selection: normalised, infeasible solutions pushed beyond the feasible box
inline Vec pn(const Problem& P, const Ind& a) {
    Vec v = P.norm(a.f);
    if (g_deb && a.cv > 0) for (double& t : v) t += 10.0 * a.cv + 1.0;
    return v;
}
inline double dist(const Vec& a, const Vec& b) { double s = 0; for (size_t i = 0; i < a.size(); i++) s += (a[i]-b[i])*(a[i]-b[i]); return std::sqrt(s); }

inline void sbx(Vec& c1, Vec& c2, const Problem& P, double eta = 15) {
    for (int i = 0; i < P.n; i++) {
        if (urand() > 0.5 || std::fabs(c1[i] - c2[i]) < 1e-14) continue;
        double y1 = std::min(c1[i], c2[i]), y2 = std::max(c1[i], c2[i]), lo = P.lo[i], hi = P.hi[i], r = urand();
        double beta = 1 + 2 * (y1 - lo) / (y2 - y1), alpha = 2 - std::pow(beta, -(eta + 1)), bq;
        bq = r <= 1 / alpha ? std::pow(r * alpha, 1 / (eta + 1)) : std::pow(1 / (2 - r * alpha), 1 / (eta + 1));
        double a = 0.5 * ((y1 + y2) - bq * (y2 - y1));
        beta = 1 + 2 * (hi - y2) / (y2 - y1); alpha = 2 - std::pow(beta, -(eta + 1));
        bq = r <= 1 / alpha ? std::pow(r * alpha, 1 / (eta + 1)) : std::pow(1 / (2 - r * alpha), 1 / (eta + 1));
        double b = 0.5 * ((y1 + y2) + bq * (y2 - y1));
        a = std::min(hi, std::max(lo, a)); b = std::min(hi, std::max(lo, b));
        if (urand() < 0.5) { c1[i] = a; c2[i] = b; } else { c1[i] = b; c2[i] = a; }
    }
}
inline void polyMut(Vec& x, const Problem& P, double eta = 20) {
    double pm = 1.0 / P.n;
    for (int i = 0; i < P.n; i++) {
        if (urand() > pm) continue;
        double lo = P.lo[i], hi = P.hi[i], d1 = (x[i] - lo) / (hi - lo), d2 = (hi - x[i]) / (hi - lo), r = urand(), dq;
        if (r < 0.5) dq = std::pow(2 * r + (1 - 2 * r) * std::pow(1 - d1, eta + 1), 1 / (eta + 1)) - 1;
        else dq = 1 - std::pow(2 * (1 - r) + 2 * (r - 0.5) * std::pow(1 - d2, eta + 1), 1 / (eta + 1));
        x[i] = std::min(hi, std::max(lo, x[i] + dq * (hi - lo)));
    }
}
template <class Sel>
std::vector<Ind> makeOffspring(const Problem& P, const std::vector<Ind>& pop, int n, Sel sel, double pc = 0.9) {
    std::vector<Ind> off;
    while ((int)off.size() < n) {
        Vec c1 = pop[sel()].x, c2 = pop[sel()].x;
        if (urand() < pc) sbx(c1, c2, P);
        polyMut(c1, P); polyMut(c2, P);
        off.push_back(P.make(c1)); if ((int)off.size() < n) off.push_back(P.make(c2));
    }
    return off;
}
inline std::vector<Ind> initPop(const Problem& P, int N) { std::vector<Ind> p; for (int i = 0; i < N; i++) p.push_back(randInd(P)); return p; }
inline std::vector<Ind> rankFront0(const std::vector<Ind>& pop) {
    std::vector<Ind> out; if (pop.empty()) return out;
    auto fr = fastNonDominatedSort(pop); for (int i : fr[0]) out.push_back(pop[i]); return out;
}

// ---- NSGA-II core ----
inline void assignRankCrowd(std::vector<Ind>& pop) {
    auto fr = fastNonDominatedSort(pop);
    for (size_t r = 0; r < fr.size(); r++) {
        std::vector<Vec> F; for (int i : fr[r]) F.push_back(pop[i].f);
        auto d = crowdingDistance(F);
        for (size_t k = 0; k < fr[r].size(); k++) { pop[fr[r][k]].rank = (int)r; pop[fr[r][k]].crowd = d[k]; }
    }
}
inline int tourCrowded(const std::vector<Ind>& pop) {
    int a = irand((int)pop.size()), b = irand((int)pop.size());
    if (pop[a].rank != pop[b].rank) return pop[a].rank < pop[b].rank ? a : b;
    if (pop[a].crowd != pop[b].crowd) return pop[a].crowd > pop[b].crowd ? a : b;
    return urand() < 0.5 ? a : b;
}
inline std::vector<Ind> nsga2Select(const std::vector<Ind>& R, int N) {
    auto fr = fastNonDominatedSort(R); std::vector<Ind> nxt;
    for (auto& f : fr) {
        if ((int)(nxt.size() + f.size()) <= N) { for (int i : f) nxt.push_back(R[i]); continue; }
        std::vector<Vec> F; for (int i : f) F.push_back(R[i].f);
        auto d = crowdingDistance(F); std::vector<int> idx(f.size()); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return d[a] > d[b]; });
        for (size_t k = 0; (int)nxt.size() < N; k++) nxt.push_back(R[f[idx[k]]]);
        break;
    }
    assignRankCrowd(nxt); return nxt;
}
// (14) NSGA-II. With g_deb=true it is the Constrained NSGA-II (21) because cdom implements Deb's rules.
inline std::vector<Ind> run_nsga2(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); assignRankCrowd(pop);
    while (P.nev + N <= budget) {
        auto off = makeOffspring(P, pop, N, [&] { return tourCrowded(pop); });
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); pop = nsga2Select(R, N);
    }
    return rankFront0(pop);
}

// ---- niche counts / sharing helpers (objective space, normalised) ----
inline double shareF(double d, double sigma) { return d < sigma ? 1 - (d / sigma) * (d / sigma) : 0.0; }
inline std::vector<double> nicheCounts(const std::vector<Vec>& V, const std::vector<int>& idx, double sigma) {
    std::vector<double> m(idx.size(), 0.0);
    for (size_t a = 0; a < idx.size(); a++) for (size_t b = 0; b < idx.size(); b++) m[a] += shareF(dist(V[idx[a]], V[idx[b]]), sigma);
    return m;
}
inline int tourMax(const std::vector<double>& fit) {   // binary tournament, higher fitness wins
    int a = irand((int)fit.size()), b = irand((int)fit.size()); return fit[a] >= fit[b] ? a : b;
}

// (10) VEGA: M sub-populations, sub-population k is filled by tournament on objective k only; shuffled; recombined.
inline std::vector<Ind> run_vega(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N);
    while (P.nev + N <= budget) {
        std::vector<Vec> V; for (auto& a : pop) V.push_back(pn(P, a));
        std::vector<Ind> mating; 
        for (int k = 0; k < P.M; k++) for (int s = 0; s < N / P.M; s++) {
            int a = irand(N), b = irand(N); mating.push_back(V[a][k] <= V[b][k] ? pop[a] : pop[b]); }
        while ((int)mating.size() < N) mating.push_back(pop[irand(N)]);
        std::shuffle(mating.begin(), mating.end(), rng);
        pop = makeOffspring(P, mating, N, [&] { return irand(N); });
    }
    return rankFront0(pop);
}
// (11) MOGA: rank = 1 + #dominators, rank-based fitness averaged per rank, fitness sharing, tournament; no elitism.
inline std::vector<Ind> run_moga(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); const double sigma = 0.15;
    while (P.nev + N <= budget) {
        std::vector<Vec> V; for (auto& a : pop) V.push_back(pn(P, a));
        std::vector<int> rk(N, 1);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) if (i != j && cdom(pop[j], pop[i])) rk[i]++;
        int maxr = *std::max_element(rk.begin(), rk.end()); std::vector<int> cnt(maxr + 2, 0); for (int r : rk) cnt[r]++;
        std::vector<double> raw(N), fit(N);
        for (int i = 0; i < N; i++) { double above = 0; for (int r = 1; r < rk[i]; r++) above += cnt[r]; raw[i] = N - above - 0.5 * (cnt[rk[i]] - 1); }
        std::vector<double> sumRaw(maxr + 2, 0), sumSh(maxr + 2, 0), sh(N);
        std::vector<int> all(N); std::iota(all.begin(), all.end(), 0);
        for (int i = 0; i < N; i++) { double m = 0; for (int j = 0; j < N; j++) if (rk[i] == rk[j]) m += shareF(dist(V[i], V[j]), sigma);
            sh[i] = raw[i] / std::max(m, 1e-9); sumRaw[rk[i]] += raw[i]; sumSh[rk[i]] += sh[i]; }
        for (int i = 0; i < N; i++) fit[i] = sh[i] * sumRaw[rk[i]] / std::max(sumSh[rk[i]], 1e-12);
        pop = makeOffspring(P, pop, N, [&] { return tourMax(fit); });
    }
    return rankFront0(pop);
}
// (12) NSGA: non-dominated fronts, dummy fitness shared inside each front, next front starts below previous minimum.
inline std::vector<Ind> run_nsga(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); const double sigma = 0.15;
    while (P.nev + N <= budget) {
        std::vector<Vec> V; for (auto& a : pop) V.push_back(pn(P, a));
        auto fr = fastNonDominatedSort(pop); std::vector<double> fit(N, 0.0); double dummy = N;
        for (auto& f : fr) {
            auto m = nicheCounts(V, f, sigma); double mn = 1e300;
            for (size_t k = 0; k < f.size(); k++) { fit[f[k]] = dummy / m[k]; mn = std::min(mn, fit[f[k]]); }
            dummy = mn * 0.9;
        }
        pop = makeOffspring(P, pop, N, [&] { return tourMax(fit); });
    }
    return rankFront0(pop);
}
// (13) NPGA: Pareto tournament against a random comparison set (10% of pop); tie -> smaller niche count.
inline std::vector<Ind> run_npga(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); const double sigma = 0.15; int tdom = std::max(2, N / 10);
    while (P.nev + N <= budget) {
        std::vector<Vec> V; for (auto& a : pop) V.push_back(pn(P, a));
        std::vector<int> all(N); std::iota(all.begin(), all.end(), 0);
        auto niche = nicheCounts(V, all, sigma);
        auto sel = [&] {
            int a = irand(N), b = irand(N); bool da = false, db = false;
            for (int k = 0; k < tdom; k++) { int c = irand(N); if (cdom(pop[c], pop[a])) da = true; if (cdom(pop[c], pop[b])) db = true; }
            if (da && !db) return b; if (db && !da) return a;
            return niche[a] <= niche[b] ? a : b; };
        pop = makeOffspring(P, pop, N, sel);
    }
    return rankFront0(pop);
}

// (16) SPEA2
inline std::vector<Ind> run_spea2(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P), Na = N; auto pop = initPop(P, N); std::vector<Ind> arch;
    std::vector<double> afit;
    auto fitness = [&](const std::vector<Ind>& U, std::vector<double>& F, std::vector<Vec>& V) {
        int n = (int)U.size(); V.clear(); for (auto& a : U) V.push_back(pn(P, a));
        std::vector<int> S(n, 0), R(n, 0);
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && cdom(U[i], U[j])) S[i]++;
        for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j && cdom(U[j], U[i])) R[i] += S[j];
        int k = (int)std::sqrt((double)n); F.assign(n, 0.0);
        for (int i = 0; i < n; i++) { std::vector<double> d; for (int j = 0; j < n; j++) if (i != j) d.push_back(dist(V[i], V[j]));
            std::sort(d.begin(), d.end()); double sk = d[std::min((int)d.size() - 1, k - 1)]; F[i] = R[i] + 1.0 / (sk + 2.0); }
    };
    while (P.nev + N <= budget) {
        std::vector<Ind> U = pop; U.insert(U.end(), arch.begin(), arch.end());
        std::vector<double> F; std::vector<Vec> V; fitness(U, F, V);
        std::vector<int> nd, dm; for (int i = 0; i < (int)U.size(); i++) (F[i] < 1 ? nd : dm).push_back(i);
        if ((int)nd.size() > Na) {   // truncation: repeatedly delete the point with the smallest nearest-neighbour distance
            while ((int)nd.size() > Na) {
                int worst = -1; std::vector<double> wd;
                for (int a : nd) { std::vector<double> d; for (int b : nd) if (a != b) d.push_back(dist(V[a], V[b])); std::sort(d.begin(), d.end());
                    if (worst < 0 || d < wd) { worst = a; wd = d; } }
                nd.erase(std::find(nd.begin(), nd.end(), worst));
            }
        } else {
            std::sort(dm.begin(), dm.end(), [&](int a, int b) { return F[a] < F[b]; });
            for (size_t k = 0; (int)nd.size() < Na && k < dm.size(); k++) nd.push_back(dm[k]);
        }
        arch.clear(); afit.clear(); for (int i : nd) { arch.push_back(U[i]); afit.push_back(F[i]); }
        auto off = makeOffspring(P, arch, N, [&] { int a = irand((int)arch.size()), b = irand((int)arch.size()); return afit[a] <= afit[b] ? a : b; });
        pop = off;
    }
    return rankFront0(arch);
}

// archive insertion with crowding truncation (used by PESA-II / swarm / bat / ABC / SA)
inline void archiveAdd(std::vector<Ind>& A, const Ind& c, int maxsize) {
    for (auto& a : A) if (cdom(a, c) || (a.f == c.f)) return;
    A.erase(std::remove_if(A.begin(), A.end(), [&](const Ind& a) { return cdom(c, a); }), A.end());
    A.push_back(c);
    while ((int)A.size() > maxsize) {
        std::vector<Vec> F; for (auto& a : A) F.push_back(a.f);
        auto d = crowdingDistance(F); int w = (int)(std::min_element(d.begin(), d.end()) - d.begin()); A.erase(A.begin() + w);
    }
}
// (17) PESA-II: hyper-grid on the archive, selection of a grid box (fewest occupants) then a member inside it.
inline std::vector<Ind> run_pesa2(const Problem& P, long budget) {
    P.nev = 0; const int Npop = 10, Amax = 100, DIV = 16; std::vector<Ind> arch;
    auto pop = initPop(P, Npop); for (auto& a : pop) archiveAdd(arch, a, Amax);
    while (P.nev + Npop <= budget) {
        int n = (int)arch.size(); std::vector<Vec> V; for (auto& a : arch) V.push_back(P.norm(a.f));
        Vec lo(P.M, 1e18), hi(P.M, -1e18); for (auto& v : V) for (int i = 0; i < P.M; i++) { lo[i] = std::min(lo[i], v[i]); hi[i] = std::max(hi[i], v[i]); }
        std::vector<long> box(n);
        for (int k = 0; k < n; k++) { long id = 0; for (int i = 0; i < P.M; i++) { int c = hi[i] - lo[i] < 1e-12 ? 0 : std::min(DIV - 1, (int)((V[k][i] - lo[i]) / (hi[i] - lo[i]) * DIV)); id = id * DIV + c; } box[k] = id; }
        std::vector<int> cnt(n); for (int k = 0; k < n; k++) for (int j = 0; j < n; j++) if (box[j] == box[k]) cnt[k]++;
        auto sel = [&] { int a = irand(n), b = irand(n); return cnt[a] <= cnt[b] ? a : b; };   // region-based: least crowded box wins
        auto off = makeOffspring(P, arch, Npop, sel);
        for (auto& c : off) archiveAdd(arch, c, Amax);
    }
    return arch;
}

// ---- decomposition: (18) MOEA/D, (19) MOEA/D-DE ----
inline std::vector<Ind> run_moead(const Problem& P, long budget, bool de) {
    P.nev = 0; int M = P.M; std::vector<Vec> W;
    if (M == 2) W = dasDennis(2, 99); else W = dasDennis(M, M == 3 ? 13 : 5);
    int N = (int)W.size(); for (auto& w : W) for (double& v : w) v = std::max(v, 1e-6);
    int T = std::max(5, N / 5); double delta = 0.9; int nr = de ? 2 : T;
    std::vector<std::vector<int>> B(N);
    for (int i = 0; i < N; i++) { std::vector<int> idx(N); std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return dist(W[i], W[a]) < dist(W[i], W[b]); }); B[i].assign(idx.begin(), idx.begin() + T); }
    auto pop = initPop(P, N); Vec z(M, 1e18);
    auto upd = [&](const Ind& a) { Vec v = pn(P, a); for (int i = 0; i < M; i++) z[i] = std::min(z[i], std::min(0.0, v[i])); };
    for (auto& a : pop) upd(a);
    auto g = [&](const Ind& a, const Vec& w) { Vec v = pn(P, a); double s = 0; for (int i = 0; i < M; i++) s = std::max(s, w[i] * std::fabs(v[i] - z[i])); return s; };
    while (P.nev + N <= budget) {
        for (int i = 0; i < N; i++) {
            bool nb = urand() < delta; auto pool = [&]() { return nb ? B[i][irand(T)] : irand(N); };
            Vec child;
            if (de) { int r1 = pool(), r2 = pool(), r3 = pool(); child = pop[i].x;
                for (int d = 0; d < P.n; d++) child[d] = pop[r1].x[d] + 0.5 * (pop[r2].x[d] - pop[r3].x[d]); }
            else { Vec c1 = pop[pool()].x, c2 = pop[pool()].x; sbx(c1, c2, P); child = c1; }
            polyMut(child, P); Ind c = P.make(child); upd(c);
            int replaced = 0; std::vector<int> cand = nb ? B[i] : std::vector<int>(); if (!nb) { cand.resize(N); std::iota(cand.begin(), cand.end(), 0); }
            std::shuffle(cand.begin(), cand.end(), rng);
            for (int j : cand) { if (replaced >= nr) break; if (g(c, W[j]) <= g(pop[j], W[j])) { pop[j] = c; replaced++; } }
        }
    }
    return rankFront0(pop);
}

// (15) NSGA-III
inline std::vector<Ind> nsga3Select(const Problem& P, const std::vector<Ind>& R, int N, const std::vector<Vec>& Z) {
    auto fr = fastNonDominatedSort(R); std::vector<int> St, Fl;
    for (auto& f : fr) { if ((int)(St.size() + f.size()) <= N) { St.insert(St.end(), f.begin(), f.end()); if ((int)St.size() == N) break; } else { Fl = f; break; } }
    std::vector<Ind> out; for (int i : St) out.push_back(R[i]);
    if ((int)St.size() == N) return out;
    int M = P.M; std::vector<int> all = St; all.insert(all.end(), Fl.begin(), Fl.end());
    std::vector<Vec> V(R.size()); for (int i : all) V[i] = pn(P, R[i]);
    Vec ideal(M, 1e18); for (int i : all) for (int k = 0; k < M; k++) ideal[k] = std::min(ideal[k], V[i][k]);
    for (int i : all) for (int k = 0; k < M; k++) V[i][k] -= ideal[k];
    std::vector<Vec> E(M);
    for (int k = 0; k < M; k++) { double best = 1e300; int bi = all[0];
        for (int i : all) { double s = 0; for (int j = 0; j < M; j++) s = std::max(s, V[i][j] / (j == k ? 1.0 : 1e6)); if (s < best) { best = s; bi = i; } }
        E[k] = V[bi]; }
    Vec inter(M), sol; bool ok = solveLinear(E, Vec(M, 1.0), sol);
    for (int k = 0; k < M; k++) { if (ok && sol[k] > 1e-9) inter[k] = 1.0 / sol[k]; else { double mx = 1e-9; for (int i : all) mx = std::max(mx, V[i][k]); inter[k] = mx; } }
    for (int i : all) for (int k = 0; k < M; k++) V[i][k] /= std::max(inter[k], 1e-9);
    int nz = (int)Z.size(); std::vector<int> assoc(R.size(), -1); std::vector<double> dd(R.size(), 1e18);
    for (int i : all) for (int j = 0; j < nz; j++) {
        double zn = 0, dot = 0; for (int k = 0; k < M; k++) { zn += Z[j][k]*Z[j][k]; dot += Z[j][k]*V[i][k]; }
        double proj = dot / zn, d2 = 0; for (int k = 0; k < M; k++) { double t = V[i][k] - proj * Z[j][k]; d2 += t*t; }
        d2 = std::sqrt(d2); if (d2 < dd[i]) { dd[i] = d2; assoc[i] = j; } }
    std::vector<int> rho(nz, 0); for (int i : St) rho[assoc[i]]++;
    std::vector<bool> inFl(R.size(), false); for (int i : Fl) inFl[i] = true;
    std::vector<bool> excluded(nz, false);
    while ((int)out.size() < N) {
        int mn = 1 << 30; for (int j = 0; j < nz; j++) if (!excluded[j]) mn = std::min(mn, rho[j]);
        std::vector<int> cands; for (int j = 0; j < nz; j++) if (!excluded[j] && rho[j] == mn) cands.push_back(j);
        int j = cands[irand((int)cands.size())];
        std::vector<int> mem; for (int i : Fl) if (inFl[i] && assoc[i] == j) mem.push_back(i);
        if (mem.empty()) { excluded[j] = true; continue; }
        int pick = mem[irand((int)mem.size())];
        if (rho[j] == 0) for (int i : mem) if (dd[i] < dd[pick]) pick = i;
        out.push_back(R[pick]); inFl[pick] = false; rho[j]++;
    }
    return out;
}
inline std::vector<Ind> run_nsga3(const Problem& P, long budget) {
    P.nev = 0; auto Z = dasDennis(P.M, P.M == 2 ? 99 : (P.M == 3 ? 12 : 5)); int N = popSize(P);
    auto pop = initPop(P, N);
    while (P.nev + N <= budget) {
        auto off = makeOffspring(P, pop, N, [&] { return irand(N); });
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); pop = nsga3Select(P, R, N, Z);
    }
    return rankFront0(pop);
}

// (31) RVEA
inline std::vector<Ind> run_rvea(const Problem& P, long budget) {
    P.nev = 0; int M = P.M; auto V0 = dasDennis(M, M == 2 ? 99 : (M == 3 ? 12 : 5));
    for (auto& v : V0) { double s = 0; for (double t : v) s += t*t; s = std::sqrt(s); for (double& t : v) t /= s; }
    auto Vc = V0; int nv = (int)V0.size(); int N = popSize(P); auto pop = initPop(P, N);
    long maxGen = std::max(1L, budget / N); const double alpha = 2.0, fr = 0.1;
    Vec zmin(M, 1e18), zmax(M, -1e18);
    for (long gen = 0; P.nev + N <= budget; gen++) {
        auto off = makeOffspring(P, pop, N, [&] { return irand((int)pop.size()); });
        auto Q = pop; Q.insert(Q.end(), off.begin(), off.end());
        std::vector<Vec> F; for (auto& a : Q) F.push_back(pn(P, a));
        for (auto& f : F) for (int k = 0; k < M; k++) { zmin[k] = std::min(zmin[k], f[k]); zmax[k] = std::max(zmax[k], f[k]); }
        std::vector<double> gam(nv, 1e18);
        for (int i = 0; i < nv; i++) for (int j = 0; j < nv; j++) if (i != j) { double c = 0; for (int k = 0; k < M; k++) c += Vc[i][k]*Vc[j][k]; gam[i] = std::min(gam[i], std::acos(std::min(1.0, std::max(-1.0, c)))); }
        std::vector<int> best(nv, -1); std::vector<double> bestA(nv, 1e300);
        for (size_t q = 0; q < Q.size(); q++) {
            Vec f = F[q]; double nrm = 0; for (int k = 0; k < M; k++) { f[k] -= zmin[k]; nrm += f[k]*f[k]; } nrm = std::sqrt(nrm);
            int bj = 0; double bc = -2; for (int j = 0; j < nv; j++) { double c = 0; for (int k = 0; k < M; k++) c += f[k]*Vc[j][k]; c /= std::max(nrm, 1e-12); if (c > bc) { bc = c; bj = j; } }
            double theta = std::acos(std::min(1.0, std::max(-1.0, bc)));
            double apd = (1 + M * std::pow((double)gen / maxGen, alpha) * theta / std::max(gam[bj], 1e-9)) * nrm;
            if (apd < bestA[bj]) { bestA[bj] = apd; best[bj] = (int)q; } }
        std::vector<Ind> np; for (int j = 0; j < nv; j++) if (best[j] >= 0) np.push_back(Q[best[j]]);
        pop = np; if (pop.size() < 2) pop = initPop(P, 2);
        if ((gen + 1) % std::max(1L, (long)(fr * maxGen)) == 0) {   // reference vector adaptation
            for (int j = 0; j < nv; j++) { double s = 0; for (int k = 0; k < M; k++) { Vc[j][k] = V0[j][k] * std::max(zmax[k] - zmin[k], 1e-9); s += Vc[j][k]*Vc[j][k]; } s = std::sqrt(s); for (int k = 0; k < M; k++) Vc[j][k] /= s; } }
    }
    return rankFront0(pop);
}

// (32) IBEA with the additive epsilon indicator. verbose prints every number (used by the hand-trace check).
inline double ibeaEps(const Vec& a, const Vec& b) { double m = -1e300; for (size_t i = 0; i < a.size(); i++) m = std::max(m, a[i] - b[i]); return m; }
inline void ibeaFitness(const std::vector<Vec>& F, double kappa, std::vector<double>& fit, std::vector<std::vector<double>>& I, double& c) {
    int n = (int)F.size(); I.assign(n, std::vector<double>(n, 0.0)); c = 1e-12;
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) { I[i][j] = ibeaEps(F[i], F[j]); c = std::max(c, std::fabs(I[i][j])); }
    fit.assign(n, 0.0);
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) if (i != j) fit[i] += -std::exp(-I[j][i] / (c * kappa));
}
inline void ibeaReduce(const Problem& P, std::vector<Ind>& R, int N, double kappa, bool verbose = false) {
    std::vector<Vec> F; for (auto& a : R) F.push_back(pn(P, a));
    std::vector<double> fit; std::vector<std::vector<double>> I; double c;
    ibeaFitness(F, kappa, fit, I, c);
    std::vector<bool> alive(R.size(), true); int cnt = (int)R.size();
    if (verbose) { std::cout << "  scaling c = " << c << ", kappa = " << kappa << "\n  initial fitness:"; for (double f : fit) std::cout << " " << f; std::cout << "\n"; }
    while (cnt > N) {
        int w = -1; for (size_t i = 0; i < R.size(); i++) if (alive[i] && (w < 0 || fit[i] < fit[w])) w = (int)i;
        if (verbose) std::cout << "  remove index " << w << " (fitness " << fit[w] << ")\n";
        alive[w] = false; cnt--;
        for (size_t i = 0; i < R.size(); i++) if (alive[i]) fit[i] += std::exp(-I[w][i] / (c * kappa));
    }
    std::vector<Ind> out; for (size_t i = 0; i < R.size(); i++) if (alive[i]) { R[i].fit = fit[i]; out.push_back(R[i]); } R = out;
}
inline std::vector<Ind> run_ibea(const Problem& P, long budget, double kappa = 0.05) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N);
    auto setFit = [&](std::vector<Ind>& p) { std::vector<Vec> F; for (auto& a : p) F.push_back(pn(P, a)); std::vector<double> f; std::vector<std::vector<double>> I; double c; ibeaFitness(F, kappa, f, I, c); for (size_t i = 0; i < p.size(); i++) p[i].fit = f[i]; };
    setFit(pop);
    while (P.nev + N <= budget) {
        auto off = makeOffspring(P, pop, N, [&] { int a = irand((int)pop.size()), b = irand((int)pop.size()); return pop[a].fit >= pop[b].fit ? a : b; });
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); ibeaReduce(P, R, N, kappa); pop = R;
    }
    return rankFront0(pop);
}

// ---- hypervolume contributions (exact in 2-D, Monte-Carlo for M>=3) ----
inline std::vector<double> hvContrib(const std::vector<Vec>& F, const Vec& ref, int samples = 3000) {
    int n = (int)F.size(), M = (int)ref.size(); std::vector<double> c(n, 0.0);
    if (n == 0) return c;
    if (M == 2) {
        std::vector<int> idx; for (int i = 0; i < n; i++) if (F[i][0] < ref[0] && F[i][1] < ref[1]) idx.push_back(i);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return F[a][0] < F[b][0]; });
        for (size_t k = 0; k < idx.size(); k++) { double nx = k + 1 < idx.size() ? F[idx[k+1]][0] : ref[0], pv = k > 0 ? F[idx[k-1]][1] : ref[1];
            c[idx[k]] = std::max(0.0, (nx - F[idx[k]][0]) * (pv - F[idx[k]][1])); }
        return c;
    }
    Vec lo(M, 1e18); bool any = false;
    for (auto& f : F) { bool ok = true; for (int i = 0; i < M; i++) if (f[i] >= ref[i]) ok = false; if (ok) { any = true; for (int i = 0; i < M; i++) lo[i] = std::min(lo[i], f[i]); } }
    if (!any) return c;
    double vol = 1; for (int i = 0; i < M; i++) vol *= ref[i] - lo[i];
    Vec s(M);
    for (int t = 0; t < samples; t++) {
        for (int i = 0; i < M; i++) s[i] = lo[i] + urand() * (ref[i] - lo[i]);
        int who = -1, cnt = 0;
        for (int j = 0; j < n; j++) { bool d = true; for (int i = 0; i < M; i++) if (F[j][i] > s[i]) { d = false; break; } if (d) { cnt++; who = j; if (cnt > 1) break; } }
        if (cnt == 1) c[who] += vol / samples;
    }
    return c;
}
// (33) SMS-EMOA: steady state, delete the member of the worst front with the smallest HV contribution.
inline std::vector<Ind> run_smsemoa(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); assignRankCrowd(pop); Vec ref(P.M, 1.1);
    while (P.nev + 1 <= budget) {
        Vec c1 = pop[irand(N)].x, c2 = pop[irand(N)].x; sbx(c1, c2, P); polyMut(c1, P);
        pop.push_back(P.make(c1));
        auto fr = fastNonDominatedSort(pop); auto& last = fr.back();
        std::vector<Vec> F; for (int i : last) F.push_back(pn(P, pop[i]));
        auto c = hvContrib(F, ref); int w = (int)(std::min_element(c.begin(), c.end()) - c.begin());
        pop.erase(pop.begin() + last[w]);
    }
    return rankFront0(pop);
}
// (34) HypE: fitness = Monte-Carlo estimate of hypervolume contribution with weights alpha_j (Bader & Zitzler).
inline std::vector<double> hypeFitness(const std::vector<Vec>& F, const Vec& ref, int k, int samples) {
    int n = (int)F.size(), M = (int)ref.size(); std::vector<double> fit(n, 0.0), alpha(n + 2, 0.0);
    for (int j = 1; j <= std::min(k, n); j++) { double a = 1.0 / j; for (int i = 1; i <= j - 1; i++) a *= (double)(k - i) / (double)(n - i); alpha[j] = a; }
    Vec lo(M, 1e18); for (auto& f : F) for (int i = 0; i < M; i++) lo[i] = std::min(lo[i], f[i]);
    Vec s(M);
    for (int t = 0; t < samples; t++) {
        for (int i = 0; i < M; i++) s[i] = lo[i] + urand() * (ref[i] - lo[i]);
        std::vector<int> dm; for (int j = 0; j < n; j++) { bool d = true; for (int i = 0; i < M; i++) if (F[j][i] > s[i]) { d = false; break; } if (d) dm.push_back(j); }
        if (!dm.empty() && (int)dm.size() <= k) for (int j : dm) fit[j] += alpha[dm.size()];
    }
    return fit;
}
inline std::vector<Ind> run_hype(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); Vec ref(P.M, 1.1); const int S = 1500;
    while (P.nev + N <= budget) {
        std::vector<Vec> F; for (auto& a : pop) F.push_back(pn(P, a));
        auto fr0 = fastNonDominatedSort(pop); std::vector<int> rk(N); for (size_t r = 0; r < fr0.size(); r++) for (int i : fr0[r]) rk[i] = (int)r;
        auto hf = hypeFitness(F, ref, N, S);
        auto off = makeOffspring(P, pop, N, [&] { int a = irand(N), b = irand(N); if (rk[a] != rk[b]) return rk[a] < rk[b] ? a : b; return hf[a] >= hf[b] ? a : b; });
        auto R = pop; R.insert(R.end(), off.begin(), off.end());
        auto fr = fastNonDominatedSort(R); std::vector<Ind> nxt; size_t fi = 0;
        while (fi < fr.size() && (int)(nxt.size() + fr[fi].size()) <= N) { for (int i : fr[fi]) nxt.push_back(R[i]); fi++; }
        if ((int)nxt.size() < N && fi < fr.size()) {
            std::vector<int> cur = fr[fi]; int need = (int)cur.size() - (N - (int)nxt.size());
            for (int rem = 0; rem < need; rem++) {   // iterative removal, fitness recomputed with k = remaining removals
                std::vector<Vec> G; for (int i : cur) G.push_back(pn(P, R[i]));
                auto f = hypeFitness(G, ref, need - rem, S);
                int w = (int)(std::min_element(f.begin(), f.end()) - f.begin()); cur.erase(cur.begin() + w); }
            for (int i : cur) nxt.push_back(R[i]);
        }
        pop = nxt;
    }
    return rankFront0(pop);
}

// (35) R-NSGA-II: NSGA-II where the last front is ordered by distance to the decision-maker's reference point,
// with eps-clearing so that neighbours of a chosen point are delayed.
inline std::vector<Ind> run_rnsga2(const Problem& P, long budget, double epsClear = 0.001) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); assignRankCrowd(pop);
    Vec ref = P.prefRef.empty() ? Vec(P.M, 0.0) : P.prefRef;
    auto select = [&](const std::vector<Ind>& R) {
        auto fr = fastNonDominatedSort(R); std::vector<Ind> nxt;
        for (auto& f : fr) {
            if ((int)(nxt.size() + f.size()) <= N) { for (int i : f) nxt.push_back(R[i]); continue; }
            std::vector<Vec> V; for (int i : f) V.push_back(pn(P, R[i]));
            std::vector<double> d; for (auto& v : V) d.push_back(dist(v, ref));
            std::vector<int> pending(f.size()), order; std::iota(pending.begin(), pending.end(), 0);
            while (!pending.empty()) {
                std::vector<int> delayed, rest = pending; std::sort(rest.begin(), rest.end(), [&](int a, int b) { return d[a] < d[b]; });
                std::vector<bool> gone(f.size(), false);
                for (int a : rest) { if (gone[a]) continue; order.push_back(a);
                    for (int b : rest) if (b != a && !gone[b] && dist(V[a], V[b]) < epsClear) { gone[b] = true; delayed.push_back(b); } }
                pending = delayed;
            }
            for (size_t k = 0; (int)nxt.size() < N; k++) nxt.push_back(R[f[order[k]]]);
            break;
        }
        assignRankCrowd(nxt); return nxt; };
    while (P.nev + N <= budget) {
        auto off = makeOffspring(P, pop, N, [&] { return tourCrowded(pop); });
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); pop = select(R);
    }
    return rankFront0(pop);
}

// (30) MOMA: NSGA-II global search + hill-climbing local search (random-weight scalarisation) on the best few.
inline std::vector<Ind> run_moma(const Problem& P, long budget) {
    P.nev = 0; int N = popSize(P); auto pop = initPop(P, N); assignRankCrowd(pop); const int nLS = 5, steps = 6;
    while (P.nev + N + nLS * steps <= budget) {
        auto off = makeOffspring(P, pop, N, [&] { return tourCrowded(pop); });
        std::vector<Ind> improved;
        for (int l = 0; l < nLS; l++) {
            Ind cur = off[irand(N)]; Vec w(P.M); double s = 0; for (double& t : w) { t = urand(); s += t; } for (double& t : w) t /= s;
            auto sc = [&](const Ind& a) { Vec v = pn(P, a); double z = 0; for (int i = 0; i < P.M; i++) z += w[i] * v[i]; return z; };
            double cs = sc(cur);
            for (int t = 0; t < steps; t++) { Vec x = cur.x; for (int i = 0; i < P.n; i++) if (urand() < 0.5) x[i] += 0.02 * (P.hi[i] - P.lo[i]) * nrand();
                Ind c = P.make(x); double v = sc(c); if (v < cs) { cur = c; cs = v; } }
            improved.push_back(cur); }
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); R.insert(R.end(), improved.begin(), improved.end()); pop = nsga2Select(R, N);
    }
    return rankFront0(pop);
}
}  // namespace moo
