// swarm.hpp - 22 MOPSO, 23 MOGWO, 24 MOABC, 26 MOSA, 27 MODE, 28 MO-TLBO, 29 MOBA
#pragma once
#include "evo.hpp"
namespace moo {

inline int pickLeader(const std::vector<Ind>& A) {   // crowding-based leader: less crowded archive member preferred
    std::vector<Vec> F; for (auto& a : A) F.push_back(a.f);
    auto d = crowdingDistance(F); int x = irand((int)A.size()), y = irand((int)A.size()); return d[x] >= d[y] ? x : y;
}
// (22) MOPSO
inline std::vector<Ind> run_mopso(const Problem& P, long budget) {
    P.nev = 0; int N = 100, Amax = 100; auto sw = initPop(P, N); auto pb = sw; std::vector<Vec> vel(N, Vec(P.n, 0.0));
    std::vector<Ind> arch; for (auto& a : sw) archiveAdd(arch, a, Amax);
    long maxIt = std::max(1L, budget / N); long it = 0;
    while (P.nev + N <= budget) {
        double w = 0.5, c1 = 1.5, c2 = 1.5;
        for (int i = 0; i < N; i++) {
            const Ind& g = arch[pickLeader(arch)];
            Vec x = sw[i].x;
            for (int d = 0; d < P.n; d++) {
                vel[i][d] = w * vel[i][d] + c1 * urand() * (pb[i].x[d] - x[d]) + c2 * urand() * (g.x[d] - x[d]);
                x[d] += vel[i][d];
                if (x[d] < P.lo[d] || x[d] > P.hi[d]) { vel[i][d] = -vel[i][d]; x[d] = std::min(P.hi[d], std::max(P.lo[d], x[d])); } }
            if (urand() < std::pow(1.0 - (double)it / maxIt, 5.0)) polyMut(x, P);   // shrinking mutation probability
            Ind c = P.make(x); sw[i] = c;
            if (cdom(c, pb[i]) || (!cdom(pb[i], c) && urand() < 0.5)) pb[i] = c;
            archiveAdd(arch, c, Amax);
        }
        it++;
    }
    return arch;
}
// (23) MOGWO
inline std::vector<Ind> run_mogwo(const Problem& P, long budget) {
    P.nev = 0; int N = 100, Amax = 100; auto pop = initPop(P, N); std::vector<Ind> arch; for (auto& a : pop) archiveAdd(arch, a, Amax);
    long maxIt = std::max(1L, budget / N), it = 0;
    while (P.nev + N <= budget) {
        double a = 2.0 - 2.0 * (double)it / maxIt;
        for (int i = 0; i < N; i++) {
            int ia = pickLeader(arch), ib = pickLeader(arch), id = pickLeader(arch);   // alpha, beta, delta from the archive
            const Ind* L[3] = {&arch[ia], &arch[ib], &arch[id]}; Vec x(P.n, 0.0);
            for (int d = 0; d < P.n; d++) { for (int k = 0; k < 3; k++) { double A = 2 * a * urand() - a, C = 2 * urand();
                    x[d] += (L[k]->x[d] - A * std::fabs(C * L[k]->x[d] - pop[i].x[d])) / 3.0; } }
            Ind c = P.make(x); pop[i] = c; archiveAdd(arch, c, Amax);
        }
        it++;
    }
    return arch;
}
// (24) MOABC: employed, onlooker, scout phases; greedy choice by dominance.
inline std::vector<Ind> run_moabc(const Problem& P, long budget) {
    P.nev = 0; int SN = 50, Amax = 100, limit = 20; auto food = initPop(P, SN); std::vector<int> trial(SN, 0);
    std::vector<Ind> arch; for (auto& a : food) archiveAdd(arch, a, Amax);
    auto tryMove = [&](int i) {
        int k; do k = irand(SN); while (k == i); Vec v = food[i].x;
        for (int d = 0; d < P.n; d++) if (urand() < 0.5 || d == 0) v[d] += (urand() * 2 - 1) * (food[i].x[d] - food[k].x[d]);
        Ind c = P.make(v); archiveAdd(arch, c, Amax);
        if (cdom(c, food[i]) || (!cdom(food[i], c) && urand() < 0.5)) { food[i] = c; trial[i] = 0; } else trial[i]++; };
    while (P.nev + 2 * SN <= budget) {
        for (int i = 0; i < SN; i++) tryMove(i);
        auto fr = fastNonDominatedSort(food); std::vector<double> p(SN); double tot = 0;
        for (size_t r = 0; r < fr.size(); r++) for (int i : fr[r]) { p[i] = 1.0 / (1 + r); tot += p[i]; }
        for (int cnt = 0, i = 0; cnt < SN; i = (i + 1) % SN) if (urand() < p[i] / tot * SN * 0.5) { tryMove(i); cnt++; }
        for (int i = 0; i < SN; i++) if (trial[i] > limit && P.nev < budget) { food[i] = randInd(P); trial[i] = 0; archiveAdd(arch, food[i], Amax); }
    }
    return arch;
}
// (26) MOSA: archive-based multi-objective SA (Suppapitnarm style), geometric cooling, several chains.
inline std::vector<Ind> run_mosa(const Problem& P, long budget) {
    P.nev = 0; const int Amax = 100, chains = 20; std::vector<Ind> arch; long steps = std::max(10L, budget / chains - 1);
    const double T0 = 0.5, T1 = 1e-3; double cool = std::pow(T1 / T0, 1.0 / steps);
    for (int c = 0; c < chains && P.nev < budget; c++) {
        Ind cur = (arch.empty() || urand() < 0.5) ? randInd(P) : arch[irand((int)arch.size())]; archiveAdd(arch, cur, Amax); double T = T0;
        for (long s = 0; s < steps && P.nev < budget; s++) {
            Vec x = cur.x; double sc = 0.02 + 0.3 * std::sqrt(T / T0);
            int forced = irand(P.n); for (int d = 0; d < P.n; d++) if (d == forced || urand() < 1.0 / P.n) x[d] += sc * (P.hi[d] - P.lo[d]) * nrand();
            Ind nw = P.make(x); Vec a = pn(P, cur), b = pn(P, nw); bool dominatedByArch = false;
            for (auto& q : arch) if (cdom(q, nw) || q.f == nw.f) { dominatedByArch = true; break; }
            if (!dominatedByArch) { archiveAdd(arch, nw, Amax); cur = nw; }
            else { double e = 0; for (int i = 0; i < P.M; i++) e += (b[i] - a[i]) / T; if (urand() < std::min(1.0, std::exp(-e))) cur = nw; }
            T *= cool;
        }
    }
    return arch;
}
// (27) MODE: DE/rand/1/bin, trial replaces target only if it dominates it.
inline std::vector<Ind> run_mode(const Problem& P, long budget) {
    P.nev = 0; int N = 100; auto pop = initPop(P, N); const double F = 0.5, CR = 0.9;
    while (P.nev + N <= budget) {
        for (int i = 0; i < N; i++) {
            int a, b, c; do a = irand(N); while (a == i); do b = irand(N); while (b == i || b == a); do c = irand(N); while (c == i || c == a || c == b);
            Vec t = pop[i].x; int jr = irand(P.n);
            for (int d = 0; d < P.n; d++) if (urand() < CR || d == jr) t[d] = pop[a].x[d] + F * (pop[b].x[d] - pop[c].x[d]);
            Ind tr = P.make(t); if (cdom(tr, pop[i])) pop[i] = tr; }
    }
    return rankFront0(pop);
}
// (28) MO-TLBO: teacher phase (teacher from the non-dominated set) + learner phase, teaching factor shrinking 2 -> 1.
inline std::vector<Ind> run_motlbo(const Problem& P, long budget) {
    P.nev = 0; int N = 50, Amax = 100; auto pop = initPop(P, N); std::vector<Ind> arch; for (auto& a : pop) archiveAdd(arch, a, Amax);
    long maxIt = std::max(1L, budget / (2 * N)), it = 0;
    auto accept = [&](Ind& old, const Ind& c) { archiveAdd(arch, c, Amax); if (cdom(c, old) || (!cdom(old, c) && urand() < 0.5)) old = c; };
    while (P.nev + 2 * N <= budget) {
        double TF = 2.0 - (double)it / maxIt; Vec mean(P.n, 0.0); for (auto& a : pop) for (int d = 0; d < P.n; d++) mean[d] += a.x[d] / N;
        auto nd = rankFront0(pop); const Ind& teacher = nd[irand((int)nd.size())];
        for (int i = 0; i < N; i++) { Vec x = pop[i].x; for (int d = 0; d < P.n; d++) x[d] += urand() * (teacher.x[d] - TF * mean[d]); accept(pop[i], P.make(x)); }
        for (int i = 0; i < N; i++) { int j; do j = irand(N); while (j == i); Vec x = pop[i].x; bool iBetter = cdom(pop[i], pop[j]);
            for (int d = 0; d < P.n; d++) x[d] += urand() * (iBetter ? pop[i].x[d] - pop[j].x[d] : pop[j].x[d] - pop[i].x[d]);
            accept(pop[i], P.make(x)); }
        it++;
    }
    return arch;
}
// (29) MOBA: frequency-driven velocity towards an archive leader, local walk with probability 1-r, loudness A, pulse rate r.
inline std::vector<Ind> run_moba(const Problem& P, long budget) {
    P.nev = 0; int N = 100, Amax = 100; auto bat = initPop(P, N); std::vector<Vec> vel(N, Vec(P.n, 0.0));
    Vec A(N, 0.9), r0(N, 0.5), r(N, 0.0); std::vector<Ind> arch; for (auto& a : bat) archiveAdd(arch, a, Amax);
    const double alpha = 0.9, gam = 0.9; int it = 1;
    while (P.nev + N <= budget) {
        double Amean = 0; for (double a : A) Amean += a / N;
        for (int i = 0; i < N; i++) {
            const Ind& L = arch[pickLeader(arch)]; double fr = urand(); Vec x = bat[i].x;
            if (urand() > r[i]) { const Ind& q = arch[irand((int)arch.size())]; for (int d = 0; d < P.n; d++) x[d] = q.x[d] + 0.05 * Amean * (P.hi[d] - P.lo[d]) * nrand(); }
            else for (int d = 0; d < P.n; d++) { vel[i][d] += (L.x[d] - bat[i].x[d]) * fr; x[d] += vel[i][d]; }
            Ind c = P.make(x); archiveAdd(arch, c, Amax);
            if (urand() < A[i] && (cdom(c, bat[i]) || !cdom(bat[i], c))) { bat[i] = c; A[i] *= alpha; r[i] = r0[i] * (1 - std::exp(-gam * it)); }
        }
        it++;
    }
    return arch;
}
}  // namespace moo
