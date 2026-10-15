// main.cpp - driver: registry of the 37 algorithms, metrics, unit tests, experiments for Tasks A-D.
#include "swarm.hpp"
#include "aco.hpp"
#include "grad.hpp"
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
using namespace moo;
namespace fs = std::filesystem;

struct Algo { int id; std::string name, family; std::vector<std::string> probs; };
static const std::vector<Algo> CATALOGUE = {
 {1,"Weighted Sum","Classical",{"P1","P3","P4"}}, {2,"eps-Constraint","Classical",{"P2","P3","P5"}},
 {3,"Benson's Method","Classical",{"P1","P7"}}, {4,"Weight Metric (Tchebycheff)","Classical",{"P1","P3"}},
 {5,"Value Function","Classical",{"P4","P14"}}, {6,"Goal Programming","Classical",{"P11","P7"}},
 {7,"Lexicographic","Classical",{"P11","P7"}}, {8,"TOPSIS","MCDM ranking",{"P12"}}, {9,"VIKOR","MCDM ranking",{"P12"}},
 {10,"VEGA","Non-elitist EA",{"P2","P3"}}, {11,"MOGA","Non-elitist EA",{"P3","P5"}}, {12,"NSGA","Non-elitist EA",{"P5","P7"}},
 {13,"NPGA","Non-elitist EA",{"P3","P5"}}, {14,"NSGA-II","Elitist EA",{"P3","P5","P6"}}, {15,"NSGA-III","Elitist many-obj",{"P10_3","P10_5"}},
 {16,"SPEA2","Elitist EA",{"P4","P5"}}, {17,"PESA-II","Elitist EA",{"P3","P5"}}, {18,"MOEA/D","Decomposition",{"P3","P10_3"}},
 {19,"MOEA/D-DE","Decomposition",{"P3","P5"}}, {20,"MOEA/D-ACO","Decomposition (discrete)",{"P8"}},
 {21,"Constrained NSGA-II","Constrained",{"P5","P6","P7"}}, {22,"MOPSO","Swarm",{"P4","P5"}}, {23,"MOGWO","Swarm",{"P5"}},
 {24,"MOABC","Swarm",{"P3","P6"}}, {25,"MOACO","Swarm (discrete)",{"P8"}}, {26,"MOSA","Physics-inspired",{"P3","P5"}},
 {27,"MODE","Evolutionary",{"P3","P5"}}, {28,"MO-TLBO","Human-inspired",{"P3","P5"}}, {29,"MOBA","Swarm",{"P3","P4"}},
 {30,"MOMA","Hybrid",{"P5"}}, {31,"RVEA","Many-objective",{"P10_5"}}, {32,"IBEA","Indicator",{"P3","P14"}},
 {33,"SMS-EMOA","Indicator",{"P3","P14"}}, {34,"HypE","Indicator",{"P10_3"}}, {35,"R-NSGA-II","Preference",{"P5","P14"}},
 {36,"MGDA","Gradient",{"P13"}}, {37,"GradNorm","Gradient",{"P13"}}};

static std::vector<Ind> dispatch(int id, const Problem& P, long budget) {
    g_deb = (id != 14);   // NSGA-II (14) is the plain version: no feasibility rules. Everything else uses Deb rules / penalties.
    std::vector<Ind> r;
    switch (id) {
        case 1: r = run_weighted_sum(P, budget); break;
        case 2: r = run_eps_constraint(P, budget); break;
        case 3: r = run_benson(P, budget); break;
        case 4: r = run_weight_metric(P, budget, 0); break;
        case 5: r = run_value_function(P, budget); break;
        case 6: { Vec t = P.name == "P11" ? Vec{-2400, -600} : Vec(P.fmin); Vec a1(P.M, 1.0), a2(P.M);
                  for (int i = 0; i < P.M; i++) a2[i] = 1.0 / std::max(1e-9, std::fabs(t[i]) > 1e-9 ? std::fabs(t[i]) : (P.fmax[i] - P.fmin[i]));
                  r.push_back(goal_programming(P, t, a1, budget / 2)); r.push_back(goal_programming(P, t, a2, budget / 2)); break; }
        case 7: r = lexicographic(P, {0, 1}, 1e-3, budget); break;
        case 10: r = run_vega(P, budget); break;
        case 11: r = run_moga(P, budget); break;
        case 12: r = run_nsga(P, budget); break;
        case 13: r = run_npga(P, budget); break;
        case 14: r = run_nsga2(P, budget); break;
        case 15: r = run_nsga3(P, budget); break;
        case 16: r = run_spea2(P, budget); break;
        case 17: r = run_pesa2(P, budget); break;
        case 18: r = run_moead(P, budget, false); break;
        case 19: r = run_moead(P, budget, true); break;
        case 21: r = run_nsga2(P, budget); break;
        case 22: r = run_mopso(P, budget); break;
        case 23: r = run_mogwo(P, budget); break;
        case 24: r = run_moabc(P, budget); break;
        case 26: r = run_mosa(P, budget); break;
        case 27: r = run_mode(P, budget); break;
        case 28: r = run_motlbo(P, budget); break;
        case 29: r = run_moba(P, budget); break;
        case 30: r = run_moma(P, budget); break;
        case 31: r = run_rvea(P, budget); break;
        case 32: r = run_ibea(P, budget); break;
        case 33: r = run_smsemoa(P, budget); break;
        case 34: r = run_hype(P, budget); break;
        case 35: r = run_rnsga2(P, budget); break;
        default: throw std::runtime_error("no generic runner for algorithm " + std::to_string(id));
    }
    g_deb = true; return r;
}

struct Met { double hv = 0, igd = NAN, sp = 0, feas = 1; int n = 0; };
static Met metrics(const Problem& P, const std::vector<Ind>& A, const std::vector<Vec>* refFront = nullptr) {
    Met m; int total = (int)A.size(); std::vector<Vec> F;
    for (auto& a : A) if (a.cv <= 1e-9) F.push_back(P.norm(a.f));
    m.feas = total ? (double)F.size() / total : 0; F = nondominatedSubset(F); m.n = (int)F.size();
    m.hv = hypervolume(F, Vec(P.M, 1.1)); m.sp = spacing(F);
    const std::vector<Vec>* T = refFront ? refFront : (P.trueFront.empty() ? nullptr : &P.trueFront);
    if (T) { std::vector<Vec> Tn; for (auto& t : *T) Tn.push_back(P.norm(t)); m.igd = igd(Tn, F); }
    return m;
}
static void writeFront(const std::string& path, const Problem& P, const std::vector<Ind>& A) {
    std::ofstream o(path); o << "f1";
    for (int i = 1; i < P.M; i++) o << ",f" << i + 1; o << ",cv\n";
    for (auto& a : A) { for (int i = 0; i < P.M; i++) o << a.f[i] << ","; o << a.cv << "\n"; }
}
static void writeTrue(const Problem& P) { std::ofstream o("results/true_" + P.name + ".csv"); for (auto& f : P.trueFront) { for (int i = 0; i < P.M; i++) o << f[i] << (i + 1 < P.M ? "," : "\n"); } }
static double ms(std::chrono::steady_clock::time_point a) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count(); }
static double mean(const std::vector<double>& v) { double s = 0; int c = 0; for (double x : v) if (!std::isnan(x)) { s += x; c++; } return c ? s / c : NAN; }
static double sd(const std::vector<double>& v) { double m = mean(v), s = 0; int c = 0; for (double x : v) if (!std::isnan(x)) { s += (x-m)*(x-m); c++; } return c > 1 ? std::sqrt(s / (c - 1)) : 0; }

// ---------------- Task D3 / P8 helper ----------------
static TSP g_tsp; static std::vector<Vec> g_tspRef; static bool g_tspReady = false;
static void ensureTSP() {
    if (g_tspReady) return; g_tsp = TSP::make(15, 7); g_tspRef = tspReferenceFront(g_tsp); g_tspReady = true;
    std::ofstream o("results/P8_coordinates.csv"); o << "city,x,y\n"; for (int i = 0; i < g_tsp.n; i++) o << i << "," << g_tsp.xy[i][0] << "," << g_tsp.xy[i][1] << "\n";
    std::ofstream r("results/P8_reference_front.csv"); for (auto& f : g_tspRef) r << f[0] << "," << f[1] << "\n";
}
static Met tspMetrics(const std::vector<Tour>& A) {
    std::vector<Vec> F; for (auto& t : A) F.push_back(g_tsp.normf(t.f)); F = nondominatedSubset(F); Met m; m.n = (int)F.size();
    m.hv = hv2d(F, {1.1, 1.1}); m.sp = spacing(F); std::vector<Vec> T; for (auto& f : g_tspRef) T.push_back(g_tsp.normf(f)); m.igd = igd(T, F); return m;
}
static std::vector<Tour> runTSP(int id, long budget) {   // id: 20, 25, 14 (perm NSGA-II), 30 (perm MOMA)
    if (id == 20) return run_moeadaco(g_tsp, budget); if (id == 25) return run_moaco(g_tsp, budget);
    return run_nsga2_perm(g_tsp, budget, id == 30);
}

// ---------------- unit tests ----------------
static int g_fail = 0;
#define CHECK(cond, msg) do { bool ok_ = (cond); std::cout << (ok_ ? "  [PASS] " : "  [FAIL] ") << msg << "\n"; if (!ok_) g_fail++; } while (0)
static bool near(double a, double b, double t = 1e-9) { return std::fabs(a - b) <= t; }
static void unitTests() {
    std::cout << "== moo_utils tests ==\n"; seed(1);
    CHECK(dominates({1,2},{2,3}) && !dominates({1,3},{2,2}) && !dominates({1,2},{1,2}), "dominance test");
    std::vector<Vec> F = {{1,5},{2,3},{3,4},{4,1},{5,5}}; auto fr = fastNonDominatedSortF(F);
    CHECK(fr.size() == 3 && fr[0].size() == 3, "non-dominated sort on 5 points: 3 fronts, first has 3 points");
    auto cd = crowdingDistance({{1,5},{2,3},{4,1}}); CHECK(cd[0] > 1e17 && cd[2] > 1e17, "crowding: boundary points infinite");
    CHECK(near(cd[1], (4.0-1)/3 + (5.0-1)/4, 1e-12), "crowding middle point = (4-1)/3 + (5-1)/4");
    CHECK(near(hv2d({{1,3},{2,2},{3,1}}, {4,4}), 6.0), "hv2d = 6 for (1,3),(2,2),(3,1), ref (4,4)");
    CHECK(near(hv2d({{1,1}}, {2,2}), 1.0) && near(hv2d({{3,3}}, {2,2}), 0.0), "hv2d single point and out-of-reference point");
    CHECK(near(hvMC({{1,3},{2,2},{3,1}}, {4,4}, 200000), 6.0, 0.06), "hvMC close to the exact value 6");
    CHECK(near(igd({{0,0},{1,1}}, {{0,0}}), std::sqrt(2.0) / 2), "IGD = (0 + sqrt2)/2");
    CHECK(near(spacing({{0,0},{1,1},{2,2},{3,3}}), 0.0), "spacing of an evenly spread front = 0");
    CHECK((int)dasDennis(3, 12).size() == 91 && (int)dasDennis(5, 6).size() == 210, "Das-Dennis counts C(14,2)=91 and C(10,4)=210");
    Vec sol; CHECK(solveLinear({{2,1},{1,3}}, {3,5}, sol) && near(sol[0], 0.8) && near(sol[1], 1.4), "Gaussian elimination");
    std::cout << "== Deb constrained dominance on P6 ==\n";
    Problem p6 = P6(); g_deb = true;
    Ind A = p6.make({0,0}), B = p6.make({5,3}), C = p6.make({0,3}), D = p6.make({0.5,3});
    std::cout << "  A=(0,0) f=(" << A.f[0] << "," << A.f[1] << ") cv=" << A.cv << " | B=(5,3) f=(" << B.f[0] << "," << B.f[1] << ") cv=" << B.cv << "\n  C=(0,3) cv=" << C.cv << " | D=(0.5,3) cv=" << D.cv << "\n";
    CHECK(near(A.f[0], 0) && near(A.f[1], 50) && near(B.f[0], 136) && near(B.f[1], 4) && A.cv == 0 && B.cv == 0, "P6 sanity: (0,0)->(0,50), (5,3)->(136,4), both feasible");
    CHECK(cdom(A, C) && !cdom(C, A), "rule 1: feasible beats infeasible");
    CHECK(cdom(D, C) && !cdom(C, D), "rule 2: both infeasible -> smaller violation wins");
    CHECK(!cdom(A, B) && !cdom(B, A), "rule 3: both feasible -> Pareto dominance (A,B incomparable)");
    std::cout << "== problem sanity checks from the sheet ==\n";
    { Problem p = P1(); Ind a = p.make({5}); CHECK(near(a.f[0], 4) && near(a.f[1], 4), "P1: x=5 -> (4,4)"); }
    { Problem p = P2(); Ind a = p.make({2}); Ind b = p.make({10}); CHECK(near(a.f[0], 50) && near(a.f[1], 4) && near(b.f[0], 10) && near(b.f[1], 100), "P2: x=2 -> (50,4), x=10 -> (10,100)"); }
    { Problem p = P4(); Ind a = p.make({0,0,1,0}); CHECK(near(a.f[0], 0.09) && near(a.f[1], -0.15), "P4: w=(0,0,1,0) -> risk 0.09, return 0.15"); }
    { Problem p = P5(); Ind a = p.make({0.01, 0.05}); CHECK(near(a.f[0], 3.925, 1e-9) && near(a.f[1], 0.016, 1e-9) && a.cv == 0, "P5: (0.01,0.05) -> 3.925 kg, 0.016 m, feasible"); }
    { Problem p = P7(); Ind a = p.make({10, 80}), b = p.make({70, 20}); CHECK(near(a.f[0], -8) && near(a.f[1], -96) && near(b.f[0], -56) && near(b.f[1], -24) && a.cv == 0, "P7: line ends (8,96) and (56,24)"); }
    { Problem p = P11(); Ind a = p.make({300, 200}); CHECK(near(a.f[0], -2400) && near(a.f[1], -500) && a.cv == 0, "P11: (300,200) -> profit 2400, units 500"); }
    { Problem p = P14(); Ind a = p.make({13, 40}), b = p.make({13, 90}); CHECK(near(a.f[0], 1.53, 1e-9) && near(-a.f[1], 9.0, 0.05) && near(b.f[0], 2.13, 1e-9) && near(-b.f[1], 20.2, 0.05), "P14: (13,40)->(1.53kg,~9h), (13,90)->(2.13kg,~20.2h)"); }
    { Problem p = P10(5); Vec x(p.n, 0.5); for (int i = 0; i < p.M - 1; i++) x[i] = 0.3 + 0.1 * i; Ind a = p.make(x); double s = 0; for (double f : a.f) s += f*f; CHECK(near(s, 1.0, 1e-9), "P10: g=0 -> sum f_m^2 = 1"); }
    { Problem p = P3(); Vec x(10, 0.0); x[0] = 0.5; Ind a = p.make(x); CHECK(near(a.f[1], 0.75), "P3: on the front f2 = 1 - f1^2"); }
    { TSP t = TSP::make(15, 7); std::vector<int> id(15); std::iota(id.begin(), id.end(), 0); Vec f = t.eval(id); CHECK(f[0] > 0 && f[1] > 0 && t.d[2][5] == t.d[5][2], "P8: symmetric distance, positive tour cost"); }

    std::cout << "== per-algorithm unit tests (smoke + key rule) ==\n";
    Problem p1 = P1(), p3 = P3();
    for (auto& a : CATALOGUE) {
        if (a.id == 8 || a.id == 9 || a.id == 20 || a.id == 25 || a.id == 36 || a.id == 37) continue;
        seed(5); Problem P = (a.id == 6 || a.id == 7) ? P11() : (a.id == 15 || a.id == 31 || a.id == 34 ? P10(3) : p1);
        auto r = dispatch(a.id, P, 3000);
        bool inb = true; for (auto& s : r) for (int i = 0; i < P.n; i++) if (s.x[i] < P.lo[i] - 1e-9 || s.x[i] > P.hi[i] + 1e-9) inb = false;
        bool nd = true; if (a.id > 7) { for (size_t i = 0; i < r.size(); i++) for (size_t j = 0; j < r.size(); j++) if (i != j && dominates(r[i].f, r[j].f)) nd = false; }
        Met m = metrics(P, r);
        double need = (a.id <= 7) ? 0.0 : (P.M == 2 ? 0.45 : 0.1);
        CHECK(!r.empty() && inb && nd && m.hv >= need && P.nev <= 3200, std::string("alg ") + std::to_string(a.id) + " " + a.name + ": bounds ok, returned set non-dominated, HV=" + std::to_string(m.hv) + ", evals=" + std::to_string(P.nev));
    }
    // key-rule tests
    { Problem P = P11(); seed(3); Vec dev; Ind r = goal_programming(P, {-2400, -600}, {1, 1}, 6000, &dev);
      CHECK(dev[1] >= 99.0 && dev[1] <= 110.0 && dev[0] <= 15.0, "alg 6 GP on P11: unit shortfall >= 100 (goal 600 unreachable), profit goal ~met; d+=(" + std::to_string(dev[0]) + "," + std::to_string(dev[1]) + ")");
      auto lx = lexicographic(P, {0, 1}, 1e-3, 9000); CHECK(-lx[0].f[0] >= 2390 && -lx[1].f[1] >= 495 && -lx[1].f[0] >= 2390, "alg 7 Lexicographic: profit 2400 first, then units 500 at (300,200)"); }
    { std::vector<Vec> X = {{120,8,10,92},{100,7,14,88},{140,9,7,95},{90,6,18,80},{130,8,9,90}}; Vec w = {0.3,0.3,0.2,0.2}; std::vector<bool> ben = {false,true,false,true};
      auto t = topsis(X, w, ben); auto X2 = X; for (auto& r : X2) r[0] *= 1000; auto t2 = topsis(X2, w, ben);
      CHECK(t.order == t2.order, "alg 8 TOPSIS: ranking unchanged when one criterion is rescaled (vector normalisation)");
      auto v = vikor(X, w, ben, 1.0); auto v5 = vikor(X, w, ben, 0.5); bool in01 = true; for (double q : v5.score) if (q < -1e-12 || q > 1 + 1e-12) in01 = false;
      CHECK(in01 && near(*std::min_element(v.score.begin(), v.score.end()), 0.0, 1e-12) && near(*std::max_element(v.score.begin(), v.score.end()), 1.0, 1e-9), "alg 9 VIKOR: Q in [0,1] for v=0.5; with v=1 (pure group utility S) best Q = 0 and worst Q = 1"); }
    { Problem P = P13(); auto g = run_mgda(P, {0, 5}, 400); Vec x = g.x.back(); CHECK(std::fabs(x[0] - x[1]) < 0.05 && std::fabs(x[0]) <= 3.05, "alg 36 MGDA from (0,5) ends on the Pareto segment x1=x2, |x|<=3"); 
      auto gn = run_gradnorm(P, {0, 5}, 400); CHECK(gn.loss.back()[0] < gn.loss.front()[0] && gn.loss.back()[1] < gn.loss.front()[1], "alg 37 GradNorm: each task loss ends below its own start value"); }
    { ensureTSP(); seed(9); auto a = run_moeadaco(g_tsp, 2000); auto b = run_moaco(g_tsp, 2000);
      bool okA = !a.empty(), okB = !b.empty(); for (auto& t : a) { auto q = t.p; std::sort(q.begin(), q.end()); for (int i = 0; i < 15; i++) if (q[i] != i) okA = false; }
      for (auto& t : b) { auto q = t.p; std::sort(q.begin(), q.end()); for (int i = 0; i < 15; i++) if (q[i] != i) okB = false; }
      CHECK(okA, "alg 20 MOEA/D-ACO returns valid permutations (every city exactly once)"); CHECK(okB, "alg 25 MOACO returns valid permutations (every city exactly once)"); }
    { Problem P = P6(); seed(2); g_deb = true; auto r = run_nsga2(P, 4000); Met m = metrics(P, r); CHECK(m.feas >= 0.999, "alg 21 Constrained NSGA-II on P6: all returned points feasible"); }
    std::cout << (g_fail == 0 ? "ALL TESTS PASSED\n" : "SOME TESTS FAILED: " + std::to_string(g_fail) + "\n");
}

// ---------------- Task A3: implementation log ----------------
static void cmdLog(long budget, uint64_t sd0) {
    fs::create_directories("results/fronts"); std::ofstream log("results/implementation_log.csv");
    log << "number,algorithm,problem,budget,evals_used,n_points,HV,IGD,Spacing,feasible_fraction,runtime_ms,status\n";
    std::cout << "== Task A3 implementation log (budget " << budget << ", seed " << sd0 << ") ==\n";
    for (auto& a : CATALOGUE) {
        if (a.id == 8 || a.id == 9 || a.id == 36 || a.id == 37) continue;   // handled by cmdMcdm / cmdMtl
        for (auto& pk : a.probs) {
            if (pk == "P12") continue;
            seed(sd0 + a.id); auto t0 = std::chrono::steady_clock::now(); Met m; long ev = 0; std::string status = "works";
            try {
                if (pk == "P8") { ensureTSP(); seed(sd0 + a.id); t0 = std::chrono::steady_clock::now(); auto r = runTSP(a.id, budget); m = tspMetrics(r); ev = budget;
                    std::ofstream o("results/fronts/alg" + std::to_string(a.id) + "_P8.csv"); for (auto& t : r) o << t.f[0] << "," << t.f[1] << "\n"; }
                else { Problem P = makeProblem(pk); writeTrue(P); auto r = dispatch(a.id, P, budget); ev = P.nev; m = metrics(P, r);
                    writeFront("results/fronts/alg" + std::to_string(a.id) + "_" + pk + ".csv", P, r); }
            } catch (std::exception& e) { status = std::string("not done: ") + e.what(); }
            double t = ms(t0);
            if (a.id == 14 && (pk == "P5" || pk == "P6") && m.feas < 0.99) status = "works (no feasibility rules: infeasible points returned)";
            if (a.id <= 7) status += " (sweep/preference method: HV of a few points)";
            char buf[400]; snprintf(buf, sizeof buf, "%d,%s,%s,%ld,%ld,%d,%.5f,%.5f,%.5f,%.3f,%.1f,%s\n", a.id, a.name.c_str(), pk.c_str(), budget, ev, m.n, m.hv, m.igd, m.sp, m.feas, t, status.c_str());
            log << buf; std::printf("%2d %-28s %-6s n=%3d HV=%.4f IGD=%.4f SP=%.4f feas=%.2f %7.0f ms\n", a.id, a.name.c_str(), pk.c_str(), m.n, m.hv, m.igd, m.sp, m.feas, t);
        }
    }
    std::cout << "(algorithms 8, 9 -> './moo mcdm'; 36, 37 -> './moo mtl')\n";
}

// ---------------- Task D: head-to-head ----------------
struct Row { int id; std::string name; double hv, hvs, sp, sps, igd, igds; };
static void printRank(const std::string& title, std::vector<Row>& rows, std::ofstream& csv, const std::string& tag) {
    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.hv > b.hv; });
    std::cout << "\n--- " << title << " (ranked by mean HV) ---\n rank id algorithm                    HV (mean±sd)        Spacing (mean±sd)   IGD (mean±sd)\n";
    int r = 1; for (auto& w : rows) { std::printf("  %2d  %2d %-26s %.4f ± %.4f    %.4f ± %.4f    %.4f ± %.4f\n", r, w.id, w.name.c_str(), w.hv, w.hvs, w.sp, w.sps, w.igd, w.igds);
        csv << tag << "," << r << "," << w.id << "," << w.name << "," << w.hv << "," << w.hvs << "," << w.sp << "," << w.sps << "," << w.igd << "," << w.igds << "\n"; r++; }
}
static void cmdCompare(int seeds, long budget, uint64_t sd0) {
    fs::create_directories("results"); std::ofstream csv("results/taskD_tables.csv"); csv << "table,rank,id,algorithm,HV_mean,HV_sd,Spacing_mean,Spacing_sd,IGD_mean,IGD_sd\n";
    std::cout << "== Task D: " << seeds << " seeds, budget " << budget << " evaluations ==\n";
    // reference fronts for problems without analytic front: pooled non-dominated union of all runs
    for (std::string pk : {"P3", "P5"}) {
        Problem P = makeProblem(pk); std::vector<int> ids = {1,2,3,4,10,11,12,13,14,15,16,17,18,19,21,22,23,24,26,27,28,29,30,31,32,33,34,35};
        std::map<int, std::vector<std::vector<Ind>>> runs; std::vector<Vec> pool;
        for (int id : ids) for (int s = 0; s < seeds; s++) { Problem Q = makeProblem(pk); seed(sd0 + 1000 * id + s); auto r = dispatch(id, Q, budget); runs[id].push_back(r); }
        if (P.trueFront.empty()) { for (auto& kv : runs) for (auto& r : kv.second) for (auto& a : r) if (a.cv <= 1e-9) pool.push_back(P.norm(a.f)); pool = nondominatedSubset(pool); for (auto& v : pool) for (int i = 0; i < P.M; i++) v[i] = v[i] * (P.fmax[i] - P.fmin[i]) + P.fmin[i]; }
        std::vector<Row> rows;
        for (int id : ids) { std::vector<double> H, S, I; for (auto& r : runs[id]) { Met m = metrics(P, r, P.trueFront.empty() ? &pool : nullptr); H.push_back(m.hv); S.push_back(m.sp); I.push_back(m.igd); }
            std::string nm; for (auto& a : CATALOGUE) if (a.id == id) nm = a.name; rows.push_back({id, nm, mean(H), sd(H), mean(S), sd(S), mean(I), sd(I)}); }
        printRank("D1 continuous bi-objective on " + pk + (P.trueFront.empty() ? " (IGD vs pooled reference front)" : " (IGD vs true front)"), rows, csv, "D1_" + pk);
    }
    { // D2: many objective P10, M = 5
        std::vector<int> ids = {14, 15, 18, 31, 34}; std::vector<Row> rows;
        for (int id : ids) { std::vector<double> H, S, I; for (int s = 0; s < seeds; s++) { Problem Q = makeProblem("P10_5"); seed(sd0 + 7000 + 100 * id + s); auto r = dispatch(id, Q, budget); Met m = metrics(Q, r); H.push_back(m.hv); S.push_back(m.sp); I.push_back(m.igd); }
            std::string nm; for (auto& a : CATALOGUE) if (a.id == id) nm = a.name; rows.push_back({id, nm, mean(H), sd(H), mean(S), sd(S), mean(I), sd(I)}); }
        printRank("D2 many-objective DTLZ2 M=5 (HV by Monte-Carlo, IGD vs true front)", rows, csv, "D2_P10_M5");
    }
    { // D3: discrete P8
        ensureTSP(); std::vector<int> ids = {20, 25, 14, 30}; std::vector<Row> rows; const char* nm[] = {"MOEA/D-ACO", "MOACO", "NSGA-II (perm)", "MOMA (perm)"};
        for (int q = 0; q < 4; q++) { std::vector<double> H, S, I; for (int s = 0; s < seeds; s++) { seed(sd0 + 9000 + 100 * ids[q] + s); Met m = tspMetrics(runTSP(ids[q], budget)); H.push_back(m.hv); S.push_back(m.sp); I.push_back(m.igd); }
            rows.push_back({ids[q], nm[q], mean(H), sd(H), mean(S), sd(S), mean(I), sd(I)}); }
        printRank("D3 discrete bi-objective TSP (P8, 15 cities)", rows, csv, "D3_P8");
        std::vector<int> id1 = {0}; std::vector<int> nnD = nearestNeighbour(g_tsp, g_tsp.d), nnT = nearestNeighbour(g_tsp, g_tsp.t);
        Vec fd = g_tsp.eval(nnD), ft = g_tsp.eval(nnT); std::printf("  nearest-neighbour baselines: distance-NN tour (%.1f, %.1f), time-NN tour (%.1f, %.1f)\n", fd[0], fd[1], ft[0], ft[1]);
    }
}

// ---------------- Task D (algs 8, 9): MCDM ----------------
static void cmdMcdm() {
    std::cout << "== P12 supplier selection (criteria: cost min, quality max, delivery min, reliability max) ==\n";
    std::vector<Vec> X = {{120,8,10,92},{100,7,14,88},{140,9,7,95},{90,6,18,80},{130,8,9,90}}; Vec w = {0.3,0.3,0.2,0.2}; std::vector<bool> ben = {false,true,false,true};
    std::vector<Vec> G; for (auto& r : X) G.push_back({r[0], -r[1], r[2], -r[3]});
    std::cout << "Pareto-optimal suppliers:"; for (int i = 0; i < 5; i++) { bool dom = false; for (int j = 0; j < 5; j++) if (i != j && dominates(G[j], G[i])) dom = true; if (!dom) std::cout << " S" << i + 1; } std::cout << "\n";
    auto show = [&](const char* tag, const Ranking& r) { std::cout << "  " << tag << " ranking:"; for (int i : r.order) std::cout << " S" << i + 1 << "(" << r.score[i] << ")"; std::cout << "\n"; };
    show("TOPSIS", topsis(X, w, ben, true));
    for (double v : {0.25, 0.5, 0.75}) { std::cout << " VIKOR v=" << v << "\n"; show("VIKOR", vikor(X, w, ben, v, true)); }
    std::cout << "-- weight shift of 0.1 (w1: 0.3->0.4, w2: 0.3->0.2) --\n"; Vec w2 = {0.4,0.2,0.2,0.2};
    show("TOPSIS", topsis(X, w2, ben)); for (double v : {0.25, 0.5, 0.75}) show(("VIKOR v=" + std::to_string(v).substr(0,4)).c_str(), vikor(X, w2, ben, v));
    std::cout << "\n== TOPSIS / VIKOR applied to the NSGA-II front of P4 (risk min, return max), weights (0.5,0.5) ==\n";
    Problem P = P4(); seed(11); auto fr = dispatch(14, P, 10000); std::vector<Vec> Y; std::vector<Ind> kept; for (auto& a : fr) { Y.push_back({a.f[0], -a.f[1]}); kept.push_back(a); }
    Vec wp = {0.5, 0.5}; std::vector<bool> bp = {false, true};
    auto top3 = [&](const char* tag, const Ranking& r) { std::cout << "  " << tag << " top 3 (risk, return):"; for (int k = 0; k < 3; k++) { int i = r.order[k]; std::cout << " [#" << i << " (" << Y[i][0] << "," << Y[i][1] << ")]"; } std::cout << "\n"; };
    auto r1 = topsis(Y, wp, bp); top3("TOPSIS", r1); auto r2 = vikor(Y, wp, bp, 0.5); top3("VIKOR v=0.5", r2);
    Vec wq = {0.6, 0.4}; auto r3 = topsis(Y, wq, bp); auto r4 = vikor(Y, wq, bp, 0.5); top3("TOPSIS (w=0.6/0.4)", r3); top3("VIKOR  (w=0.6/0.4)", r4);
    std::cout << "  TOPSIS best changed after the 0.1 shift: " << (r1.order[0] != r3.order[0] ? "YES" : "no") << ";  VIKOR best changed: " << (r2.order[0] != r4.order[0] ? "YES" : "no") << "\n";
}
// ---------------- P11: goal programming & lexicographic ----------------
static void cmdP11() {
    Problem P = P11(); std::cout << "== P11: Goal Programming and Lexicographic ==\n"; seed(21);
    Vec t = {-2400, -600}, dev;
    Ind a = goal_programming(P, t, {1, 1}, 12000, &dev); std::printf("GP unit weights      : xA=%.1f xB=%.1f profit=%.1f units=%.1f  d+ profit=%.2f d+ units=%.2f\n", a.x[0], a.x[1], -a.f[0], -a.f[1], dev[0], dev[1]);
    Ind b = goal_programming(P, t, {1.0/2400, 1.0/600}, 12000, &dev); std::printf("GP normalised weights: xA=%.1f xB=%.1f profit=%.1f units=%.1f  d+ profit=%.2f d+ units=%.2f\n", b.x[0], b.x[1], -b.f[0], -b.f[1], dev[0], dev[1]);
    auto lx = lexicographic(P, {0, 1}, 1e-3, 12000);
    std::printf("Lexicographic step 1 (profit only): xA=%.1f xB=%.1f profit=%.1f units=%.1f  <- any point on the edge (300,200)-(0,400) is optimal\n", lx[0].x[0], lx[0].x[1], -lx[0].f[0], -lx[0].f[1]);
    std::printf("Lexicographic step 2 (units | profit kept): xA=%.1f xB=%.1f profit=%.1f units=%.1f\n", lx[1].x[0], lx[1].x[1], -lx[1].f[0], -lx[1].f[1]);
}
// ---------------- P13: multi-task toy ----------------
static void cmdMtl() {
    fs::create_directories("results"); Problem P = P13(); std::cout << "== P13(a) toy multi-task problem, start (0,5), 300 steps, step size 0.05 ==\n";
    auto m = run_mgda(P, {0, 5}, 300), g = run_gradnorm(P, {0, 5}, 300), f = run_fixed(P, {0, 5}, 300);
    std::ofstream o("results/P13_trajectories.csv"); o << "iter,fixed_L1,fixed_L2,mgda_L1,mgda_L2,gradnorm_L1,gradnorm_L2,mgda_a1,gradnorm_w1,gradnorm_w2,mgda_x1,mgda_x2\n";
    for (int t = 0; t < 300; t++) o << t << "," << f.loss[t][0] << "," << f.loss[t][1] << "," << m.loss[t][0] << "," << m.loss[t][1] << "," << g.loss[t][0] << "," << g.loss[t][1] << "," << m.w[t][0] << "," << g.w[t][0] << "," << g.w[t][1] << "," << m.x[t][0] << "," << m.x[t][1] << "\n";
    auto line = [&](const char* n, const GDTrace& r) { std::printf("  %-9s start L=(%.2f,%.2f)  step 50 L=(%.3f,%.3f)  final L=(%.4f,%.4f)  final x=(%.4f,%.4f)\n", n, r.loss[0][0], r.loss[0][1], r.loss[50][0], r.loss[50][1], r.loss.back()[0], r.loss.back()[1], r.x.back()[0], r.x.back()[1]); };
    line("fixed 1/2", f); line("MGDA", m); line("GradNorm", g);
    std::cout << "  Pareto set is x1=x2 with -3<=x1<=3: MGDA end point on it: " << (std::fabs(m.x.back()[0] - m.x.back()[1]) < 0.05 ? "yes" : "no") << "\n";
    // one task with much larger gradients: scale task 2 by 50
    Problem Q = P13(); Q.grad = [&](const Vec& x, int k) { Vec g0 = P.grad(x, k); if (k == 1) for (double& v : g0) v *= 50; return g0; };
    Q.fn = [&](const Vec& x, Vec& ff, double& cv) { P.fn(x, ff, cv); ff[1] *= 50; };
    auto f2 = run_fixed(Q, {0, 5}, 300, 0.001), m2 = run_mgda(Q, {0, 5}, 300, 0.001), g2 = run_gradnorm(Q, {0, 5}, 300, 0.001);
    std::cout << "  -- task 2 scaled x50 (large gradients), step 0.001 --\n"; line("fixed 1/2", f2); line("MGDA", m2); line("GradNorm", g2);
    std::ofstream o2("results/P13_front_mgda.csv"); for (double sx = -9; sx <= 9; sx += 3) for (double sy : {-9.0, 0.0, 5.0, 9.0}) { auto r = run_mgda(P, {sx, sy}, 600); o2 << r.loss.back()[0] << "," << r.loss.back()[1] << "\n"; }
}
// ---------------- investigations from Section 4 ----------------
static void cmdInvestigate() {
    std::cout << "== investigations (evidence for Task B/E) ==\n"; fs::create_directories("results");
    { Problem P = P1(); seed(1); auto ws = run_weighted_sum(P, 4000, 20); auto key = [&](std::vector<Ind>& r) { std::set<long> s; for (auto& a : r) s.insert(std::lround(a.f[0] * 100)); return s.size(); };
      std::cout << "P1 weighted sum, 20 equally spaced weights -> distinct points (0.01 resolution): " << key(ws) << "; f1 values:"; for (auto& a : ws) std::printf(" %.2f", a.f[0]); std::cout << "\n";
      for (double p : {1.0, 2.0, 0.0}) { seed(1); auto r = run_weight_metric(P, 4000, p, 20); std::printf("P1 weight metric p=%s -> distinct points %zu, f1 range [%.2f, %.2f]\n", p == 0 ? "inf" : (p == 1 ? "1" : "2"), key(r), std::min_element(r.begin(), r.end(), [](const Ind& a, const Ind& b) { return a.f[0] < b.f[0]; })->f[0], std::max_element(r.begin(), r.end(), [](const Ind& a, const Ind& b) { return a.f[0] < b.f[0]; })->f[0]); } }
    { Problem P = P3(); seed(1); auto ws = run_weighted_sum(P, 6000, 20); int mid = 0; for (auto& a : ws) if (a.f[0] > 0.15 && a.f[0] < 0.85) mid++;
      seed(1); auto tc = run_weight_metric(P, 6000, 0, 20); int midT = 0; for (auto& a : tc) if (a.f[0] > 0.15 && a.f[0] < 0.85) midT++;
      std::cout << "P3 (ZDT2, non-convex): weighted-sum points with 0.15<f1<0.85: " << mid << " of 20 (expected ~0);  Tchebycheff: " << midT << " of 20\n"; }
    { Problem P = P2(); seed(1); auto r = dispatch(14, P, 5000); Problem Q = P2(true); seed(1); auto ri = dispatch(14, Q, 5000); std::set<long> a, b; for (auto& s : r) a.insert(std::lround(s.f[0] * 10)); for (auto& s : ri) b.insert(std::lround(s.f[0] * 10));
      double lo = 1e9, hi = -1e9; for (auto& s : r) { lo = std::min(lo, s.x[0]); hi = std::max(hi, s.x[0]); }
      std::cout << "P2 NSGA-II real x: " << r.size() << " front points, x in [" << lo << "," << hi << "]; integer x: " << b.size() << " distinct front points (max possible 10)\n"; }
    { Problem P = P5(); seed(1); int inf = 0, T = 100000; for (int i = 0; i < T; i++) { Ind a = randInd(P); if (a.cv > 0) inf++; } std::cout << "P5: fraction of random designs violating the stress limit = " << (double)inf / T << "\n";
      g_deb = true; seed(2); auto r1 = run_nsga2(P, 10000); g_deb = false; seed(2); auto r2 = run_nsga2(P, 10000); g_deb = true; Met m1 = metrics(P, r1), m2 = metrics(P, r2);
      std::printf("P5 NSGA-II with Deb rules: feasible fraction %.2f, HV %.4f;  ignoring constraints: feasible fraction %.2f, HV(feasible only) %.4f\n", m1.feas, m1.hv, m2.feas, m2.hv);
      double mn = 1e9; for (auto& a : r2) mn = std::min(mn, a.f[0]); double mn1 = 1e9; for (auto& a : r1) mn1 = std::min(mn1, a.f[0]); std::printf("   lightest design: with rule %.2f kg, ignoring stress limit %.2f kg\n", mn1, mn); }
    { Problem P = P6(); g_deb = false; seed(4); auto r = run_nsga2(P, 10000); Met m = metrics(P, r); g_deb = true; std::printf("P6 NSGA-II without feasibility rules: %.0f%% of returned points are feasible\n", 100 * m.feas); }
    for (int M : {3, 5}) { Problem P = P10(M); int tot = 0, nd = 0; for (int rep = 0; rep < 50; rep++) { seed(100 + rep); std::vector<Vec> F; for (int i = 0; i < 100; i++) F.push_back(randInd(P).f); nd += (int)nondominatedSubset(F).size(); tot += 100; }
      std::printf("P10 M=%d: fraction of a random population of 100 that is non-dominated = %.3f\n", M, (double)nd / tot); }
    { Problem P = P14(); seed(3); auto r = dispatch(14, P, 5000); double mx = 0, mn = 1e9; for (auto& a : r) { mx = std::max(mx, a.x[0]); mn = std::min(mn, a.x[0]); }
      std::printf("P14 NSGA-II: x1 over the front in [%.3f, %.3f] (theory: all 13); reference (1.8kg,12h) is behind the front: best runtime at 1.8 kg is %.2f h\n", mn, mx, 0.2 * ((1.8 - 1.05) / 0.012) / 0.89);
      for (int c = 0; c < 2; c++) { Problem Q = P14(); if (c == 1) Q.prefRef = {(1.8 - 1.53) / 0.6, (-12.0 + 20.2) / 11.2}; else Q.prefRef = {0.0, 0.0}; seed(3); auto rr = dispatch(35, Q, 5000);
        double lo = 1e9, hi = -1e9; for (auto& a : rr) { lo = std::min(lo, a.f[0]); hi = std::max(hi, a.f[0]); } std::printf("   R-NSGA-II with %s reference: weight range of returned points [%.2f, %.2f] kg\n", c ? "(1.8 kg, 12 h)" : "ideal-corner", lo, hi); } }
    { Problem P = P7(); Vec pr = {(-40.0 + 56) / 48, (-50.0 + 96) / 72}; P.prefRef = pr; seed(5); auto r = dispatch(35, P, 5000); double a1 = 0, a2 = 0; for (auto& s : r) { a1 += -s.f[0]; a2 += -s.f[1]; } std::printf("P7 R-NSGA-II with reference (40,50): mean returned point (agri yield %.1f, industry output %.1f); on-line check max |x1+x2-90| = ", a1 / r.size(), a2 / r.size());
      double mx = 0; for (auto& s : r) mx = std::max(mx, std::fabs(s.x[0] + s.x[1] - 90)); std::printf("%.4f\n", mx); }
}
// ---------------- Task C: hand-trace helpers & sensitivity ----------------
static void cmdTrace() {
    std::cout << "================ TOPSIS (alg 8) on P12 : every number ================\n";
    std::vector<Vec> X = {{120,8,10,92},{100,7,14,88},{140,9,7,95},{90,6,18,80},{130,8,9,90}}; topsis(X, {0.3,0.3,0.2,0.2}, {false,true,false,true}, true);
    std::cout << "\n================ IBEA (alg 32): tiny instance, population 4 on P1, two iterations, kappa=0.05 ================\n";
    Problem P = P1(); seed(42); g_deb = true; auto pop = initPop(P, 4); const double kappa = 0.05;
    for (int it = 1; it <= 2; it++) {
        std::cout << " iteration " << it << " population (x -> normalised f):\n"; for (auto& a : pop) { Vec v = P.norm(a.f); std::printf("   x=%.4f  f=(%.4f,%.4f)  fn=(%.4f,%.4f)\n", a.x[0], a.f[0], a.f[1], v[0], v[1]); }
        std::vector<Vec> F; for (auto& a : pop) F.push_back(P.norm(a.f)); std::vector<double> fit; std::vector<std::vector<double>> I; double c; ibeaFitness(F, kappa, fit, I, c);
        std::cout << "  I_eps(i,j) = max_k (fn_i,k - fn_j,k):\n"; for (auto& r : I) { std::cout << "   "; for (double v : r) std::printf(" %8.4f", v); std::cout << "\n"; }
        std::cout << "  fitness F(x)=sum_j -exp(-I(j,x)/(c*kappa)):"; for (double f : fit) std::printf(" %.4f", f); std::cout << "\n";
        std::vector<Ind> off; for (int k = 0; k < 4; k++) { int a = irand(4), b = irand(4); Vec c1 = pop[fit[a] >= fit[b] ? a : b].x, c2 = pop[irand(4)].x; sbx(c1, c2, P); polyMut(c1, P); off.push_back(P.make(c1)); }
        std::cout << "  offspring x:"; for (auto& a : off) std::printf(" %.4f", a.x[0]); std::cout << "\n  environmental selection (8 -> 4):\n";
        auto R = pop; R.insert(R.end(), off.begin(), off.end()); ibeaReduce(P, R, 4, kappa, true); pop = R;
    }
    std::cout << "\n================ MOEA/D-ACO (alg 20): 5-city instance, 2 subproblems, 2 iterations ================\n";
    std::vector<Vec> d = {{0,3,4,2,5},{3,0,2,6,3},{4,2,0,3,4},{2,6,3,0,2},{5,3,4,2,0}};
    std::vector<Vec> t = {{0,6,4,3,5},{6,0,4,3,6},{4,4,0,6,2},{3,3,6,0,4},{5,6,2,4,0}};
    TSP s = TSP::fromMatrices(d, t); s.fmin = {10, 10}; s.fmax = {25, 30}; ACOParams pr; pr.N = 2; pr.nbr = 2; pr.rho = 0.1; seed(7);
    std::cout << " distance matrix d and time matrix t are symmetric (given in the code); alpha=1 beta=2 rho=0.1, 2 neighbours, normalisation fmin=(10,10) fmax=(25,30)\n";
    run_moeadaco(s, 2 * 2 + 2, pr, true);   // 2 initial tours + 2 iterations x 2 ants
}
static void cmdSens(int seeds, long budget, uint64_t sd0) {
    fs::create_directories("results"); std::ofstream csv("results/taskC_sensitivity.csv"); csv << "algorithm,parameter,value,HV_mean,HV_sd,IGD_mean,IGD_sd\n";
    std::cout << "== Task C sensitivity (" << seeds << " seeds, budget " << budget << ") ==\n";
    { Problem P = P3(); for (double kappa : {0.01, 0.05, 0.25}) { std::vector<double> H, I; for (int s = 0; s < seeds; s++) { Problem Q = P3(); seed(sd0 + s); g_deb = true; auto r = run_ibea(Q, budget, kappa); Met m = metrics(Q, r); H.push_back(m.hv); I.push_back(m.igd); }
        std::printf("IBEA kappa=%.2f on P3: HV %.4f ± %.4f, IGD %.4f ± %.4f\n", kappa, mean(H), sd(H), mean(I), sd(I)); csv << "IBEA,kappa," << kappa << "," << mean(H) << "," << sd(H) << "," << mean(I) << "," << sd(I) << "\n"; } }
    { ensureTSP(); for (double rho : {0.02, 0.1, 0.5}) { std::vector<double> H, I; for (int s = 0; s < seeds; s++) { seed(sd0 + s); ACOParams pr; pr.rho = rho; Met m = tspMetrics(run_moeadaco(g_tsp, budget, pr)); H.push_back(m.hv); I.push_back(m.igd); }
        std::printf("MOEA/D-ACO rho=%.2f on P8: HV %.4f ± %.4f, IGD %.4f ± %.4f\n", rho, mean(H), sd(H), mean(I), sd(I)); csv << "MOEA/D-ACO,rho," << rho << "," << mean(H) << "," << sd(H) << "," << mean(I) << "," << sd(I) << "\n"; } }
    { std::vector<Vec> X = {{120,8,10,92},{100,7,14,88},{140,9,7,95},{90,6,18,80},{130,8,9,90}}; std::vector<bool> ben = {false,true,false,true};
      std::cout << "TOPSIS: most important parameter = the weight vector. Shifting weight 0.1 from criterion 2 to criterion 1 (three settings):\n";
      for (double w1 : {0.3, 0.4, 0.5}) { Vec w = {w1, 0.6 - w1, 0.2, 0.2}; auto r = topsis(X, w, ben); std::cout << "  w=(" << w[0] << "," << w[1] << ",0.2,0.2): order"; for (int i : r.order) std::cout << " S" << i + 1; std::cout << "\n"; } }
}
int main(int argc, char** argv) {
    std::string cmd = argc > 1 ? argv[1] : "help"; long budget = 10000; int seeds = 10;
    for (int i = 2; i + 1 < argc; i++) { if (std::string(argv[i]) == "--budget") budget = std::atol(argv[i + 1]); if (std::string(argv[i]) == "--seeds") seeds = std::atoi(argv[i + 1]); }
    fs::create_directories("results");
    if (cmd == "unit") { unitTests(); return g_fail ? 1 : 0; }
    if (cmd == "log") { cmdLog(budget, 2026); return 0; }
    if (cmd == "compare") { cmdCompare(seeds, budget, 5000); return 0; }
    if (cmd == "mcdm") { cmdMcdm(); return 0; }
    if (cmd == "p11") { cmdP11(); return 0; }
    if (cmd == "mtl") { cmdMtl(); return 0; }
    if (cmd == "investigate") { cmdInvestigate(); return 0; }
    if (cmd == "trace") { cmdTrace(); return 0; }
    if (cmd == "sens") { cmdSens(seeds, budget, 777); return 0; }
    if (cmd == "all") { unitTests(); cmdLog(budget, 2026); cmdP11(); cmdMcdm(); cmdMtl(); cmdInvestigate(); cmdTrace(); cmdSens(seeds, budget, 777); cmdCompare(seeds, budget, 5000); return 0; }
    std::cout << "usage: ./moo <unit|log|compare|mcdm|p11|mtl|investigate|trace|sens|all> [--budget N] [--seeds K]\n"; return 0;
}
