#include "moo_utils.hpp"
namespace moo {
std::mt19937_64 rng(1);
bool g_deb = true;
void seed(uint64_t s) { rng.seed(s); }
double urand() { return std::uniform_real_distribution<double>(0.0, 1.0)(rng); }
double urand(double a, double b) { return a + (b - a) * urand(); }
int irand(int n) { return (int)(rng() % (uint64_t)n); }
double nrand() { return std::normal_distribution<double>(0.0, 1.0)(rng); }

bool dominates(const Vec& a, const Vec& b) {
    bool strict = false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] > b[i]) return false;
        if (a[i] < b[i]) strict = true;
    }
    return strict;
}

// Deb's rules: (1) feasible beats infeasible, (2) both infeasible: smaller violation wins,
// (3) both feasible: ordinary Pareto dominance.
bool cdom(const Ind& a, const Ind& b) {
    if (!g_deb) return dominates(a.f, b.f);
    bool fa = a.cv <= 1e-12, fb = b.cv <= 1e-12;
    if (fa && !fb) return true;
    if (!fa && fb) return false;
    if (!fa && !fb) return a.cv < b.cv;
    return dominates(a.f, b.f);
}

std::vector<std::vector<int>> ndsGeneric(int n, const std::function<bool(int, int)>& dom) {
    std::vector<std::vector<int>> S(n), fronts;
    std::vector<int> cnt(n, 0);
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            if (dom(i, j)) { S[i].push_back(j); cnt[j]++; }
            else if (dom(j, i)) { S[j].push_back(i); cnt[i]++; }
        }
    std::vector<int> cur;
    for (int i = 0; i < n; i++) if (cnt[i] == 0) cur.push_back(i);
    while (!cur.empty()) {
        fronts.push_back(cur);
        std::vector<int> nxt;
        for (int i : cur) for (int j : S[i]) if (--cnt[j] == 0) nxt.push_back(j);
        cur = nxt;
    }
    return fronts;
}
std::vector<std::vector<int>> fastNonDominatedSort(const std::vector<Ind>& P) {
    return ndsGeneric((int)P.size(), [&](int i, int j) { return cdom(P[i], P[j]); });
}
std::vector<std::vector<int>> fastNonDominatedSortF(const std::vector<Vec>& F) {
    return ndsGeneric((int)F.size(), [&](int i, int j) { return dominates(F[i], F[j]); });
}
std::vector<Vec> nondominatedSubset(const std::vector<Vec>& F) {
    std::vector<Vec> out;
    if (F.empty()) return out;
    auto fr = fastNonDominatedSortF(F);
    for (int i : fr[0]) out.push_back(F[i]);
    return out;
}

std::vector<double> crowdingDistance(const std::vector<Vec>& F) {
    int n = (int)F.size();
    const double INF = 1e18;
    std::vector<double> d(n, 0.0);
    if (n == 0) return d;
    if (n <= 2) { std::fill(d.begin(), d.end(), INF); return d; }
    int M = (int)F[0].size();
    std::vector<int> idx(n);
    for (int m = 0; m < M; m++) {
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) { return F[a][m] < F[b][m]; });
        double lo = F[idx[0]][m], hi = F[idx[n - 1]][m];
        d[idx[0]] = d[idx[n - 1]] = INF;
        if (hi - lo < 1e-14) continue;
        for (int k = 1; k < n - 1; k++) d[idx[k]] += (F[idx[k + 1]][m] - F[idx[k - 1]][m]) / (hi - lo);
    }
    return d;
}

double hv2d(const std::vector<Vec>& pts, const Vec& ref) {
    std::vector<Vec> p;
    for (auto& v : pts) if (v[0] < ref[0] && v[1] < ref[1]) p.push_back(v);
    std::sort(p.begin(), p.end(), [](const Vec& a, const Vec& b) { return a[0] < b[0]; });
    double hv = 0, prev = ref[1];
    for (auto& v : p) if (v[1] < prev) { hv += (ref[0] - v[0]) * (prev - v[1]); prev = v[1]; }
    return hv;
}
double hvMC(const std::vector<Vec>& pts, const Vec& ref, int samples, uint64_t sd) {
    int M = (int)ref.size();
    std::vector<Vec> p;
    for (auto& v : pts) { bool ok = true; for (int i = 0; i < M; i++) if (v[i] >= ref[i]) ok = false; if (ok) p.push_back(v); }
    if (p.empty()) return 0.0;
    Vec lo(M, 1e18);
    for (auto& v : p) for (int i = 0; i < M; i++) lo[i] = std::min(lo[i], v[i]);
    double vol = 1; for (int i = 0; i < M; i++) vol *= (ref[i] - lo[i]);
    std::mt19937_64 g(sd);
    std::uniform_real_distribution<double> U(0, 1);
    int hit = 0; Vec s(M);
    for (int k = 0; k < samples; k++) {
        for (int i = 0; i < M; i++) s[i] = lo[i] + U(g) * (ref[i] - lo[i]);
        for (auto& v : p) { bool d = true; for (int i = 0; i < M; i++) if (v[i] > s[i]) { d = false; break; } if (d) { hit++; break; } }
    }
    return vol * hit / samples;
}
double hypervolume(const std::vector<Vec>& pts, const Vec& ref) {
    return ref.size() == 2 ? hv2d(pts, ref) : hvMC(pts, ref);
}
double igd(const std::vector<Vec>& T, const std::vector<Vec>& A) {
    if (A.empty() || T.empty()) return std::numeric_limits<double>::infinity();
    double s = 0;
    for (auto& t : T) {
        double best = 1e18;
        for (auto& a : A) { double d = 0; for (size_t i = 0; i < t.size(); i++) d += (t[i] - a[i]) * (t[i] - a[i]); best = std::min(best, d); }
        s += std::sqrt(best);
    }
    return s / T.size();
}
double spacing(const std::vector<Vec>& A) {
    int n = (int)A.size();
    if (n < 2) return 0.0;
    std::vector<double> d(n, 1e18);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) if (i != j) { double s = 0; for (size_t k = 0; k < A[i].size(); k++) s += std::fabs(A[i][k] - A[j][k]); d[i] = std::min(d[i], s); }
    double mean = 0; for (double v : d) mean += v; mean /= n;
    double var = 0; for (double v : d) var += (v - mean) * (v - mean);
    return std::sqrt(var / (n - 1));
}

static void ddRec(Vec& cur, int idx, int left, int M, int p, std::vector<Vec>& out) {
    if (idx == M - 1) { cur[idx] = (double)left / p; out.push_back(cur); return; }
    for (int i = 0; i <= left; i++) { cur[idx] = (double)i / p; ddRec(cur, idx + 1, left - i, M, p, out); }
}
std::vector<Vec> dasDennis(int M, int p) {
    std::vector<Vec> out; Vec cur(M, 0.0);
    ddRec(cur, 0, p, M, p, out);
    return out;
}
bool solveLinear(std::vector<Vec> A, Vec b, Vec& x) {
    int n = (int)b.size();
    for (int c = 0; c < n; c++) {
        int piv = c;
        for (int r = c + 1; r < n; r++) if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
        if (std::fabs(A[piv][c]) < 1e-12) return false;
        std::swap(A[piv], A[c]); std::swap(b[piv], b[c]);
        for (int r = c + 1; r < n; r++) {
            double f = A[r][c] / A[c][c];
            for (int k = c; k < n; k++) A[r][k] -= f * A[c][k];
            b[r] -= f * b[c];
        }
    }
    x.assign(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double s = b[i];
        for (int k = i + 1; k < n; k++) s -= A[i][k] * x[k];
        x[i] = s / A[i][i];
    }
    return true;
}
}  // namespace moo
