// moo_utils.hpp - shared toolkit (Task A1): RNG, dominance, Deb's constrained dominance,
// fast non-dominated sorting, crowding distance, hypervolume (exact 2-D, Monte-Carlo M>=3),
// IGD, Spacing, Das-Dennis reference points, small linear solver.
// Everything here is hand-written with the C++17 standard library only.
#pragma once
#include <vector>
#include <random>
#include <functional>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <numeric>
#include <iostream>
#include <cstdint>
#include <stdexcept>

namespace moo {
using Vec = std::vector<double>;

extern std::mt19937_64 rng;
void seed(uint64_t s);
double urand();                    // U[0,1)
double urand(double a, double b);  // U[a,b)
int irand(int n);                  // integer in [0,n)
double nrand();                    // N(0,1)

extern bool g_deb;  // true: use Deb's feasibility rules / penalties; false: ignore constraints

// one candidate solution. All objectives are MINIMISED (maximisation is negated in problems.hpp)
struct Ind {
    Vec x, f;
    double cv = 0;      // total (normalised) constraint violation, 0 = feasible
    int rank = 0;
    double crowd = 0;
    double fit = 0;
};

bool dominates(const Vec& a, const Vec& b);   // Pareto dominance (minimisation)
bool cdom(const Ind& a, const Ind& b);        // constrained dominance (Deb's 3 rules when g_deb)

std::vector<std::vector<int>> ndsGeneric(int n, const std::function<bool(int, int)>& dom);
std::vector<std::vector<int>> fastNonDominatedSort(const std::vector<Ind>& P);   // uses cdom
std::vector<std::vector<int>> fastNonDominatedSortF(const std::vector<Vec>& F);  // plain dominance
std::vector<Vec> nondominatedSubset(const std::vector<Vec>& F);
std::vector<double> crowdingDistance(const std::vector<Vec>& F);                  // F = one front

double hv2d(const std::vector<Vec>& pts, const Vec& ref);
double hvMC(const std::vector<Vec>& pts, const Vec& ref, int samples = 20000, uint64_t sd = 12345);
double hypervolume(const std::vector<Vec>& pts, const Vec& ref);  // exact for M=2, MC otherwise
double igd(const std::vector<Vec>& trueFront, const std::vector<Vec>& A);
double spacing(const std::vector<Vec>& A);   // Schott spacing (std of nearest-neighbour L1 distances)

std::vector<Vec> dasDennis(int M, int p);    // simplex-lattice weights / reference points
bool solveLinear(std::vector<Vec> A, Vec b, Vec& x);  // Gaussian elimination with pivoting
}  // namespace moo
