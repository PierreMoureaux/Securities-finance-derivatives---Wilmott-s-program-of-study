#include "FiniteDifferences.h"
#include <numeric>
#include <cmath>

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N - 1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

template<typename Iter_T> 
double vectorNorm(Iter_T first, Iter_T last) {
    return sqrt(inner_product(first, last, first, 0.0L));
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

FiniteDifferences::FiniteDifferences(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
    double _rmax, int _M, int _N)
	: r0(_r0), T(_T), sigma(_sigma), alpha(_alpha), beta(_beta), rmin(_rmin), rmax(_rmax), M(_M), N(_N)
{
    boundary_conds = linspace<double>(rmin, rmax, M + 1);

    //dr definition
    auto s = boundary_conds.size();
    std::vector<double> vecnorm(s);
    std::transform(boundary_conds.begin(), boundary_conds.end(), vecnorm.begin(),[&](double x) {return alpha * (beta - x); });
    auto testdr1 = (rmax - rmin) / static_cast<double>(M);
    auto testdr2 = (pow(sigma, 2) / (vectorNorm(vecnorm.begin(), vecnorm.end())));
    dr = (testdr1 > testdr2) ? testdr2 : testdr1;

    //dt definition
    auto testdt1 = T / static_cast<double>(N);
    auto testdt2 = pow(dr, 2) / pow(sigma, 2);
    dt = (testdt1 > testdt2) ? testdt2 : testdt1;

    j_values = arange<double>(0, N, 1);

    grid = Eigen::MatrixXd(M + 1, N + 1);
}

void FiniteDifferences::_setup_coefficients_()
{
    auto s = boundary_conds.size();
    a.resize(s);
    b.resize(s);
    c.resize(s);
    a0.resize(s);
    b0.resize(s);
    aM.resize(s);
    bM.resize(s);
    std::transform(boundary_conds.begin(), boundary_conds.end(), a.begin(), [&](double x) {return
        dt / 2 * (pow((sigma / dr),2) - alpha * (beta - x) / dr); });
    std::transform(boundary_conds.begin(), boundary_conds.end(), b.begin(), [&](double x) {return
        -(pow(sigma,2)) * dt / (pow(dr,2)) - x * dt + 1; });
    std::transform(boundary_conds.begin(), boundary_conds.end(), c.begin(), [&](double x) {return
        dt / 2 * (pow((sigma / dr),2) + alpha * (beta - x) / dr); });
    std::transform(boundary_conds.begin(), boundary_conds.end(), a0.begin(), [&](double x) {return
        dt * (-(pow((sigma / dr),2)) / 3 - 1 / dr * alpha * (beta - x) - x) + 1; });
    std::transform(boundary_conds.begin(), boundary_conds.end(), b0.begin(), [&](double x) {return
        dt / dr * ((pow(sigma,2)) / (2 * dr) + alpha * (beta - x)); });
    c0 = -dt * (pow(sigma,2)) / (6 * (pow(dr,2)));
    std::transform(boundary_conds.begin(), boundary_conds.end(), aM.begin(), [&](double x) {return
        dt * ((pow((sigma / dr), 2)) / 3 + 1 / dr * alpha * (beta - x) - x) + 1; });
    std::transform(boundary_conds.begin(), boundary_conds.end(), bM.begin(), [&](double x) {return
        dt / dr * (-(pow(sigma, 2)) / (2 * dr) - alpha * (beta - x)); });
    cM = dt * (pow(sigma, 2)) / (6 * (pow(dr, 2)));
}

double FiniteDifferences::_interpolate_() const
{
    auto downPrice = grid(Eigen::all, Eigen::seq(0, 0));
    auto lower = lower_bound(boundary_conds.begin(), boundary_conds.end(), r0) - boundary_conds.begin() -1;
    auto upper = upper_bound(boundary_conds.begin(), boundary_conds.end(), r0) - boundary_conds.begin();
    double x1 = grid(1,lower);
    double x2 = grid(1, upper);
    return std::lerp(x1,x2,r0);
}

double FiniteDifferences::price()
{
    _setup_boundary_conditions_();
    _setup_coefficients_();
    _traverse_grid_();
    return _interpolate_();
}

int FiniteDifferences::getN() const
{
    return N;
}

double FiniteDifferences::getT() const
{
    return T;
}

Eigen::MatrixXd FiniteDifferences::getGrid() const
{
    return grid;
}