#include "FDExplicitBond.h"

FDExplicitBond::FDExplicitBond(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
    double _rmax, int _M, int _N, const std::vector<double>& _couponSchedule, double _coupon)
    : FiniteDifferences(_r0, _T, _sigma, _alpha, _beta, _rmin, _rmax, _M, _N)
{
    couponSchedule = _couponSchedule;
    coupon = _coupon;
}

void FDExplicitBond::_setup_boundary_conditions_()
{
    Eigen::VectorXd v(N + 1);
    for (auto row = 0; row != (N+1); row++)
    {
        v(row) = 1.0 + coupon;
    }
    grid(Eigen::all, Eigen::last) = v;
}

void FDExplicitBond::_traverse_grid_()
{
    for (auto j = N-1; j != 0; j--)
    {
        if (std::find(couponSchedule.begin(), couponSchedule.end(), j) != couponSchedule.end())
        {
            for (auto i = 0; i != (M+1); i++)
            {
                if (i == 0)
                {
                    grid(i, j) = a0[i] * grid(i, j + 1) + b0[i] * grid(i + 1, j + 1) + c0 * grid(i + 3, j + 1) + coupon;
                }
                else if (i == M)
                {
                    grid(i, j) = aM[i] * grid(i, j + 1) + bM[i] * grid(i - 1, j + 1) + cM * grid(i - 3, j + 1) + coupon;
                }
                else
                {
                    grid(i, j) = a[i] * grid(i-1, j + 1) + b[i] * grid(i, j + 1) + c[i] * grid(i + 1, j + 1) + coupon;
                }
            }
        }
        else
        {
            for (auto i = 0; i != (M + 1); i++)
            {
                if (i == 0)
                {
                    grid(i, j) = a0[i] * grid(i, j + 1) + b0[i] * grid(i + 1, j + 1) + c0 * grid(i + 3, j + 1);
                }
                else if (i == M)
                {
                    grid(i, j) = aM[i] * grid(i, j + 1) + bM[i] * grid(i - 1, j + 1) + cM * grid(i - 3, j + 1);
                }
                else
                {
                    grid(i, j) = a[i] * grid(i - 1, j + 1) + b[i] * grid(i, j + 1) + c[i] * grid(i + 1, j + 1);
                }
            }
        }
    }
}

std::vector<double> FDExplicitBond::getCouponSchedule() const
{
    return couponSchedule;
}

double FDExplicitBond::getCoupon() const
{
    return coupon;
}
