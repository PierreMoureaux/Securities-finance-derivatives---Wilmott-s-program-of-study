#include "FDExplicitTRS.h"

FDExplicitTRS::FDExplicitTRS(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
    double _rmax, double _M, double _N,
    std::shared_ptr<FDExplicitBond> _Bond, double _K, double _sTRS, bool _is_perfreceiver)
    : FiniteDifferences(_r0, _T, _sigma, _alpha, _beta, _rmin, _rmax, _M, _N)
{
    Bond = _Bond;
    K = _K;
    sTRS = _sTRS;
    is_perfreceiver = _is_perfreceiver;
    Bond->price();
}

void FDExplicitTRS::_setup_boundary_conditions_()
{
    auto index = static_cast<int>(Bond->getN() * T / Bond->getT());
    if (is_perfreceiver)
    {
        for (auto x = 0;x!=(N+1);x++)
        {
            grid(x, N) = Bond->getGrid()(x, index) - K - sTRS * T;
        }
    }
    else
    {
        for (auto x = 0; x != (N + 1); x++)
        {
            grid(x, N) = -Bond->getGrid()(x, index) + K + sTRS * T;
        }
    }
}

void FDExplicitTRS::_traverse_grid_()
{
    auto cps = Bond->getCouponSchedule();
    auto cp = Bond->getCoupon();
    for (auto j = N - 1; j != 0; j--)
    {
        if (std::find(cps.begin(), cps.end(), j) != cps.end())
        {
            for (auto i = 0; i != (M + 1); i++)
            {
                if (i == 0)
                {
                    grid(i, j) = a0[i] * grid(i, j + 1) + b0[i] * grid(i + 1, j + 1) + c0 * grid(i + 3, j + 1) + cp;
                }
                else if (i == M)
                {
                    grid(i, j) = aM[i] * grid(i, j + 1) + bM[i] * grid(i - 1, j + 1) + cM * grid(i - 3, j + 1) + cp;
                }
                else
                {
                    grid(i, j) = a[i] * grid(i - 1, j + 1) + b[i] * grid(i, j + 1) + c[i] * grid(i + 1, j + 1) + cp;
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
void FDExplicitTRS::substitution(std::shared_ptr<FDExplicitBond> _Bond)
{
    Bond = _Bond;
    Bond->price();
}
