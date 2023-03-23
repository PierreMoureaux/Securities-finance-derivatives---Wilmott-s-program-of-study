#pragma once
#include "FiniteDifferences.h"
class FDExplicitBond :
    public FiniteDifferences
{
public:
	FDExplicitBond(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
		double _rmax, int _M, int _N, const std::vector<double>& _couponSchedule, double _coupon);
	virtual void _setup_boundary_conditions_();
	virtual void _traverse_grid_();
	virtual std::vector<double> getCouponSchedule() const;
	virtual double  getCoupon() const;

protected:
	std::vector<double> couponSchedule;
	double coupon;
};

