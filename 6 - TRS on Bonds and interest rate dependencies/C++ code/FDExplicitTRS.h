#pragma once
#include "FiniteDifferences.h"
#include "FDExplicitBond.h"

class FDExplicitTRS :
    public FiniteDifferences
{
public:
	FDExplicitTRS(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
		double _rmax, double _M, double _N, 
		std::shared_ptr<FDExplicitBond> _Bond, double _K, double _sTRS, bool _is_perfreceiver = true);
	virtual void _setup_boundary_conditions_();
	virtual void _traverse_grid_();
	virtual void substitution(std::shared_ptr<FDExplicitBond> _Bond);

protected:
	//shared pointer because the Bond underlying can of course be used by another TRS
	std::shared_ptr<FDExplicitBond> Bond;
	double K;
	double sTRS;
	bool is_perfreceiver;
};

