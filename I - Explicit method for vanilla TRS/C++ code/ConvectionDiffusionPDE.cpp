
#include "pch.h"
#include "ConvectionDiffusionPDE.h"
#include <math.h>

BlackScholesPDE::BlackScholesPDE(std::shared_ptr<VanillaTRS> _TRS) : TRS(_TRS) {}

// Convection coefficient
double BlackScholesPDE::conv_coeff(double t, double x) const {
	return (TRS->r)*x;  // rS
}

// Zero-term coefficient
double BlackScholesPDE::zero_coeff(double t, double x) const {
	return -(TRS->r);  // -r
}

// Initial condition (vanilla TRS)
double BlackScholesPDE::init_cond(double x, double y) const {
	return TRS->pay_off->operator()(x,y);
}

std::shared_ptr<VanillaTRS> BlackScholesPDE::getTRS() const {
	return this->TRS;
}
