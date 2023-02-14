#ifndef __PDE_H
#define __PDE_H

#include "VanillaTRS.h"
#include <memory>

// Convection Diffusion Equation - Second-order PDE with simplified delta one coefficients
class ConvectionDiffusionPDE {
public:
	// PDE Coefficients 
	virtual double conv_coeff(double t, double x) const = 0;
	virtual double zero_coeff(double t, double x) const = 0;

	// Boundary and initial conditions
	virtual double init_cond(double x, double y) const = 0;

	virtual std::shared_ptr<VanillaTRS> getTRS() const = 0;
};

// Black-Scholes PDE
class BlackScholesPDE : public ConvectionDiffusionPDE {
public:
	std::shared_ptr<VanillaTRS> TRS;
	BlackScholesPDE(std::shared_ptr<VanillaTRS> _TRS);

	double conv_coeff(double t, double x) const;
	double zero_coeff(double t, double x) const;

	double init_cond(double x, double y) const;

	std::shared_ptr<VanillaTRS> getTRS() const;
};

#endif