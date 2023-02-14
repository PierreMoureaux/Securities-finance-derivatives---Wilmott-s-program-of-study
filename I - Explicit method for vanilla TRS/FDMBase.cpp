
#include "pch.h"
#include <fstream>
#include "FDMBase.h"
#include <cmath>

FDMBase::FDMBase(double _x_dom, unsigned long _J,
	double _t_dom, unsigned long _N,
	std::shared_ptr<ConvectionDiffusionPDE> _pde)
	: x_dom(_x_dom), J(_J), t_dom(_t_dom), N(_N), pde(_pde) {}

FDMEulerExplicit::FDMEulerExplicit(double _x_dom, unsigned long _J,
	double _t_dom, unsigned long _N,
	std::shared_ptr<ConvectionDiffusionPDE> _pde)
	: FDMBase(_x_dom, _J, _t_dom, _N, _pde) {
	calculate_step_sizes();
	set_initial_conditions();
}

void FDMEulerExplicit::calculate_step_sizes() {
	dx = x_dom / static_cast<double>(J - 1);
	dt = t_dom / static_cast<double>(N - 1);
}

void FDMEulerExplicit::set_initial_conditions() {
	// Spatial settings
	double cur_spot = 0.0;
	double init_rate = pde->getTRS()->getPayOff()->getr0();

	old_result.resize(J, 0.0);
	new_result.resize(J, 0.0);
	x_values.resize(J, 0.0);

	for (unsigned long j = 0; j < J; j++) {
		cur_spot = static_cast<double>(j)*dx;
		old_result[j] = pde->init_cond(cur_spot, init_rate);
		x_values[j] = cur_spot;
	}

	// Temporal settings
	prev_t = 0.0;
	cur_t = 0.0;
}

void FDMEulerExplicit::calculate_inner_domain() {
	for (unsigned long j = 0; j < J; j++) {

		// Differencing coefficients (see a and b in text)
		a = dt * (pde->conv_coeff(prev_t, x_values[j]));
		b = 1 - pde->zero_coeff(prev_t, x_values[j])*dt;

		// Update inner values of spatial discretisation grid (Explicit Euler)
		new_result[j] = a + b*old_result[j];
	}
}

void FDMEulerExplicit::step_march() {
	std::ofstream fdm_out("fdm.csv");

	while (cur_t < t_dom) {
		cur_t = prev_t + dt;
		calculate_inner_domain();
		for (int j = 0; j < J; j++) {
			fdm_out << x_values[j] << " " << prev_t << " " << new_result[j] << std::endl;
		}

		old_result = new_result;
		prev_t = cur_t;
	}

	fdm_out.close();
}
