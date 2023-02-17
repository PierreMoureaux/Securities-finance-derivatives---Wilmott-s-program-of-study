#include "pch.h"
#include "PayOff.h"
#include "VanillaTRS.h"
#include "ConvectionDiffusionPDE.h"
#include "FDMBase.h"
#include <memory>

int main(int argc, char **argv) {
	// Create the option parameters
	double K = 50;  // Strike price
	double r0 = 0.04;  // agreed TRS financing rate
	double r = 0.01;   // Risk-free rate (1%)
	double T = 1.00;    // One year until expiry

	// FDM discretisation parameters
	double x_dom = 100;       // Spot goes from [0.0, 200]
	unsigned long J = 20;
	double t_dom = T;         // Time period as for the TRS
	unsigned long N = 20;

	// Create the PayOff and Option objects
	std::shared_ptr<PayOff> Pay_Off_Perf_Receiver = std::make_shared<PayOffPerfReceiver>(K, r0);
	std::shared_ptr<VanillaTRS> Vanilla_TRS = std::make_shared<VanillaTRS>(K, r, T, Pay_Off_Perf_Receiver);

	// Create the PDE and FDM objects
	std::shared_ptr<BlackScholesPDE> bs_pde = std::make_shared<BlackScholesPDE>(Vanilla_TRS);
	FDMEulerExplicit fdm_euler(x_dom, J, t_dom, N, bs_pde);

	// Run the FDM solver
	fdm_euler.step_march();

	return 0;
}