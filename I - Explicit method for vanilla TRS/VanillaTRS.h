#ifndef __VANILLA_OPTION_H
#define __VANILLA_OPTION_H

#include "PayOff.h"
#include <memory>

class VanillaTRS {
public:
	std::shared_ptr<PayOff> pay_off;

	double K;
	double r;
	double T;

	VanillaTRS();
	VanillaTRS(double _K, double _r, double _T, std::shared_ptr<PayOff> _pay_off);

	std::shared_ptr<PayOff> getPayOff() const;
};

#endif