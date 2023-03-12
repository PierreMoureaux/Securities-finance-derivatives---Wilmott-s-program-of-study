
#include "pch.h"
#include "VanillaTRS.h"

VanillaTRS::VanillaTRS() {}

VanillaTRS::VanillaTRS(double _K, double _r, double _T, std::shared_ptr<PayOff> _pay_off) :
	K(_K), r(_r), T(_T), pay_off(std::move(_pay_off)) {}

std::shared_ptr<PayOff> VanillaTRS::getPayOff() const {
	return this->pay_off;
}