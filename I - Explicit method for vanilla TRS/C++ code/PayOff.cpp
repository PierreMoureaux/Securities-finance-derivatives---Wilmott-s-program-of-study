
#include "pch.h"
#include "PayOff.h"

PayOff::PayOff() {}

// ==========
// PayOffPerfReceiver
// ==========

PayOffPerfReceiver::PayOffPerfReceiver(const double& _K, const double& _r0) {
	K = _K;
	r0 = _r0;
}

double PayOffPerfReceiver::operator() (const double& S, const double& T) const {
	return S - K - r0*K*T; // Standard performance receiver TRS pay-off
}

double PayOffPerfReceiver::getr0() const {
	return this->r0;
}

// =========
// PayOffPerfPayer
// =========

PayOffPerfPayer::PayOffPerfPayer(const double& _K, const double& _r0) {
	K = _K;
	r0 = _r0;
}

double PayOffPerfPayer::operator() (const double& S, const double& T) const {
	return -(S - K - r0 * K*T); // Standard performance payer TRS pay-off
}

double PayOffPerfPayer::getr0() const {
	return this->r0;
}
