#ifndef __PAY_OFF_HPP
#define __PAY_OFF_HPP

#include "pch.h"

class PayOff {
public:
	PayOff();
	virtual ~PayOff() {};
	virtual double operator() (const double& S, const double& T) const = 0;
	virtual double getr0() const = 0;
};

class PayOffPerfReceiver : public PayOff {
private:
	double K; // Strike price (also called settlement price)
	double r0; // agreed TRS financing rate
public:
	PayOffPerfReceiver(const double& _K, const double& _r0);
	virtual ~PayOffPerfReceiver() {};
	virtual double operator() (const double& S, const double& T) const;
	double getr0() const;
};

class PayOffPerfPayer : public PayOff {
private:
	double K; // Strike (also called settlement price)
	double r0; // agreed TRS financing rate
public:
	PayOffPerfPayer(const double& _K, const double& _r0);
	virtual ~PayOffPerfPayer() {};
	virtual double operator() (const double& S, const double& T) const;
	double getr0() const;
};

#endif

