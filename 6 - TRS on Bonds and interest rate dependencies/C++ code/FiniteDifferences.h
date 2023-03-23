#pragma once
#include <vector>
#include <Eigen>

class FiniteDifferences
{
public:
	FiniteDifferences(double _r0, double _T, double _sigma, double _alpha, double _beta, double _rmin,
		double _rmax, int _M, int _N);
	virtual void _setup_boundary_conditions_() = 0;
	virtual void _setup_coefficients_() final;
	virtual void _traverse_grid_() = 0;
	virtual double _interpolate_() const final;
	virtual double price() final;
	virtual int getN() const final;
	virtual double getT() const final;
	virtual Eigen::MatrixXd getGrid() const final;

protected:
	double r0;
	double T;
	double sigma;
	double alpha;
	double beta;
	double rmin;
	double rmax;
	int M;
	int N;
	std::vector<double> boundary_conds;
	double dr;
	double dt;
	std::vector<double> j_values;
	Eigen::MatrixXd grid;
	std::vector<double> a;
	std::vector<double> b;
	std::vector<double> c;
	std::vector<double> a0;
	std::vector<double> b0;
	double c0;
	std::vector<double> aM;
	std::vector<double> bM;
	double cM;
};

