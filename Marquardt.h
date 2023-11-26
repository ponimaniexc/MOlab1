#pragma once
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

class Marquardt
{
private:
	double(*f)(vector<double>);
	vector<vector<double>> H;
	vector<double> x, grad, d;
	double lambda, eps, h = 1e-3, M;
	int k = 0;
	size_t dim;

public:
	Marquardt(size_t, double(*)(vector<double>), vector<double>, double, double, double);
	void Gradient();
	void Hessian();
	void Find_d();
	vector<double> Calculate();
	void Print();
};

