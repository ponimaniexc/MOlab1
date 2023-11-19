#pragma once
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

const int dim = 2; //конфликт с другим классом, если называть "n"

class Marquardt
{
private:
	vector<vector<double>> H;
	vector<double> x, grad, d;
	double lambda, eps, h = 1e-3;
	int k = 0, M;

public:
	Marquardt(vector<double>, double, double, int);
	double Function(vector<double>);
	void Gradient();
	void Hessian();
	void Find_d();
	vector<double> Calculate();
	void Print();
};

