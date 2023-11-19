#pragma once
#include <iostream>
#include <vector>
#include <math.h>

using namespace std;

const int n = 2;

class GradientDescent
{
private:
	vector<double> x, grad;
	double t, eps, h = 1e-5;
	int k = 0, M;

public:
	GradientDescent(vector<double>, double, double, int);
	double Function(vector<double>);
	void Gradient();
	bool StopCriteria(vector<double>);
	vector<double> Calculate();
	void Print();
};

