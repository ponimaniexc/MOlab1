#include <iostream>
#include <iomanip>
#include "GradientDescent.h"
#include "Marquardt.h"

using namespace std;

double Function(vector<double> u)
{
	return (u[0] + 6. * u[1]) * (u[0] + 6. * u[1]) + (u[0] + 2.) * (u[0] + 2.);
}

int main()
{
	vector<double> x = { 1., 2. };
	double t = 1., eps = 1e-10, lambda = 1e3;
	double M = 1e5;

	cout << setprecision(5);

	GradientDescent result_1(x.size(), Function, x, t, eps, M);
	result_1.Calculate();
	result_1.Print();

	Marquardt result_2(x.size(), Function, x, lambda, eps, M);
	result_2.Calculate();
	result_2.Print();

	return 0;
};