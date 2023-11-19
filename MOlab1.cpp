#include <iostream>
#include <iomanip>
#include "GradientDescent.h"
#include "Marquardt.h"

using namespace std;

int main()
{
	vector<double> x = { 1., 2. };
	double t = 1., eps = 1e-10, lambda = 1e3;
	int M = 1e5;

	cout << setprecision(5);

	GradientDescent result_1(x, t, eps, M);
	result_1.Calculate();
	result_1.Print();

	Marquardt result_2(x, lambda, eps, M);
	result_2.Calculate();
	result_2.Print();

	return 0;
};