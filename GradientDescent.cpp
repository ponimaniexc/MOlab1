#include "GradientDescent.h"

GradientDescent::GradientDescent(vector<double> x, double t, double eps, int M)
{
	this->x = x;
	this->t = t;
	this->eps = eps;
	this->M = M;
}

double GradientDescent::Function(vector<double> q)
{
	return (q[0] + 6. * q[1]) * (q[0] + 6. * q[1]) + (q[0] + 2.) * (q[0] + 2.);
}

void GradientDescent::Gradient()
{
	grad.resize(2);

	for (int i = 0; i < n; i++)
	{
		vector<double> x_i_next = x, x_i_prev = x;
		x_i_next[i] += h;
		x_i_prev[i] -= h;
		grad[i] = (Function(x_i_next) - Function(x_i_prev)) / (2 * h);
	}
}

bool GradientDescent::StopCriteria(vector<double> next)
{
	double x_norma = 0.;
	double max_grad = 0.;
	double f_norma = 0.;

	for (int i = 0; i < n; i++)
	{
		x_norma += (next[i] - x[i]) * (next[i] - x[i]);

		if (max_grad < abs(grad[i])) max_grad = abs(grad[i]);

		f_norma += grad[i] * grad[i];
	}

	if (sqrt(x_norma) <= eps * eps)
		return true;

	if (abs(Function(next) - Function(x)) <= eps)
		return true;

	if (max_grad < eps)
		return true;

	if (f_norma < eps)
		return true;

	return false;
}

vector<double> GradientDescent::Calculate()
{
	vector<double> x_next(n);

	while (true)
	{
		Gradient();

		for (int i = 0; i < n; i++)
			x_next[i] = x[i] - t * grad[i];

		if (Function(x_next) - Function(x) < 0)
		{
			if (k + 1 > M || StopCriteria(x_next))
			{
				x = x_next;
				break;
			}
			else
			{
				x = x_next;
				k++;
			}
		}
		else
		{
			t /= 2.;
		}
	}

	return x;
}

void GradientDescent::Print()
{
	cout << "GradientDescent:\t" << "k = " << k << ";\tx = " << '(';
	for (int i = 0; i < n; i++)
	{
		cout << x[i];
		if (i != n - 1)
			cout << ", ";
	}
	cout << ')' << endl;
}
