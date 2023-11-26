#include "Marquardt.h"

Marquardt::Marquardt(size_t dim, double(*f)(vector<double>), vector<double> x, double lambda, double eps, double M)
{
	this->dim = dim;
	this->f = f;
	this->x = x;
	this->lambda = lambda;
	this->eps = eps;
	this->M = M;
}

void Marquardt::Gradient()
{
	grad.resize(dim);

	for (int i = 0; i < dim; i++)
	{
		vector<double> x_i_next = x, x_i_prev = x;
		x_i_next[i] += h;
		x_i_prev[i] -= h;
		grad[i] = (f(x_i_next) - f(x_i_prev)) / (2 * h);
	}
}

void Marquardt::Hessian()
{
	H.resize(dim, vector<double>(dim));

	//требует доработки для общего случая (тут для n = 2)

	vector<double> x_1 = x, x_2 = x, x_3 = x, x_4 = x;

	x_1[0] += h;
	x_2[0] -= h;
	H[0][0] = (f(x_1) - 2 * f(x) + f(x_2)) / (h * h);

	x_3[1] += h;
	x_4[1] -= h;
	H[1][1] = (f(x_3) - 2 * f(x) + f(x_4)) / (h * h);

	x_1 = x;  x_1[0] += h; x_1[1] += h;
	x_2 = x;  x_2[0] += h; x_2[1] -= h;
	x_3 = x;  x_3[0] -= h; x_3[1] += h;
	x_4 = x;  x_4[0] -= h; x_4[1] -= h;

	H[0][1] = (f(x_1) - f(x_2) - f(x_3) + f(x_4)) / (4 * h * h);
	H[1][0] = H[0][1];
}

void Marquardt::Find_d()
{
	vector<vector<double>> temp, invers;
	d.resize(dim);
	temp.resize(dim, vector<double>(dim));
	invers.resize(dim, vector<double>(dim));

	for (int i = 0; i < dim; i++)
	{
		for (int j = 0; j < dim; j++)
		{
			temp[i][j] = H[i][j];
			if (i == j) temp[i][j] += lambda;
		}
	}

	//требует доработки для общего случая (тут для n = 2)

	double det = temp[0][0] * temp[1][1] - temp[0][1] * temp[1][0];

	invers[0][0] = temp[1][1] / det;
	invers[0][1] = -temp[0][1] / det;
	invers[1][0] = -temp[1][0] / det;
	invers[1][1] = temp[0][0] / det;

	for (int i = 0; i < dim; i++)
	{
		double sum = 0.;

		for (int j = 0; j < dim; j++)
			sum += invers[i][j] * grad[j];

		d[i] = -sum;
	}
}

vector<double> Marquardt::Calculate()
{
	vector<double> x_next(dim);

	while (k < M)
	{
		double f_norma = 0.;

		Gradient();
		Hessian();

		for (int i = 0; i < dim; i++)
			f_norma += grad[i] * grad[i];

		if (sqrt(f_norma) < eps) break;

		while (true)
		{
			Find_d();

			for (int i = 0; i < dim; i++)
			{
				x_next[i] = x[i] + d[i];
			}

			if (f(x_next) - f(x) < 0)
			{
				x = x_next;
				lambda /= 2.;
				k++;
				break;
			}
			else
			{
				lambda *= 2.;
			}
		}
	}

	return x;
}

void Marquardt::Print()
{
	cout << "Marquardt:\t\t" << "k = " << k << ";\t\tx = " << '(';
	for (int i = 0; i < dim; i++)
	{
		cout << x[i];
		if (i != dim - 1)
			cout << ", ";
	}
	cout << ')' << endl;
}
