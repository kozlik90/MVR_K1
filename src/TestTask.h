#pragma once

class TestTask
{
private:
	int n, m, iter;
	double hx, ky, eps, omega, maxDiff, nevyazkaMax, maxError;
	unsigned int Nmax;
	double** u;
	double** Error;
	const double PI = 3.14159265358979323846;

public:
	TestTask(int, int, double, double, unsigned int);
	void compute();
	void printInfo();

private:
	double u_func(double x, double y);
	double f_func(double x, double y);
	void set_GU();
	void set_inter();
	void calculateIst();
	void calculateChisl();
	void calculateError();
};

