#pragma once
#include <vector>

class Drawing;

class TestTask
{
private:
	int n, m, iter;
	double hx, ky, eps, omega, maxDiff, nevyazkaMax, maxError;
	unsigned int Nmax;
	std::vector<std::vector<double>> u;
	std::vector<std::vector<double>> Error;
	

public:
	TestTask(int, int, double, double, unsigned int);
	void compute();
	void printInfo();
	friend class Drawing;

private:
	static double u_func(double x, double y);
	double f_func(double x, double y);
	void set_GU();
	void set_inter();
	void calculateIst();
	void calculateChisl();
	void calculateError();
};

