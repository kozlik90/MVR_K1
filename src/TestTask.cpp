#include "TestTask.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

TestTask::TestTask(int n_, int m_, double eps_, double omega_, unsigned int Nmax_) :  n(n_),
m(m_), u(n+1, std::vector<double>(m+1, 0)), Error(n+1, std::vector<double>(m+1, 0)), eps(eps_), omega(omega_), Nmax(Nmax_), hx(1.0/n), ky(1.0/m),
 iter(0), maxDiff(0.0), nevyazkaMax(0.0), maxError(0.0){
	
	if (omega >= 2) {
		omega = 2.0 / (1 + 2 * sin(M_PI * hx / 2));
	}

}

double TestTask::u_func(double x, double y) {
	return exp(pow(sin(M_PI * x * y), 2));
}

double TestTask::f_func(double x, double y) {
	return -2 * M_PI * M_PI * exp(pow(sin(M_PI * x * y), 2)) * (y * y * pow(cos(M_PI * x * y), 2) + x * x * pow(cos(M_PI * x * y), 2) + pow(sin(M_PI * x * y), 2) * cos(2 * M_PI * x * y) * (y * y + x * x));
}

void TestTask::set_GU() {
	
	for (int i = 0; i <= n; i++) {
		double x = i * hx;
		u[i][0] = u_func(x, 0.0);
		u[i][m] = u_func(x, 1.0);
	}
	
	for (int j = 0; j <= m; j++) {
		double y = j * ky;
		u[0][j] = u_func(0.0, y);
		u[n][j] = u_func(1.0, y);
	}
}

void TestTask::set_inter() {
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			double x = i * hx;
			// linear interpolation x:
			u[i][j] = u[0][j] * (1 - x) + u[n][j] * x; //
		}
	}
}

void TestTask::calculateIst() {
	std::ofstream outFile1("1istinnoe.csv");
	if (!outFile1.is_open()) {
		cout << "File is not open.\n";
		//return 1;
	}
	outFile1 << ";";
	for (int i = 0; i <= n; i++) {
		outFile1 << "x" << i;
		if (i < n)
			outFile1 << ";";
	}
	outFile1 << "\n";

	for (int j = 0; j < m + 1; j++) {
		outFile1 << "y" << j << ";";
		for (int i = 0; i < n + 1; i++) {
			double x = j * hx;
			double y = i * ky;
			outFile1 << std::setprecision(30) << u_func(x, y);
			if (i < n)
				outFile1 << ";";

		}
		outFile1 << "\n";
	}
	outFile1.close();
	std::string empty(10, ' ');
	
	cout << "File 1istinnoe.csv is ready\n";
}

void TestTask::calculateChisl() {
	double pow_h2 = 1.0 / (hx * hx);
	double pow_k2 = 1.0 / (ky * ky);
	double denom = 2 * pow_h2 + 2 * pow_k2;
	do {
		maxDiff = 0.0;

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < m; j++) {

				double x = i * hx;
				double y = j * ky;

				double rhs = f_func(x, y) + pow_h2 * (u[i - 1][j] + u[i + 1][j]) + pow_k2 * (u[i][j - 1] + u[i][j + 1]);
				double u_old = u[i][j];
				double u_new = (1 - omega) * u_old + (omega / denom) * rhs;
				u[i][j] = u_new;
				double diff = fabs(u_new - u_old);
				if (diff > maxDiff) maxDiff = diff;

			}
		}
		iter++;

		//std::cout << "\r" << iter << "     " << maxDiff << std::flush;
	} while (maxDiff > eps && iter < Nmax);

	std::ofstream outFile2("2chislennoe.csv");
	if (!outFile2.is_open()) {
		cout << "File is not open.\n";
		//return 1;
	}

	outFile2 << ";";
	for (int i = 0; i <= n; i++) {
		outFile2 << "x" << i;
		if (i < n)
			outFile2 << ";";
	}
	outFile2 << "\n";

	for (int j = 0; j < m + 1; j++) {
		outFile2 << "y" << j << ";";
		for (int i = 0; i < n + 1; i++) {
			outFile2 << std::setprecision(30) << u[i][j];
			if (i < n)
				outFile2 << ";";

		}
		outFile2 << "\n";
	}
	outFile2.close();
	//cout << "\r" << empty << "\r" << flush;
	cout << "File 2chislennoe.csv is ready\n";
}

void TestTask::calculateError() {
	double pow_h2 = 1.0 / (hx * hx);
	double pow_k2 = 1.0 / (ky * ky);
	double denom = 2 * pow_h2 + 2 * pow_k2;
	std::ofstream outFile3("3pogreshnost.csv");
	if (!outFile3.is_open()) {
		cout << "File is not open.\n";
		//return 1;
	}
	outFile3 << ";";
	for (int i = 0; i <= n; i++) {
		outFile3 << "x" << i;
		if (i < n)
			outFile3 << ";";
	}
	outFile3 << "\n";

	
	double xMax = 0.0;
	double yMax = 0.0;
	int xIndex = 0, yIndex = 0;


	for (int j = 0; j < m + 1; j++) {
		outFile3 << "y" << j << ";";
		double y = j * ky;
		for (int i = 0; i < n + 1; i++) {
			double x = i * hx;
			double err = fabs(u[i][j] - u_func(x, y));
			Error[i][j] = err;
			outFile3 << std::setprecision(30) << err;
			if (i < n)
				outFile3 << ";";

			if (err > maxError) {
				maxError = err;
				xMax = x;
				xIndex = i;
				yMax = y;
				yIndex = j;
			}
		}
		outFile3 << "\n";
	}

	outFile3.close();
	cout << "File 3pogreshnost.csv is ready\n";
	
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			double x = i * hx;
			double y = j * ky;
			double r = -pow_h2 * u[i - 1][j] - pow_k2 * u[i][j - 1]
				+ (2 * pow_h2 + 2 * pow_k2) * u[i][j]
				- pow_h2 * u[i + 1][j] - pow_k2 * u[i][j + 1]
				- f_func(x, y);
			if (fabs(r) > nevyazkaMax) {
				nevyazkaMax = fabs(r);
			}
		}
	}
}
void TestTask::compute() {
	set_GU();
	set_inter();
	calculateIst();
	calculateChisl();
	calculateError();
}

void TestTask::printInfo() {
	std::cout << "The number of partitions by x: n = " << n << " and number of partitions by y: m = " << m << "\n";
	std::cout << "SOR with the parameter omega = " << omega << ", criteria for stopping by accuracy eps = " << eps << " and by the number of iterations Nmax = " << Nmax << "\n";
	std::cout << "iterations spent N = " << iter << "; the achieved accuracy of the iterative method eps^(N) = " << maxDiff << "\n";
	std::cout << "The discrepancy ||R^(N)|| = " << nevyazkaMax << "\n";
	std::cout << "The test task should be solved with an error of no more than eps = 0.5*10^(-6); the task was solved with an error of eps1 = " << maxError << "\n";
	std::cout << "x interpolation is used as an initial approximation.\n";
}
