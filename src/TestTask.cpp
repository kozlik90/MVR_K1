#include "TestTask.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

TestTask::TestTask(int n_, int m_, double eps_, double omega_, unsigned int Nmax_) : n(n_),
m(m_), eps(eps_), omega(omega_), Nmax(Nmax_), hx(1.0/n), ky(1.0/m), iter(0), maxDiff(0.0), nevyazkaMax(0.0), maxError(0.0){
	u = new double* [n + 1];
	Error = new double* [n + 1];
	for (int i = 0; i < n + 1; i++) {
		u[i] = new double[m + 1];
		Error[i] = new double[m + 1];
	}

	if (omega == 2) {
		omega = 2.0 / (1 + 2 * sin(PI * hx / 2));
	}
}

double TestTask::u_func(double x, double y) {
	return exp(pow(sin(PI * x * y), 2));
}

double TestTask::f_func(double x, double y) {
	return -2 * PI * PI * exp(pow(sin(PI * x * y), 2)) * (y * y * pow(cos(PI * x * y), 2) + x * x * pow(cos(PI * x * y), 2) + pow(sin(PI * x * y), 2) * cos(2 * PI * x * y) * (y * y + x * x));
}

void TestTask::set_GU() {
	// по горизонталям
	for (int i = 0; i <= n; i++) {
		double x = i * hx;
		u[i][0] = u_func(x, 0.0);
		u[i][m] = u_func(x, 1.0);
	}
	// по вертикалям
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
			// Линейная интерполяция по оси x:
			u[i][j] = u[0][j] * (1 - x) + u[n][j] * x; //
		}
	}
}

void TestTask::calculateIst() {
	std::ofstream outFile1("1istinnoe.csv");
	if (!outFile1.is_open()) {
		cout << "Не удалось открыть файл для записи.\n";
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
	//qDebug() << "\r" << empty << "\r" << flush;
	cout << "Файл 1istinnoe.csv готов\n";
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
		cout << "Не удалось открыть файл для записи.\n";
		//return 1;
	}
	cout << "Файл 2 открыт\n";
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
	cout << "Файл 2chislennoe.csv готов\n";
}

void TestTask::calculateError() {
	double pow_h2 = 1.0 / (hx * hx);
	double pow_k2 = 1.0 / (ky * ky);
	double denom = 2 * pow_h2 + 2 * pow_k2;
	std::ofstream outFile3("3pogreshnost.csv");
	if (!outFile3.is_open()) {
		cout << "Не удалось открыть файл для записи.\n";
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
	cout << "Файл 3pogreshnost.csv готов\n";
	
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
	std::cout << "Для решения тестовой задачи использованы сетка с числом разбиений по x: n = " << n << " и числом разбиений по y: m = " << m << "\n";
	std::cout << "Метод верхней релаксации с параметром omega = " << omega << ", применены критерии остановки по точности eps = " << eps << " и по числу итераций Nmax = " << Nmax << "\n";
	std::cout << "На решение схемы (СЛАУ) затрачено итераций N = " << iter << " и достигнута точность итерационного метода eps^(N) = " << maxDiff << "\n";
	std::cout << "Схема (СЛАУ) решена с невязкой ||R^(N)|| = " << nevyazkaMax << " для невязки СЛАУ использована норма max" << "\n";
	std::cout << "Тестовая задача должна быть решена с погрешностью не более eps = 0.5*10^(-6); задача решена с погрешностью eps1 = " << maxError << "\n";
	//std::cout << "Максимальное отклонение точного и численного решений наблюдается в узле x = " << xMax << "; y = " << yMax << "(x" << xIndex << ",y" << yIndex << ")\n";
	std::cout << "В качестве начального приближения использована интерполяция по x\n";
}
