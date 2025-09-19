#pragma once
#include <vector>

class Drawing;

class MainTask
{
private:
    int n, n2, m, m2, iter, iter2;
    std::vector<std::vector<double>> v1, v2, Error;
    double eps, omega, maxDiff, maxDiff2, nevyazkaMax, nevyazkaMax2, maxError;
    double eps2, omega2;
	unsigned int Nmax, Nmax2;
    const double PI = 3.14159265358979323846;
    void set_GU(std::vector<std::vector<double>>&, int, int);
    void set_inter(std::vector<std::vector<double>>&, int, int);
    void calculateChisl(std::vector<std::vector<double>>&, int, int, double, double, unsigned int, int&, double&);
    void calculateError();
    double f_func(double, double);
    void writeChisl1();
    void writeChisl2();
public:
    MainTask(int, int, double, double, unsigned int);
    void compute();
    void printInfo();
    friend class Drawing;
};

