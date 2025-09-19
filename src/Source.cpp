#include "TestTask.h"
#include "MainTask.h"
#include "Drawing.h"
#include <iostream>
#include <iomanip>



using namespace std;

int main(int argc, char** argv) {
	setlocale(LC_ALL, "Rus");
	cout << "Choose task:\n1.Test\n2.Main" << endl;
	int answer;
	do {
		cin >> answer;
		if(answer < 1 || answer > 2)
			cout << "Wrong number, try again:" << endl;
	} while (answer < 1 || answer > 2);

	switch (answer) {
	case 1: {
		int n, m;
		double eps, omega;
		unsigned int Nmax;
		cout << "n: ";
		do {
			cin >> n;
			if (n < 0)
				cout << "n must be > 0, try again:" << endl;
		} while (n < 0);
		cout << "m: ";
		do {
			cin >> m;
			if (m < 0)
				cout << "m must be > 0, try again:" << endl;
		} while (m < 0);
		cout << "eps: ";
		do {
			cin >> eps;
			if (eps < 0)
				cout << "eps must be > 0, try again:" << endl;
		} while (eps < 0);
		cout << "omega: ";
		do {
			cin >> omega;
			if (omega < 0 || omega > 2)
				cout << "omega must be > 0 and < 2, try again:" << endl;
		} while (omega < 0 || omega > 2);
		cout << "Nmax: ";
		do {
			cin >> Nmax;
			if (Nmax < 0)
				cout << "Nmax must be > 0, try again:" << endl;
		} while (Nmax < 0);

		TestTask task(n, m, eps, omega, Nmax);
		task.compute();
		task.printInfo();
		Drawing draw1(&task);
		draw1.draw(&argc, argv);
		break;
	}
	case 2: {
		int n, m;
		double eps, omega;
		unsigned int Nmax;
		cout << "n: ";
		do {
			cin >> n;
			if (n < 0)
				cout << "n must be > 0, try again:" << endl;
		} while (n < 0);
		cout << "m: ";
		do {
			cin >> m;
			if (m < 0)
				cout << "m must be > 0, try again:" << endl;
		} while (m < 0);
		cout << "eps: ";
		do {
			cin >> eps;
			if (eps < 0)
				cout << "eps must be > 0, try again:" << endl;
		} while (eps < 0);
		cout << "omega: ";
		do {
			cin >> omega;
			if (omega < 0 || omega > 2)
				cout << "omega must be > 0 and < 2, try again:" << endl;
		} while (omega < 0 || omega > 2);
		cout << "Nmax: ";
		do {
			cin >> Nmax;
			if (Nmax < 0)
				cout << "Nmax must be > 0, try again:" << endl;
		} while (Nmax < 0);
		MainTask task1(n, m, eps, omega, Nmax);
		task1.compute();
		task1.printInfo();
	}
	}

	int a;
	cin >> a;

	
}