#include "TestTask.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main() {
	setlocale(LC_ALL, "Rus");
	cout << "Выберите задачу:\n1.Тестовая\n2.Основная" << endl;
	int answer;
	do {
		cin >> answer;
		if(answer < 1 || answer > 2)
			cout << "Неверный номер задачи, попробуйте ещё раз:" << endl;
	} while (answer < 1 || answer > 2);

	switch (answer) {
	case 1: {
		int n, m;
		double eps, omega;
		unsigned int Nmax;
		cout << "Введите n: ";
		do {
			cin >> n;
			if (n < 0)
				cout << "Число должно быть больше 0:" << endl;
		} while (n < 0);
		cout << "Введите m: ";
		do {
			cin >> m;
			if (m < 0)
				cout << "Число должно быть больше 0:" << endl;
		} while (m < 0);
		cout << "Введите eps: ";
		do {
			cin >> eps;
			if (eps < 0)
				cout << "Число должно быть больше 0:" << endl;
		} while (eps < 0);
		cout << "Введите omega: ";
		do {
			cin >> omega;
			if (omega < 0 || omega > 2)
				cout << "Число должно быть > 0 и < 2:" << endl;
		} while (omega < 0 || omega > 2);
		cout << "Введите Nmax: ";
		do {
			cin >> Nmax;
			if (Nmax < 0)
				cout << "Число должно быть > 0 и < 2:" << endl;
		} while (Nmax < 0);

		TestTask task(n, m, eps, omega, Nmax);
		task.compute();
		task.printInfo();
	}
	}

	
}