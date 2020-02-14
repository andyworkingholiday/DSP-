#include <iostream>
#include <cmath>
#include <string>
#include "complex.h"
using namespace std;
#define PHI 3.141592

void print_complex(complex x, string a); // 복소수를 "a+bj" 또는 "a-bj"로 출력하기 위한 함수 
void square_equation(double a, double b, double c); // 2차 방정식의 계수를 입력하면 그 해를 구하기 위한 함수

int main() {

	// #2번문제 //

	cout << "---------- 2번문제 ----------" << endl << endl;
	complex x, y, a, b, c, d, e, f;
	a = complex(2, 1);
	b = complex(3, -2);
	c = complex(3, 0);
	d = complex(3, -5);
	e = complex(2, 4);
	f = complex(1, 1);

	x = (e * d - f * b) / (a * d - b * c);
	y = (a * f - c * e) / (a * d - b * c);


	print_complex(x, "x");
	print_complex(y, "y");
	cout << "\n";

	// #3번 문제 //

	cout << "---------- 3번문제 ----------" << endl << endl;

	int i;
	int m = 3, n = 3, N = 8;
	complex com1, com2;
	double result1;
	double sum = 0;
	for (i = 0; i < N; i++) {
		com1 = (complex((2. * PHI*m*i) / N)) / sqrt(N);
		com2 = (complex((-2. * PHI*n*i) / N)) / sqrt(N);
		result1 = (com1.re)*(com2.re) - (com1.im)*(com2.im);
		sum += result1;
	}

	cout << sum << endl << "\n";


	// #4번 문제 //

	cout << "---------- 4번문제 ----------" << endl << endl;
	square_equation(3, 2, 7);
	cout << "\n";

}

void print_complex(complex x, string a) {
	if ((x.im) > 0)
		cout << a << " = " << x.re << "+" << x.im << "j" << endl;

	else if ((x.im) < 0)
		cout << a << " = " << x.re << x.im << "j" << endl;

	else
		cout << a << " = " << x.re << endl;
}

void square_equation(double a, double b, double c) {
	double pan = b * b - 4 * a*c;
	double result1, result2;
	complex answer1, answer2;

	if (pan > 0) {
		result1 = (-b + sqrt(pan)) / (2 * a);
		result2 = (-b - sqrt(pan)) / (2 * a);

		cout << "첫 번째 근 = " << result1 << endl;
		cout << "두 번째 근 = " << result2 << endl;
	}

	else if (pan < 0) {
		answer1.re = -b / (2 * a);
		answer1.im = sqrt(-pan) / (2 * a);
		answer2.re = -b / (2 * a);
		answer2.im = -answer1.im; //켤레복소수

		print_complex(answer1, "첫 번째 근");
		print_complex(answer2, "두 번째 근");
	}

	else {
		//중근일 때
		result1 = -b / (2 * a);
		cout << "중근 = " << result1 << endl;
	}
}