#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
#define e 2.71828

int main() {

	//1번문제//
	int x[5] = { 1,1,1,1,1 };
	int h1[7] = { 1,0,0,0,0,0,1 };
	int h2[5] = { 1,2,3,2,1 };

	int result1 = 0;
	int y1[11] = { 0, };

	for (int n = 0; n < 11; n++) {
		for (int k = 0; k <= 6; k++) {
			if (n - k < 0)
				result1 += 0;
			else if (n - k > 4)
				result1 += 0;
			else
				result1 += h1[k] * x[n - k];
		}

		y1[n] = result1;
		result1 = 0;
	}

	int result2 = 0;
	int y2[15] = { 0, };

	for (int n = 0; n < 15; n++) {
		for (int k = 0; k <= 10; k++) {
			if (n - k < 0)
				result2 += 0;
			else if (n - k > 4)
				result2 += 0;
			else
				result2 += y1[k] * h2[n - k];
		}

		y2[n] = result2;
		result2 = 0;
	}
	

	for (int i = 0; i < 15; i++)
		cout << y2[i] << endl;

	ofstream fout1("Q1.txt");
	for (int i = 0; i < 15; i++)
		fout1 << i << "\t" << y2[i] << endl;

	cout << endl;

	//2번문제//

	double cosine[10];
	for (int i = 0; i < 10; i++)
		cosine[i] = cos(i); 
	double exponent[10];
	for (int i = 0; i < 10; i++)
		exponent[i] = pow(e, -i);

	double result3 = 0;
	double answers[19] = { 0, };

	for (int n = 0; n < 19; n++) {
		for (int k = 0; k <= 9; k++) {
			if (n - k < 0)
				result3 += 0;
			else if (n - k > 9)
				result3 += 0;
			else
				result3 += cosine[k] * exponent[n - k];
		}

		answers[n] = result3;
		result3 = 0;
	}

	for (int i = 0; i < 19; i++)
		cout << answers[i] << endl;

	ofstream fout2("Q2.txt");
	for (int i = 0; i < 19; i++)
		fout2 << i << "\t" << answers[i] << endl;

	system("pause");

	return 0;

}