#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;
#define PHI 3.141592

int matrixMultiplication(double *, double *, double *, int, int, int);
void swaprow(double *, double *, int, int);
void swapcol(double *, double *, int, int);
void swap(double& a, double& b);
void printMatrix(double *);

int main() {

	double *A, *B, *C;
	int m = 8, n = 8, p = 8;


	A = new double[m*n];
	B = new double[m*n];
	C = new double[m*n];

	const double root2 = sqrt(2);

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i == 0) {
				*(A + j) = (0.5 / root2)*cos((PHI*(2 * j + 1)*i) / 16);
			}

			else {
				*(A + i * m + j) = 0.5*cos((PHI*(2 * j + 1)*i) / 16);
			}
		}
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (j == 0) {
				*(B + i * m + j) = (0.5 / root2)*cos((PHI*(2 * i + 1)*j) / 16);
			}

			else {
				*(B + i * m + j) = 0.5*cos((PHI*(2 * i + 1)*j) / 16);
			}
		}
	}




	matrixMultiplication(A, B, C, m, n, p);
	cout << "Matirx C" << endl;
	printMatrix(C);

	double *Av, *Bv, *Cv;

	Av = new double[m*n];
	Bv = new double[m*n];
	Cv = new double[m*n];


	swaprow(A, Av, 2, 5);
	swapcol(B, Bv, 3, 4);

	cout << "A row를 바꾼 Matrix Av" << endl;
	printMatrix(Av);
	cout << "B row를 바꾼 Matrix Bv" << endl;
	printMatrix(Bv);

	matrixMultiplication(Av, Bv, Cv, m, n, p);
	cout << "Av와 Bv의 곱인 Matrix Cv" << endl;
	printMatrix(Cv);


	return 0;

}






int matrixMultiplication(double *a, double *b, double *c, int m, int n, int p) {
	double dum;
	for (int j = 0; j < m; j++) {
		for (int k = 0; k < p; k++) {
			dum = 0;
			for (int l = 0; l < n; l++) dum += a[j*n + l] * b[l*p + k];
			c[j*p + k] = dum;
		}
	}

	return 1;


}

void swaprow(double *a, double *b, int row1, int row2) {

	for (int i = 0; i < 8; i++) {
		swap(a[8 * (row1 - 1) + i], a[8 * (row2 - 1) + i]);
	}

	for (int i = 0; i < 8 * 8; i++)
		b[i] = a[i];
}

void swapcol(double *a, double *b, int col1, int col2) {

	for (int i = 0; i < 8; i++) {
		swap(a[8 * i + (col1 - 1)], a[8 * i + (col2 - 1)]);
	}


	for (int i = 0; i < 8 * 8; i++)
		b[i] = a[i];
}

void swap(double& a, double& b) {
	double temp = 0;
	temp = a;
	a = b;
	b = temp;
}

void printMatrix(double *a) {
	int length = sizeof(a[0]);

	cout << endl;

	for (int i = 0; i < length; i++) {
		for (int j = 0; j < length; j++) {
			cout << setw(14) << a[i*length + j];
		}
		cout << "\n" << endl;
	}

	cout << "\n";
}