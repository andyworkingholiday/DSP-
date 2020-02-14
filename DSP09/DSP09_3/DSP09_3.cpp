#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;
#define PHI 3.141592
#define W 256
#define H 256
#define uchar unsigned char
#define Datalength 256

void DCT_2D(uchar**, double**);
void IDCT_2D(double**, uchar**);
void make_zero1(double**);
void make_zero2(double**);
void make_zero3(double**);



int main() {

	ifstream noise;
	noise.open("twin_.bmp", ios::binary);

	const int headersize = 54;
	char* header = new char[headersize];
	noise.read((char*)header, headersize);

	/*unsigned char**A;
	A = new unsigned char*[H];
	for (int i = 0; i < H; i++)
		A[i] = new unsigned char[W];*/

	unsigned char** R;
	unsigned char** G;
	unsigned char** B;
	unsigned char** RGB;

	double** dct;
	R = new uchar*[H];
	G = new uchar*[H];
	B = new uchar*[H];
	RGB = new uchar*[H]; 

	dct = new double*[H];

	for (int i = 0; i < H; i++) {
		R[i] = new uchar[W];
		G[i] = new uchar[W];
		B[i] = new uchar[W];
		RGB[i] = new uchar[3*W];
		dct[i] = new double[W];
	}

	for (int i = 0; i < H; i++)
		noise.read((char*)RGB[i], 3 * W);


	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			R[i][j] = RGB[i][jj];
			G[i][j] = RGB[i][jj + 1];
			B[i][j] = RGB[i][jj + 2];
		}
	}

	DCT_2D(R, dct);

	ofstream fff;
	fff.open("dct.txt", ios::binary);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			fff << dct[i][j] << "\t";
		}
		fff << endl;
	}
	

	make_zero3(dct);

	ofstream fff3;
	fff3.open("zerodct.txt", ios::binary);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			fff3 << dct[i][j] << "\t";
		}
		fff3 << endl;
	}

	unsigned char** result_data;
	result_data = new unsigned char*[H];
	for (int i = 0; i < H; i++)
		result_data[i] = new unsigned char[W];


	

	IDCT_2D(dct, result_data);

	ofstream fff2;
	fff2.open("idct.txt", ios::binary);
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			fff2 << result_data[i][j] << "\t";
		}
		fff2 << endl;
	}
	fff2.close();
	

	unsigned char** output_data;
	output_data = new unsigned char*[H];
	for (int i = 0; i < H; i++)
		output_data[i] = new unsigned char[3 * W];


	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			output_data[i][jj] = result_data[i][j];
			output_data[i][jj + 1] = result_data[i][j];
			output_data[i][jj + 2] = result_data[i][j];

		}
	}

	
	
	ofstream clean;
	clean.open("output_twin.bmp", ios::binary);
	clean.write((char*)header, 54);

	
	for (int i = 0; i < H; i++)
		clean.write((char*)output_data[i], 3*W);

	system("pause");
	return 0;





 }

void DCT_2D(uchar** data, double** dct) {
	int M = 8;
	int N = 8;
	int mcrNb = Datalength / N;
	
	double sum = 0;

	int u, v, n, m = 0;
	for (int mcr1 = 0; mcr1 < mcrNb; mcr1++) {
		for (int mcr2 = 0; mcr2 < mcrNb; mcr2++) {
			for (int k = 0; k < M; k++) {
				for (int l = 0; l < N; l++) {
					u = mcr1 * M + k;
					v = mcr2 * N + l;
					sum = 0;

					for (int i = 0; i < M; i++) {
						for (int j = 0; j < N; j++) {
							m = mcr1 * M + i;
							n = mcr2 * N + j;
							double theta1 = (double)(2.0*i + 1)*k*PHI / (2.0*M);
							double theta2 = (double)(2.0*j + 1)*l*PHI / (2.0*N);
							sum += (double)cos(theta1)*(double)cos(theta2)*(double)data[m][n];
						}
					}

					double ck = (k) ? 0.5 : sqrt((double)1.0 / (double)M);
					double cl = (l) ? 0.5 : sqrt((double)1.0 / (double)N);
					dct[u][v] = ck * cl *sum;
				}
			}
		}
	}
}

void IDCT_2D(double** dct, uchar** result_data) {
	int M = 8;
	int N = 8;
	int mcrNb = Datalength / N;

	double sum = 0;

	int u, v, m, n = 0;
	for (int mcr1 = 0; mcr1 < mcrNb; mcr1++) {
		for (int mcr2 = 0; mcr2 < mcrNb; mcr2++) {
			for (int i = 0; i < M; i++) {
				for (int j = 0; j < N; j++) {
					m = mcr1 * M + i;
					n = mcr2 * N + j;
					sum = 0;
					for (int k = 0; k < M; k++) {
						for (int l = 0; l < N; l++) {
							u = mcr1 * M + k;
							v = mcr2 * N + l;
							double theta1 = (double)(2.0*i + 1)*k*PHI / (double)(2.0*M);
							double theta2 = (double)(2.0*j + 1)*l*PHI / (double)(2.0*N);
							double ck = (k) ? 0.5 : sqrt((double)1.0 / (double)M);
							double cl = (l) ? 0.5 : sqrt((double)1.0 / (double)N);
							sum += ck * cl *(double)cos(theta1)*(double)cos(theta2)*dct[u][v];
						}
					}

					
					result_data[m][n] = (uchar)sum;
				}
			}
		}
	}
}

void make_zero1(double** dct) {
	int M = 8;
	int N = 8;
	int mcrNb = Datalength / N;
	int u, v;
	for (int mcr1 = 0; mcr1 < mcrNb; mcr1++) {
		for (int mcr2 = 0; mcr2 < mcrNb; mcr2++) {
			for (int k = 0; k < M; k++) {
				for (int l = 0; l < N; l++) {
					u = mcr1 * M + k;
					v = mcr2 * N + l;
					if (k > 4 || l > 4) 
						dct[u][v] = 0;						
				}
			}
		}
	}
}

void make_zero2(double** dct) {
	int M = 8;
	int N = 8;
	int mcrNb = Datalength / N;
	int u, v;
	for (int mcr1 = 0; mcr1 < mcrNb; mcr1++) {
		for (int mcr2 = 0; mcr2 < mcrNb; mcr2++) {
			for (int k = 0; k < M; k++) {
				for (int l = 0; l < N; l++) {
					u = mcr1 * M + k;
					v = mcr2 * N + l;
					if (k > 1 && l > 1)
						dct[u][v] = 0;
				}
			}
		}
	}

}

void make_zero3(double** dct) {
	int M = 8;
	int N = 8;
	int mcrNb = Datalength / N;
	int u, v;
	for (int mcr1 = 0; mcr1 < mcrNb; mcr1++) {
		for (int mcr2 = 0; mcr2 < mcrNb; mcr2++) {
			for (int k = 0; k < M; k++) {
				for (int l = 0; l < N; l++) {
					u = mcr1 * M + k;
					v = mcr2 * N + l;
					if (k < 3 && l < 3)
						dct[u][v] = 0;
				}
			}
		}
	}

}
