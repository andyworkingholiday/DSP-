#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "complex.h"
#include "pch.h"
using namespace std;
#define PHI 3.141592
#define W 64
#define H 64
#define uchar unsigned char

int main() {

	ifstream noise;
	noise.open("twin_noise_64.bmp", ios::binary);

	int headersize = 54;
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

	complex** dft;
	R = new uchar*[H];
	G = new uchar*[H];
	B = new uchar*[H];
	RGB = new uchar*[H]; 

	dft = new complex*[H];

	for (int i = 0; i < H; i++) {
		R[i] = new uchar[W];
		G[i] = new uchar[W];
		B[i] = new uchar[W];
		RGB[i] = new uchar[3*W];
		dft[i] = new complex[W];
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

	for (int v = 0; v < H; v++) {
		for (int u = 0; u < W; u++) {
			for (int y = 0; y < H; y++) {
				for (int x = 0; x < W; x++) {
					dft[v][u] += complex(R[y][x], 0)*complex(-2.*PHI*(((double)u*x / W) + ((double)v*y / H)));
				}
			}
		}
	}
	//흑백일때는 r값만 해도된다.



	ofstream outfile;
	outfile.open("dft.txt");
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			outfile << dft[i][j].mag() << "\t";
		}
		outfile << endl;
	}
	

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if ((i >= 25 && i <= 32) && (j >= 20 && j <= 30))
				dft[i][j] = complex(0, 0);

			if ((i >= 34 && i <= 41) && (j >= 36 && j <= 46))
				dft[i][j] = complex(0, 0);
		}
	}


	complex **result;
	result = new complex*[H];
	for (int i = 0; i < H; i++)
		result[i] = new complex[W];

	for (int y = 0; y < H; y++) {
		for (int x = 0; x < W; x++) {
			for (int v = 0; v < H; v++) {
				for (int u = 0; u < W; u++) {
					/*result[y][x] = result[y][x] + (dft[v][u].re)*cos(2.*PHI*((double)u*x / W + (double)v*y / H));
					result[y][x] = result[y][x] - (dft[v][u].im)*cos(2.*PHI*((double)u*x / W + (double)v*y / H));*/
					result[y][x] += dft[v][u] * complex(2.*PHI*(((double)u*x / W) + ((double)v*y / H)));

				}
			}
			result[y][x] = result[y][x] / (W*H);
		}
	}

	unsigned char** clean_data;
	clean_data = new unsigned char*[H];
	for (int i = 0; i < H; i++)
		clean_data[i] = new unsigned char[3*W];


	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj+=3) {
			clean_data[i][jj] = result[i][j].re;
			clean_data[i][jj+1] = result[i][j].re;
			clean_data[i][jj+2] = result[i][j].re;

		}
	}




	ofstream clean;
	clean.open("clean.bmp", ios::binary);
	clean.write((char*)header, 54);

	for (int i = 0; i < H; i++)
		clean.write((char*)clean_data[i], 3*W);

	system("pause");
	return 0;




 }