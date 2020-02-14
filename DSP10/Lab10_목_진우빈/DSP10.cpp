#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "complex.h"
using namespace std;
#define PHI 3.141592
#define W 256
#define H 256
#define uchar unsigned char
#define Datalength 256



void FFT2Radix(double* Xr, double* Xi, double* Yr, double* Yi, int nN, int binverse);
void FFT2D(uchar** img, double** OutputReal, double** OutputImag, int nw, int nH);
void FFT2D_inverse(double** InputReal, double** InputImag, double** OutputDouble, int nW, int nH);
void DNormalize2D(double **p1, uchar **p2, int nW, int nH);


int main() {
	ifstream infile;
	infile.open("twin_noise.bmp", ios::binary);
	char* header;
	uchar**RGB;
	uchar**R;
	complex**fft;
	double**mag;
	double**fftRe;
	double**fftIm;
	uchar**R_;
	uchar**RGB_;
	uchar**R_Freq_sort;
	uchar**RGB_Freq_sort;
	

	header = new char[54];
	RGB = new uchar*[H];
	R = new uchar*[H];
	fft = new complex*[H];
	mag = new double*[H];
	fftRe = new double*[H];
	fftIm = new double*[H];
	R_ = new uchar*[H];
	RGB_ = new uchar*[H];
	R_Freq_sort = new uchar*[H];
	RGB_Freq_sort = new uchar*[H];

	for (int i = 0; i < H; i++) {
		RGB[i] = new uchar[3 * W];
		R[i] = new uchar[W];
		fft[i] = new complex[W];
		mag[i] = new double[W];
		fftRe[i] = new double[W];
		fftIm[i] = new double[W];
		R_[i] = new uchar[W];
		RGB_[i] = new uchar[W * 3];
		R_Freq_sort[i] = new uchar[W];
		RGB_Freq_sort[i] = new uchar[W * 3];
	}

	infile.read((char*)header, 54);
	for (int i = 0; i < H; i++)
		infile.read((char*)RGB[i], 3 * W);

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			R[i][j] = RGB[i][jj];
		}
	}

	FFT2D(R, fftRe, fftIm, W, H);

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			fft[i][j] = complex(fftRe[i][j], fftIm[i][j]);
			mag[i][j] = 10 * log(fft[i][j].mag() + 1);
		}
	}

	DNormalize2D(mag, R_, W, H);

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			RGB_[i][jj] = R_[i][j];
			RGB_[i][jj+1] = R_[i][j];
			RGB_[i][jj+2] = R_[i][j];
		}
	}
	
	ofstream Outfile;
	Outfile.open("result.bmp", ios::binary);
	Outfile.write((char*)header, 54);
	for (int i = 0; i < H; i++)
		Outfile.write((char*)RGB_[i], 3 * W);
	Outfile.close();
	
	//sortFreq(R_Freq_sort, R_);

	const int center_i = 128;
	const int center_j = 128;

	for (int i = center_i; i < W; i++) {
		for (int j = center_j; j < H; j++) {
			R_Freq_sort[i - center_i][j - center_j] = R_[i][j];
			//4사분면을 2사분면으로

		}
	}

	for (int i = 0; i < center_i; i++) {
		for (int j = 0; j < center_j; j++) {
			R_Freq_sort[i + center_i][j + center_j] = R_[i][j];
			//2사분면을 4사분면으로

		}
	}

	for (int i = 0; i < center_i; i++) {
		for (int j = center_j; j < H; j++) {
			R_Freq_sort[i + center_i][j - center_j] = R_[i][j];
			//1사분면을 3사분면으로

		}
	}

	for (int i = center_i; i < W; i++) {
		for (int j = 0; j < center_j; j++) {
			R_Freq_sort[i - center_i][j + center_j] = R_[i][j];
			//3사분면을 1사분면으로

		}
	}


	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			RGB_Freq_sort[i][jj] = R_Freq_sort[i][j];
			RGB_Freq_sort[i][jj + 1] = R_Freq_sort[i][j];
			RGB_Freq_sort[i][jj + 2] = R_Freq_sort[i][j];
		}
	}

	//주파수 정렬 후 출력//
	ofstream Outfile2;
	Outfile2.open("Freq_sort_result.bmp", ios::binary);
	Outfile2.write((char*)header, 54);
	for (int i = 0; i < H; i++)
		Outfile2.write((char*)RGB_Freq_sort[i], 3 * W);
	Outfile2.close();

	//노이즈 확인 후 제거//

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (sqrt((i - 229)*(i - 229) + (j - 41)*(j - 41)) < 9)
				R_[i][j] = 0;

			else if (sqrt((i - 27)*(i - 27) + (j - 215)*(j - 215)) < 9)
				R_[i][j] = 0;
		
		}
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			RGB_[i][jj] = R_[i][j];
			RGB_[i][jj + 1] = R_[i][j];
			RGB_[i][jj + 2] = R_[i][j];
		}
	}

	ofstream Outfile3;
	Outfile3.open("noise_eliminated_result.bmp", ios::binary);
	Outfile3.write((char*)header, 54);


	for (int i = 0; i < H; i++)
		Outfile3.write((char*)RGB_[i], 3 * W);
	Outfile3.close();
	
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if (sqrt((i - 229)*(i - 229) + (j - 41)*(j - 41)) < 9) {
				fftRe[i][j] = 0;
				fftIm[i][j] = 0;
			}
				
				

			else if (sqrt((i - 27)*(i - 27) + (j - 215)*(j - 215)) < 9) {
				fftRe[i][j] = 0;
				fftIm[i][j] = 0;
			}

		}
	}



	/*for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			fft2[i][j] = complex(fftRe2[i][j], fftIm2[i][j]);
			mag2[i][j] = pow(10, (fft2[i][j].mag() / 10)) - 1;
			
		}
	}*/


	double** result_data;
	result_data = new double*[H];
	for (int i = 0; i < H; i++)
		result_data[i] = new double[W];

	FFT2D_inverse(fftRe, fftIm, result_data, W, H);

	uchar** output_data;
	output_data = new uchar*[H];
	for (int i = 0; i < H; i++)
		output_data[i] = new uchar[3 * W];

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			output_data[i][jj] = (uchar)result_data[i][j];
			output_data[i][jj + 1] = (uchar)result_data[i][j];
			output_data[i][jj + 2] = (uchar)result_data[i][j];

		}
	}

	


	ofstream fftout;
	fftout.open("fft_output_twin.bmp", ios::binary);
	fftout.write((char*)header, 54);


	for (int i = 0; i < H; i++)
		fftout.write((char*)output_data[i], 3 * W);

	system("pause");
	return 0;


}

void FFT2Radix(double* Xr, double* Xi, double* Yr, double* Yi, int nN, int binverse) {
	double T, Wr, Wi;

	for (int i = 0; i < nN; i++) {
		Yr[i] = Xr[i];
		Yi[i] = Xi[i];
	}
	int j = 0, k = 0;
	for (int i = 1; i < (nN - 1); i++) {
		k = nN / 2;
		while (k <= j) {
			j = j - k;
			k = k / 2;
		}
		j = j + k;
		if (i < j) {
			T = Yr[j];
			Yr[j] = Yr[i];
			Yr[i] = T;

			T = Yi[j];
			Yi[j] = Yi[i];
			Yi[i] = T;
		}
	}

	double Tr, Ti;
	int iter, j2, pos;
	k = nN >> 1;
	iter = 1;
	while (k > 0) {
		j = 0;
		j2 = 0;
		for (int i = 0; i < nN >> 1; i++) {
			Wr = cos(2.*PHI*(j2*k) / nN);
			if (binverse == 0)
				Wi = -sin(2.*PHI*(j2*k) / nN);
			else
				Wi = sin(2.*PHI*(j2*k) / nN);
			pos = j + (1 << (iter - 1));

			Tr = Yr[pos] * Wr - Yi[pos] * Wi;
			Ti = Yr[pos] * Wi + Yi[pos] * Wr;

			Yr[pos] = Yr[j] - Tr;
			Yi[pos] = Yi[j] - Ti;
			Yr[j] += Tr;
			Yi[j] += Ti;

			j += 1 << iter;
			if (j >= nN) j = ++j2;
		}
		k >>= 1;
		iter++;

	}
	if (binverse) {
		for (int i = 0; i < nN; i++) {
			Yr[i] /= nN;
			Yi[i] /= nN;
		}
	}
}

void FFT2D(uchar** img, double** OutputReal, double** OutputImag, int nW, int nH) {
	int x, y;
	double *dRealX, *dImagX;
	double *dRealY, *dImagY;

	dRealX = new double[nW];
	dImagX = new double[nW];
	dRealY = new double[nW];
	dImagY = new double[nW];

	for (y = 0; y < nH; y++) {
		for (x = 0; x < nW; x++) {
			dRealX[x] = img[y][x];
			dImagX[x] = 0.0;
		}

		FFT2Radix(dRealX, dImagX, dRealY, dImagY, nW, false);

		for (x = 0; x < nW; x++) {
			OutputReal[y][x] = dRealY[x];
			OutputImag[y][x] = dImagY[x];
		}
	}

	delete[] dRealX;
	delete[] dImagX;
	delete[] dRealY;
	delete[] dImagY;

	dRealX = new double[nH];
	dImagX = new double[nH];
	dRealY = new double[nH];
	dImagY = new double[nH];

	for (x = 0; x < nW; x++) {
		for (y = 0; y < nH; y++) {
			dRealX[y] = OutputReal[y][x];
			dImagX[y] = OutputImag[y][x];
		}

		FFT2Radix(dRealX, dImagX, dRealY, dImagY, nH, false);

		for (y = 0; y < nH; y++) {
			OutputReal[y][x] = dRealY[y];
			OutputImag[y][x] = dImagY[y];
		}
	}

	delete[] dRealX;
	delete[] dImagX;
	delete[] dRealY;
	delete[] dImagY;
}

void FFT2D_inverse(double** InputReal, double** InputImag, double** OutputDouble, int nW, int nH) {
	int x, y;
	double *dRealX, *dImagX;
	double *dRealY, *dImagY;
	double ** OutputReal, **OutputImag;

	OutputReal = new double*[nH];
	OutputImag = new double*[nH];
	for (int i = 0; i < nH; i++) {
		OutputReal[i] = new double[nW];
		OutputImag[i] = new double[nW];
	}

	dRealX = new double[nW];
	dImagX = new double[nW];
	dRealY = new double[nW];
	dImagY = new double[nW];

	for (y = 0; y < nH; y++) {
		for (x = 0; x < nW; x++) {
			dRealX[x] = InputReal[y][x];
			dImagX[x] = InputImag[y][x];
		}

		FFT2Radix(dRealX, dImagX, dRealY, dImagY, nW, true);

		for (x = 0; x < nH; x++) {
			OutputReal[y][x] = dRealY[x];
			OutputImag[y][x] = dImagY[x];
		}
	}

	delete[] dRealX;
	delete[] dImagX;
	delete[] dRealY;
	delete[] dImagY;

	dRealX = new double[nH];
	dImagX = new double[nH];
	dRealY = new double[nH];
	dImagY = new double[nH];

	for (x = 0; x < nW; x++) {
		for (y = 0; y < nH; y++) {
			dRealX[y] = OutputReal[y][x];
			dImagX[y] = OutputImag[y][x];
		}

		FFT2Radix(dRealX, dImagX, dRealY, dImagY, nW, true);

		for (y = 0; y < nH; y++) {
			OutputReal[y][x] = dRealY[y];
			OutputImag[y][x] = dImagY[y];
		}

	}

	delete[] dRealX;
	delete[] dImagX;
	delete[] dRealY;
	delete[] dImagY;

	for (y = 0; y < nH; y++) {
		for (x = 0; x < nW; x++) {
			OutputDouble[y][x] = OutputReal[y][x];
		}
	}

	for (int i = 0; i < nH; i++) {
		delete[] OutputReal[i];
		delete[] OutputImag[i];
	}

	delete[] OutputReal;
	delete[] OutputImag;


}

void DNormalize2D(double **p1, uchar **p2, int nW, int nH) {
	int x, y;
	double min = 9999;
	double max = -9999;
	double val;
	for (y=0; y<nH; y++) 
		for (x = 0; x < nW; x++) {
			val = p1[y][x];
			if (val > max) max = val;
			if (val < min) min = val;
		}

	if (max == min) {
		for (y = 0; y < nH; y++)
			for (x = 0; x < nW; x++)
				p2[y][x] = 0;
		return;
	}

	double dFactor = 255 / (max - min);
	for (y = 0; y < nH; y++)
		for (x = 0; x < nW; x++)
			p2[y][x] = (uchar)((p1[y][x] - min)*dFactor);
}



	


	

	
