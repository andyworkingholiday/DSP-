#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "complex.h"
using namespace std;
#define PHI 3.141592
#define W 256
#define H 256

#define W1 252
#define H1 252
#define uchar unsigned char
#define Datalength 256
#define WORD unsigned short
#define DWORD unsigned int

void FFT2Radix(double* Xr, double* Xi, double* Yr, double* Yi, int nN, int binverse);
void FFT2D(uchar** img, double** OutputReal, double** OutputImag, int nw, int nH);
void FFT2D_inverse(double** InputReal, double** InputImag, double** OutputDouble, int nW, int nH);
void DNormalize2D(double **p1, uchar **p2, int nW, int nH);

int main() {

	//노이즈 있는 252크기 파일//
	ifstream infile;
	infile.open("twins_noise.bmp", ios::binary);
	//256 크기 파일//
	ifstream infile2;
	infile2.open("twins_256.bmp", ios::binary);
	//노이즈 없는 깨끗한 252파일 (마스킹할 때 사용)
	ifstream infile3;
	infile3.open("twins_.bmp", ios::binary);
	char* header;
	char* header2;


	uchar**RGB;
	uchar**R;
	uchar**RGB256;
	uchar**R256;
	uchar**R252_clean;
	uchar**RGB252_clean;
	complex**fft;
	complex**temp2;
	double**mag;
	double**fftRe;
	double**fftIm;
	uchar**R_;
	uchar**RGB_;
	uchar**R_Freq_sort;
	uchar**RGB_Freq_sort;


	header2 = new char[54];
	header = new char[54];
	RGB = new uchar*[H1];
	R = new uchar*[H1];
	R252_clean = new uchar*[H1];
	RGB252_clean = new uchar*[H1];
	R256 = new uchar*[H];
	RGB256 = new uchar*[H];
	fft = new complex*[H];
	temp2 = new complex*[H];

	mag = new double*[H];
	fftRe = new double*[H];
	fftIm = new double*[H];
	R_ = new uchar*[H];
	RGB_ = new uchar*[H];
	R_Freq_sort = new uchar*[H];
	RGB_Freq_sort = new uchar*[H];

	for (int i = 0; i < H1; i++) {
		RGB[i] = new uchar[3 * W1];
		R[i] = new uchar[W1];
		R252_clean[i] = new uchar[W1];
		RGB252_clean[i] = new uchar[3 * W1];
	}

	for (int i = 0; i < H; i++) {
		RGB256[i] = new uchar[3 * W];
		R256[i] = new uchar[W];
		fft[i] = new complex[W];
		temp2[i] = new complex[W];
		mag[i] = new double[W];
		fftRe[i] = new double[W];
		fftIm[i] = new double[W];
		R_[i] = new uchar[W];
		RGB_[i] = new uchar[W * 3];
		R_Freq_sort[i] = new uchar[W];
		RGB_Freq_sort[i] = new uchar[W * 3];
	}

	infile.read((char*)header, 54);
	infile2.read((char*)header2, 54);
	infile3.read((char*)header, 54);

	//헤더정보 출력 252 크기 일 때, 256 크기 일 때 각각
	cout << "매직넘버" << "\t" << header[0] << header[1] << endl;
	cout << "BMP파일 전체크기" << "\t" << *(DWORD*)(header + 2) << endl;
	cout << "0FFSET" << "\t" << *(DWORD*)(header + 10) << endl;

	cout << "가로화소" << "\t" << *(DWORD*)(header + 18) << endl;
	cout << "세로화소" << "\t" << *(DWORD*)(header + 22) << endl;
	cout << "데이터 총 크기" << "\t" << *(DWORD*)(header + 34) << endl;
	cout << "가로 해상도 " << "\t" << *(DWORD*)(header + 38) << endl;
	cout << "세로 해상도" << "\t" << *(DWORD*)(header + 42) << endl;
	cout << endl << endl;

	cout << "매직넘버" << "\t" << header2[0] << header2[1] << endl;
	cout << "BMP파일 전체크기" << "\t" << *(DWORD*)(header2 + 2) << endl;
	cout << "0FFSET" << "\t" << *(DWORD*)(header2 + 10) << endl;

	cout << "가로화소" << "\t" << *(DWORD*)(header2 + 18) << endl;
	cout << "세로화소" << "\t" << *(DWORD*)(header2 + 22) << endl;
	cout << "데이터 총 크기" << "\t" << *(DWORD*)(header2 + 34) << endl;
	cout << "가로 해상도 " << "\t" << *(DWORD*)(header2 + 38) << endl;
	cout << "세로 해상도" << "\t" << *(DWORD*)(header2 + 42) << endl;

	//노이즈 있는 파일 불러들여오기
	for (int i = 0; i < H1; i++)
		infile.read((char*)RGB[i], 3 * W1);

	for (int i = 0; i < H1; i++) {
		for (int j = 0, jj = 0; j < W1; j++, jj += 3) {
			R[i][j] = RGB[i][jj];
		}
	}

	//노이즈 없는 파일 불러들여오기
	for (int i = 0; i < H1; i++)
		infile3.read((char*)RGB252_clean[i], 3 * W1);

	for (int i = 0; i < H1; i++) {
		for (int j = 0, jj = 0; j < W1; j++, jj += 3) {
			R252_clean[i][j] = RGB252_clean[i][jj];
		}
	}

	//초기화
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			R256[i][j] = 0;
		}
	}

	//FFT를 하기 위해 252크기를 256크기 배열로 옮겨준다
	for (int i = 0; i < H1; i++) {
		for (int j = 0; j < W1; j++) {
			R256[i][j] = R[i][j];
		}
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			RGB256[i][jj] = R256[i][j];
			RGB256[i][jj + 1] = R256[i][j];
			RGB256[i][jj + 2] = R256[i][j];

		}
	}

	
	//FFT 실행//
	FFT2D(R256, fftRe, fftIm, W, H);


	
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
			RGB_[i][jj + 1] = R_[i][j];
			RGB_[i][jj + 2] = R_[i][j];
		}
	}

	// 저주파, 고주파수 대역 보기 위해 출력
	ofstream Outfile;
	Outfile.open("result.bmp", ios::binary);
	Outfile.write((char*)header2, 54);
	for (int i = 0; i < H; i++)
		Outfile.write((char*)RGB_[i], 3 * W);
	Outfile.close();


	//주파수 정렬
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
	Outfile2.write((char*)header2, 54);
	for (int i = 0; i < H; i++)
		Outfile2.write((char*)RGB_Freq_sort[i], 3 * W);
	Outfile2.close();

	//노이즈 확인 후 제거


	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			if
				((i - 156)*(i - 156) + (j - 86)*(j - 86) < 30)
				R_[i][j] = 0;
			else if
				((i - 100)*(i - 100) + (j - 170)*(j - 170) < 30)
				R_[i][j] = 0;

			//4사분면 세로
			for (int k = 80; k < 126; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 170)*(j - 170) < 2)
					R_[i][j] = 0;

			}

			// 4사분면 세로 맨밑

			for (int k = 30; k < 70; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 170)*(j - 170) < 2)
					R_[i][j] = 0;

			}
			//2사분면 가로
			for (int k = 170; k < 200; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 170)*(j - 170) < 1)
					R_[i][j] = 0;

			}

			// 2사분면 세로
			for (int k = 129; k < 256; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 86)*(j - 86) < 2)
					R_[i][j] = 0;

			}
			// 4사분면 가로
			for (int k = 150; k < 200; k += 1)
			{
				if
					((i - 100)*(i - 100) + (j - k)*(j - k) < 2)
					R_[i][j] = 0;

			}
			// 4사분면 가로 맨오른쪽
			for (int k = 227; k < 256; k += 1)
			{
				if
					((i - 100)*(i - 100) + (j - k)*(j - k) < 2)
					R_[i][j] = 0;

			}
			// 1사분면 가로
			for (int k = 60; k < 126; k += 1)
			{
				if
					((i - 156)*(i - 156) + (j - k)*(j - k) < 2)
					R_[i][j] = 0;
			}
			// 1사분면 세로 원쪽
			for (int k = 130; k < 17; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 86)*(j - 86) < 2)
					R_[i][j] = 0;

			}

			// 3사분면 가로 왼쪽
			for (int k = 20; k < 50; k += 1)
			{
				if
					((i - 100)*(i - 100) + (j - k)*(j - k) < 2)
					R_[i][j] = 0;

			}

			// 3사분면 가로 오른쪽
			for (int k = 80; k < 126; k += 1)
			{
				if
					((i - 100)*(i - 100) + (j - k)*(j - k) < 2)
					R_[i][j] = 0;

			}

			//3사분면 세로
			for (int k = 0; k < 126; k += 1)
			{
				if
					((i - k)*(i - k) + (j - 86)*(j - 86) < 2)
					R_[i][j] = 0;

			}


		}

	}

	for (int i = 0; i < H; i++) {
		for (int j = 0, jj = 0; j < W; j++, jj += 3) {
			RGB_[i][jj] = R_[i][j];
			RGB_[i][jj + 1] = R_[i][j];
			RGB_[i][jj + 2] = R_[i][j];
		}
	}

	ofstream eli_freq;
	eli_freq.open("noise_eliminated_freq.bmp", ios::binary);
	eli_freq.write((char*)header2, 54);

	for (int i = 0; i < H; i++)
		eli_freq.write((char*)RGB_[i], 3 * W);
	eli_freq.close();

	// FFT Centralization -> Noise 좌표 알아낸대로 필터링하기 위하여
  //////////////////////////////////////////////////////
	for (int i = 0; i < H / 2; i++)
	{
		for (int j = 0; j < W / 2; j++) {

			temp2[i][j] = fft[i][j];
			fft[i][j] = fft[i + (H / 2)][j + (W / 2)];
			fft[i + (H / 2)][j + (W / 2)] = temp2[i][j];
		}
	}

	for (int i = 0; i < H / 2; i++)
	{
		for (int j = W / 2; j < W; j++)
		{
			temp2[i][j] = fft[i][j];
			fft[i][j] = fft[i + (H / 2)][j - (W / 2)];
			fft[i + (H / 2)][j - (W / 2)] = temp2[i][j];
		}
	}

	//FFT noise 필터링
	//////////////////////////////////////////////////////

	for (int i = 0; i < H; i++)
	{
		for (int j = 0, jj = 0; j < W; j++, jj += 3)
		{
			for (int k = 0; k < 127; k += 1)
			{
				if
					((i - 100) * (i - 100) + (j - k) * (j - k) < 2)
					fft[i][j] = complex(0, 0);
			}
			for (int k = 130; k < 256; k += 1)
			{
				if
					((i - 100) * (i - 100) + (j - k) * (j - k) < 2)
					fft[i][j] = complex(0, 0);
			}

			for (int k = 0; k < 126; k += 1)
			{
				if
					((i - k) * (i - k) + (j - 170) * (j - 170) < 2)
					fft[i][j] = complex(0, 0);
			}
			for (int k = 130; k < 256; k += 1)
			{
				if
					((i - k) * (i - k) + (j - 170) * (j - 170) < 2)
					fft[i][j] = complex(0, 0);
			}

			for (int k = 0; k < 126; k += 1)
			{
				if
					((i - 156) * (i - 156) + (j - k) * (j - k) < 2)
					fft[i][j] = complex(0, 0);
			}
			for (int k = 130; k < 256; k += 1)
			{
				if
					((i - 156) * (i - 156) + (j - k) * (j - k) < 2)
					fft[i][j] = complex(0, 0);
			}

			for (int k = 0; k < 126; k += 1)
			{
				if
					((i - k) * (i - k) + (j - 86) * (j - 86) < 2)
					fft[i][j] = complex(0, 0);
			}
			for (int k = 130; k < 256; k += 1)
			{
				if
					((i - k) * (i - k) + (j - 86) * (j - 86) < 2)
					fft[i][j] = complex(0, 0);
			}

			if
				((i - 100) * (i - 100) + (j - 170) * (j - 170) < 25)
				fft[i][j] = complex(0, 0);
			// 오른쪽
			else if
				((i - 156) * (i - 156) + (j - 86) * (j - 86) < 25)
				fft[i][j] = complex(0, 0);
		}
	}
	//////////////////////////////////////////////////////


	// FFT Centralization 재정렬 -> 원래대로 돌려놓기

	for (int i = 0; i < H / 2; i++)
	{
		for (int j = 0; j < W / 2; j++) {

			temp2[i][j] = fft[i][j];
			fft[i][j] = fft[i + (H / 2)][j + (W / 2)];
			fft[i + (H / 2)][j + (W / 2)] = temp2[i][j];
		}
	}

	for (int i = 0; i < H / 2; i++)
	{
		for (int j = W / 2; j < W; j++)
		{
			temp2[i][j] = fft[i][j];
			fft[i][j] = fft[i + (H / 2)][j - (W / 2)];
			fft[i + (H / 2)][j - (W / 2)] = temp2[i][j];
		}
	}
	//////////////////////////////////////////////////////


	 //fftRe, fftIm 초기화
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			fftRe[i][j] = fft[i][j].re;
			fftIm[i][j] = fft[i][j].im;
		}
	}


	
	

	//Inverse DFT 진행//
	double** result_data;
	result_data = new double*[H];
	for (int i = 0; i < H; i++)
		result_data[i] = new double[W];

	FFT2D_inverse(fftRe, fftIm, result_data, W, H);


	
	double** result_data_252;
	result_data_252 = new double*[H1];
	for (int i = 0; i < H1; i++)
		result_data_252[i] = new double[W1];


	for (int i = 0; i < H1; i++) {
		for (int j = 0; j < W1; j++) {
			result_data_252[i][j] = result_data[i][j];
		}
	}

	uchar** noise_eli;
	noise_eli = new uchar*[H1];
	for (int i = 0; i < H1; i++)
		noise_eli[i] = new uchar[3 * W1];


	for (int i = 0; i < H1; i++) {
		for (int j = 0, jj = 0; j < W1; j++, jj += 3) {
			noise_eli[i][jj] = result_data_252[i][j];
			noise_eli[i][jj + 1] = result_data_252[i][j];
			noise_eli[i][jj + 2] = result_data_252[i][j];

		}
	}


	ofstream fftinverse;
	fftinverse.open("noise_elminated_but_notmask.bmp", ios::binary);
	fftinverse.write((char*)header, 54);


	for (int i = 0; i < H1; i++)
		fftinverse.write((char*)noise_eli[i], 3 * W1);





	//마스킹 작업 시작//
	uchar* Orig;
	Orig = new uchar[252 * 252];

	for (int i = 0; i < H1; i++) {
		for (int j = 0; j < W1; j++) {
			Orig[i * 252 + j] = R252_clean[i][j];
		}
	}

	double* temp;
	temp = new double[252 * 252];

	//초기화
	for (int i = 0; i < 252 * 252; i++) {
		temp[i] = 0;
	}

	int m, n;
	double norm = 0;

	double kernel[3][3] = { {-1, -1, -1},
		  {-1,  9, -1},
		  {-1, -1, -1} };

	for (m = 0; m < 3; m++) {
		for (n = 0; n < 3; n++) {
			norm += kernel[m][n];
		}
	}
	if (norm == 0)
		norm = 1;
	else
		norm = 1 / norm;


	for (int row = 1; row < H1 - 1; row++) {
		for (int col = 1; col < W1 - 1; col++) {
			for (m = 0; m < 3; m++) {
				for (n = 0; n < 3; n++) {
					temp[row*252 + col] += (kernel[m][n] * (double)Orig[(row - 1 + m)*252 + (col - 1 + n)]) * norm;
				}
			}
		}
	}

	// 마스킹 작업 후 값이 정해지지 않은 외곽부분 메우기
	for (int row = 0; row < H1; row++) {
		for (int col = 0; col < W1; col++) {
			switch (row) {
			case 0:
				temp[row*W1 + col] = temp[(row + 1)*W1 + col];
				break;
			case 251:
				temp[row*W1 + col] = temp[(row - 1)*W1 + col];
				break;
			}
		}
	}
	for (int row = 0; row < H1; row++) {
		for (int col = 0; col < W1; col++) {
			switch (col) {
			case 0:
				temp[row*W1 + col] = temp[row*W1 + col + 1];
				break;
			case 251:
				temp[row*W1 + col] = temp[row*W1 + col - 1];
				break;
			}
		}
	}

	//윤곽검출//
	for (int Row = 0; Row < H1; Row++) {
		for (int Col = 0; Col < W1; Col++) {
			temp[Row*W1 + Col] = temp[Row*W1 + Col] - Orig[Row*W1 + Col];
		}
	}

	//최대최소//
	for (int row = 0; row < H1; row++) {
		for (int col = 0; col < W1; col++) {
			if (temp[row*W1 + col] < 0) temp[row*W1 + col] = 0;
			else if (temp[row*W1 + col] > 255) temp[row*W1 + col] = 255;
		}
	}

	//마스킹 출력//
	uchar** masking;
	masking = new uchar*[H1];
	for (int i = 0; i < H1; i++)
		masking[i] = new uchar[3 * W1];

	for (int i = 0; i < H1; i++) {
		for (int j = 0, jj = 0; j < W1; j++, jj += 3) {
			masking[i][jj] = temp[252 * i + j];
			masking[i][jj + 1] = temp[252 * i + j];
			masking[i][jj + 2] = temp[252 * i + j];

		}
	}

	ofstream mask;
	mask.open("masking.bmp", ios::binary);
	mask.write((char*)header, 54);


	for (int i = 0; i < H1; i++)
		mask.write((char*)masking[i], 3 * W1);

	//noise 제거 된 파일과 합치기//
	for (int i = 0; i < H1; i++) {
		for (int j = 0; j < W1; j++) {
			result_data_252[i][j] = result_data_252[i][j] + 0.3*temp[252 * i + j];
		}
	}

	for (int row = 0; row < H1; row++) {
		for (int col = 0; col < W1; col++) {
			if (result_data_252[row][col] < 0) result_data_252[row][col] = 0;
			else if (result_data_252[row][col] > 255) result_data_252[row][col] = 255;
		}
	}

	uchar** end_output;
	end_output= new uchar*[H1];
	for (int i = 0; i < H1; i++)
		end_output[i] = new uchar[3 * W1];
	

	for (int i = 0; i < H1; i++) {
		for (int j = 0, jj = 0; j < W1; j++, jj += 3) {
			end_output[i][jj] = result_data_252[i][j];
			end_output[i][jj + 1] = result_data_252[i][j];
			end_output[i][jj + 2] = result_data_252[i][j];

		}
	}
	


	//최종출력//
	ofstream fftout;
	fftout.open("noise_elminated_output.bmp", ios::binary);
	fftout.write((char*)header, 54);


	for (int i = 0; i < H1; i++)
		fftout.write((char*)end_output[i], 3 * W1);

	

	
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
	for (y = 0; y < nH; y++)
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











