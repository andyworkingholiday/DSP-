#include <iostream>
#include <cmath>
#include <fstream>
#include "complex.h"
using namespace std;
#define PHI 3.141592
#define FILESIZE 32000
#define SMPL 16000

void Dirichlet(int L, int fx, ofstream& out_diri, complex* Hf);

int main() {
	int fs = 16000;
	int L = 8;
	complex *DFT = new complex[fs];
	complex *Hf = new complex[fs];
	complex *FilteredDFT = new complex[fs];
	
	ifstream InputMusic("noiseMusic.wav", ios::binary);
	ofstream out_diri;
	out_diri.open("diri.txt");

	char*header = new char[44];
	short* data = new short[SMPL];

	InputMusic.read((char*)header, 44);
	InputMusic.read((char*)data, FILESIZE);

	cout << "Start DFT" << endl;
	for (int k = 0; k < fs; k++) {
		for (int n = 0; n < SMPL; n++) {
			DFT[k] += complex(data[n], 0)*complex(-2 * PHI*k*n / (double)fs);
		}
	}

	ofstream Freq("freq.txt");
	for (int k = 0; k < fs; k++)
		Freq << k << "\t" << DFT[k].mag() << endl;

	cout << "Filtering" << endl;
	Dirichlet(L, fs, out_diri, Hf);

	for (int i = 0; i < fs; i++) {
		FilteredDFT[i] = DFT[i] * Hf[i];
	}

	ofstream Filtered_Freq("Filtered_freq.txt");
	for (int k = 0; k < fs; k++)
		Filtered_Freq << k << "\t" << FilteredDFT[k].mag() << endl;

	

	complex* IDFT = new complex[SMPL];
	short* filtered_data = new short[SMPL];
	cout << "Start IDFT" << endl;
	for (int n = 0; n < SMPL; n++) {
		for (int k = 0; k < fs; k++) {
			IDFT[n] += FilteredDFT[k] * complex(2 * PHI*k*n / (double)fs);
		}
		IDFT[n] = IDFT[n] / (double)fs;
	}

	for (int t = 0; t < 16000; t++)
		filtered_data[t] = IDFT[t].re;

	ofstream musicout;
	musicout.open("new_music.wav", ios::binary);
	musicout.write((char*)header, 44);
	musicout.write((char*)filtered_data, sizeof(short) * 16000);
	system("pause");
	



	
}




void Dirichlet(int L, int fx, ofstream& out_diri, complex* H) {
	complex upper, bottom;
	double lim = 8.0;

	for (int k = 0; k < fx; k++)
	{
		bottom = complex(sin(2. * PHI * k / (double)(2. * fx)), 0.);
		if (bottom.mag() == 0.0) {
			H[k] = complex(cos(-2. * PHI * k * ((L - 1) / 2.) / (double)fx), sin(-2. * PHI * k * ((L - 1) / 2.) / (double)fx)) * lim;
		}
		else {
			upper = complex(sin(2. * PHI * k * L / (double)(2. * fx)), 0.0);
			H[k] = upper / bottom * complex(cos((double)(-2. * PHI * k * (L - 1) / 2.) / (double)fx), sin(-2. * PHI * k * ((L - 1) / 2.) / (double)fx));
		}
		H[k].re = H[k].re*0.3;
		H[k].im = H[k].im*0.3;
	}

	

	for (int k = 0; k < fx; k++) {
		out_diri << k << "\t" << H[k].mag() << endl;
	}
}



