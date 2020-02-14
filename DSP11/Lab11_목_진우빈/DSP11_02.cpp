#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "complex.h"
using namespace std;
#define PHI 3.141592
#define FILESIZE 32000
#define SMPL 16000
#define N 8000

int main() {
	int fs = 16000;
	int L = 8;
	complex *DFT = new complex[fs];
	complex *FilteredDFT = new complex[fs];

	ifstream InputMusic("noiseMusic.wav", ios::binary);

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
	
	
	double polestheta1, polestheta2, polestheta3, polestheta4, polestheta5, polestheta6;
	double zerotheta1, zerotheta2, zerotheta3, zerotheta4, zerotheta5, zerotheta6;
	polestheta1 = 6 * PHI / 16;
	polestheta2 = 8 * PHI / 16;
	polestheta3 = 10 * PHI / 16;
	polestheta4 = 22 * PHI / 16;
	polestheta5 = 24 * PHI / 16;
	polestheta6 = 26 * PHI / 16;

	zerotheta1 = 3 * PHI / 8;
	zerotheta2 = 0.;
	zerotheta3 = 5 * PHI / 8;
	zerotheta4 = 11 * PHI / 8;
	zerotheta5 = PHI;
	zerotheta6 = 13 * PHI / 8;

	complex pole1, pole2, pole3, pole4, pole5, pole6;
	complex zero1, zero2, zero3, zero4, zero5, zero6;

	pole1 = complex(polestheta1);
	pole2 = complex(polestheta2);
	pole3 = complex(polestheta3);
	pole4 = complex(polestheta4);
	pole5 = complex(polestheta5);
	pole6 = complex(polestheta6);

	zero1 = complex(zerotheta1);
	zero2 = complex(zerotheta2);
	zero3 = complex(zerotheta3);
	zero4 = complex(zerotheta4);
	zero5 = complex(zerotheta5);
	zero6 = complex(zerotheta6);

	complex* H = new complex[16000];
	complex* Z = new complex[16000];

	ofstream out;
	out.open("zeropole.txt");


	for (int k = 0; k < 16000; k++) {
		Z[k] = complex(2 * PHI * k / (double)16000);
		double a1 = Z[k].re - pole1.re;
		double b1 = Z[k].im - pole1.im;
		double a2 = Z[k].re - pole2.re;
		double b2 = Z[k].im - pole2.im;
		double a3 = Z[k].re - pole3.re;
		double b3 = Z[k].im - pole3.im;
		double a4 = Z[k].re - pole4.re;
		double b4 = Z[k].im - pole4.im;
		double a5 = Z[k].re - pole5.re;
		double b5 = Z[k].im - pole5.im;
		double a6 = Z[k].re - pole6.re;
		double b6 = Z[k].im - pole6.im;

		H[k] = (Z[k] - zero1) * (Z[k] - zero2) * (Z[k] - zero3)* (Z[k] - zero4)*(Z[k] - zero5)*(Z[k] - zero6)*
			complex(a1, -b1)* complex(a2, -b2)*complex(a3, -b3)* complex(a4, -b4)*complex(a5, -b5)* complex(a6, -b6);
		H[k] = H[k] / ((a1*a1 + b1 * b1)*(a2*a2 + b2 * b2)* (a3*a3 + b3 * b3)*(a4*a4 + b4 * b4)* (a5*a5 + b5 * b5)*(a6*a6 + b6 * b6));
		H[k].re = abs(log(abs(H[k].re) + 1)) / 3;
		H[k].im = abs(log(abs(H[k].im) + 1)) / 3;
	
		

	}

	for (int k = 0; k < 16000; k++) {
		if (k == 4000)
			H[k] = complex(0, 3);
		else if (k == 12000)
			H[k] = complex(0, -3);

		else if (k == 11000)
			H[k] = complex((1.22795/3)*cos(11 * PHI / 8), (1.22795/3)*sin(11 * PHI / 8));

		else if (k == 13000)
			H[k] = complex((1.22795/3)*cos(13 * PHI / 8), (1.22795/3)*sin(13 * PHI / 8));


		out << k << "\t" << H[k].mag() << endl;
	}
 	
	

	cout << "Filtering" << endl;


	for (int i = 0; i < fs; i++) {
		FilteredDFT[i] = DFT[i] * H[i];
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
	musicout.open("new_music2.wav", ios::binary);
	musicout.write((char*)header, 44);
	musicout.write((char*)filtered_data, sizeof(short) * 16000);
	system("pause");


	
	
}







