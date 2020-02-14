#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "complex.h"
using namespace std;
#define PHI 3.141592
#define WORD unsigned short
#define DWORD unsigned int
int main() {

	char*cMusic = new char[44];
	char*Music = new char[44];
	char*new_header = new char[44];
	ifstream fff, fff2;
	fff.open("musicA.wav", ios::binary);
	fff2.open("MixA.wav", ios::binary);

	fff.read((char*)cMusic, 44);
	fff2.read((char*)Music, 44);

	cout << "RIF:		" << Music[0] << Music[1] << Music[2] << Music[3] << endl;
	cout << "Filesize:	" << *(DWORD*)(Music + 4) << endl;
	cout << "WAVE:		" << Music[8] << Music[9] << Music[10] << Music[11] << endl;
	cout << "cksize		" << *(DWORD*)(Music + 16) << endl;
	cout << "channels	" << *(WORD*)(Music + 22) << endl;
	cout << "fs			" << *(DWORD*)(Music + 24) << endl;
	cout << "bytes/s	" << *(DWORD*)(Music + 28) << endl;
	cout << "bits/sam	" << *(WORD*)(Music + 34) << endl;
	cout << "data		" << Music[36] << Music[37] << Music[38] << Music[39] << endl;
	cout << "cksize		" << *(DWORD*)(Music + 40) << endl;

	int N = 16000;

	short *clean_data = new short[N];
	short *mixed_data = new short[N];

	fff.read((char*)clean_data, N);
	fff2.read((char*)mixed_data, N);


	complex *DFT = new complex[N];
	complex *DFT2 = new complex[N];

	for (int k = 0; k < N; k++) {
		for (int n = 0; n < N; n++) 
			DFT[k] += complex(clean_data[n], 0) * complex(cos((-2 * PHI*k*n) / (double)N), sin((-2 * PHI*k*n) / (double)N));

	}

	for (int k = 0; k < N; k++) {
		for (int n = 0; n <N; n++) 
			DFT2[k] += complex(mixed_data[n], 0) * complex(cos((-2 * PHI*k*n) / (double)N), sin((-2 * PHI*k*n) / (double)N));
	}

	ofstream clean;
	clean.open("clean.txt");
	for (int k = 0; k < N; k++)
		clean << (double)8000 / 16000 * k << "\t" << DFT[k].mag() << DFT[k].phase() << endl;


	ofstream mixed;
	clean.open("mixed.txt");
	for (int k = 0; k < N; k++)
		clean << (double)8000 / 16000 * k << "\t" << DFT2[k].mag() << DFT2[k].phase() << endl;


	for (int k = 0; k < 16000; k++) {
		if (k > 2000 && k < 14000)
			DFT2[k] = complex(0, 0);
	}

	complex *td_music = new complex[N];



	for (int n = 0; n < 16000; n++) {
		for (int k = 0; k < 16000; k++) {
			td_music[n] += DFT2[k] * complex(cos((2 * PHI*k*n) / (double)N), sin((2 * PHI*k*n) / (double)N));
		}
		td_music[n] = td_music[n] / 16000;
	}

	short* music_data = new short[16000];
	for (int t = 0; t < 16000; t++)
		music_data[t] = td_music[t].re;

	ofstream musicout;
	musicout.open("new_music.wav", ios::binary);
	musicout.write((char*)new_header, 44);
	musicout.write((char*)music_data, sizeof(short) * 16000);
	system("pause");
	return 0;



	
	/*ifstream in;
	in.open("MixA.wav", ios::binary);
	char * header = new char[44];
	short*data = new short[16000];
	in.read((char*)header, 44);
	in.read((char*)data, sizeof(short) * 16000);

	int N = 16000;
	complex *X = new complex[16000];
	complex *x = new complex[16000];

	for (int k = 0; k < N; k++) {
		for (int n = 0; n < N; n++) {
			X[k] += complex(data[n],0) * complex(cos((-2 * PHI*k*n) / (double)N), sin((-2 * PHI*k*n) / (double)N));
		}
	}

	ofstream outfile("result.txt");
	for (int k = 0; k < 16000; k++) {
		outfile << (double)8000 / 16000 * k << "\t" << X[k].mag() << endl;
	}

	for (int k = 0; k < 16000; k++) {
		if (k > 2000 && k < 14000)
			X[k] = complex(0, 0);
	}
	
	for (int n = 0; n < 16000; n++) {
		for (int k = 0; k < 16000; k++) {
			x[n]+=X[k]*complex(cos((2 * PHI*k*n) / (double)N), sin((2 * PHI*k*n) / (double)N));
		}
		x[n] = x[n] / 16000;
	}

	short* data_ = new short[16000];
	for (int i = 0; i < 16000; i++)
		data_[i] = x[i].re;
	ofstream outfiles;
	outfiles.open("remusic.wav", ios::binary);
	outfiles.write((char*)header, 44);
	outfiles.write((char*)data_, sizeof(short) * 16000);
	system("pause");
	return 0;*/

	
}