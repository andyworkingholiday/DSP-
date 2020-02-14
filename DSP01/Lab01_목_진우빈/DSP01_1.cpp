#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#define PHI 3.141592

int main() {

	ofstream fout;
	fout.open("dsp_ass1.txt", ios::out);
	float t = 0, fs = 200.0, dt = 1.0 / fs;
	//fs를 200으로 설정해 주었습니다.//
	int f0 = 5, n = 3, smp_cnt;
	//1번의 신호들의 주파수의 최대공약수인 5가 기본주파수 입니다.//
	smp_cnt = (fs / f0)*n;

	for (int i = 0; i <= smp_cnt; i++, t += dt)
		fout << t << "\t" << 2 + 4.0*cos(30.0*PHI*t - 0.2*PHI) + 3.0*sin(40.0*PHI*t) + 4.0*cos(60.0*PHI*t - (PHI / 3)) << endl;

	fout.close();
	return 0;



}