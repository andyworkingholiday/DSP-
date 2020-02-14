#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#define PHI 3.141592

int main() {

	ofstream fout;
	fout.open("dsp_ass2.txt", ios::out);
	float t = 0, fs = 8000, dt = 1.0 / fs;
	int f0 = 400;
	//fundamental frequency는 2800과 3600의 최대공약수인 400//

	for (int i = 0; i < 400; i++, t += dt)
		fout << t << "\t" << 0.5*(sin(3600.0*PHI*t) + sin(2800.0*PHI*t)) << endl;
	
	//삼각함수의 합차 공식을 사용하였습니다//

	fout.close();
	return 0;



}