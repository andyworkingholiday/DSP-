#include <iostream>
using namespace std;
#define PHI 3.141592
int main() {
	float t = 0.0, fs = 120; dt = 1.0 / fs, f0 = 5;
	int n = 3, smp_cnt = n * (fs / f0);

	for (int i = 0; i < smp_cnt; i++, t += dt) {
		cout << t << "\t" << cos(2 * 15 * PHI*t) + sin(2 * 20 * PHI*t) << endl;
	}

	return 0;
}