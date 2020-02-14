#include <iostream>
#include <fstream>
using namespace std;
#define W 400
#define H 300
#define WH 120000

int main() {
	/*int i, j, m, individual_idx, compositive_idx;
	ofstream fff("video.rgb", ios::out | ios::binary);
	unsigned char* R = new unsigned char[WH];
	unsigned char* G = new unsigned char[WH];
	unsigned char* B = new unsigned char[WH];
	unsigned char* RGB = new unsigned char[WH * 3];




	for (m = 5; m < 370; m += 5) {
		//축구장의 배경색입니다//
		for (i = 0; i < WH; i++) {
			R[i] = 41;
			G[i] = 230;
			B[i] = 15;
		}

		for (i = 0; i < WH; i += 400) {
			// 잔디의 색깔의 다름을 표현하기 위해 비슷한 색을 섞었습니다.
			for (j = 0; j < 80; j++) {
				R[i + j] = 28;
				G[i + j] = 156;
				B[i + j] = 10;
			}
		}

		for (i = 160; i < WH; i += 400) {
			//3부분정도 달라집니다
			for (j = 0; j < 80; j++) {
				R[i + j] = 28;
				G[i + j] = 156;
				B[i + j] = 10;
			}
		}

		for (i = 320; i < WH; i += 400) {
			for (j = 0; j < 80; j++) {
				R[i + j] = 28;
				G[i + j] = 156;
				B[i + j] = 10;
			}
		}

		for (i = 120; i < 180; i++) {
			//페널티 에이리어 입니다
			for (j = 0; j < 30; j++) {
				R[i * 400 + j] = 255;
				G[i * 400 + j] = 255;
				B[i * 400 + j] = 255;
			}
		}

		for (i = 120; i < 180; i++) {
			for (j = 399; j > 369; j--) {
				R[i * 400 + j] = 255;
				G[i * 400 + j] = 255;
				B[i * 400 + j] = 255;
			}
		}

		for (i = 123; i < 177; i++) {
			for (j = 0; j < 27; j++) {
				R[i * 400 + j] = 28;
				G[i * 400 + j] = 156;
				B[i * 400 + j] = 10;
			}
		}

		for (i = 123; i < 177; i++) {
			for (j = 399; j > 372; j--) {
				R[i * 400 + j] = 28;
				G[i * 400 + j] = 156;
				B[i * 400 + j] = 10;
			}
		}


		for (i = 0; i < WH; i++) {
			//중앙에 센터서클
			int q = i / 400;
			int r = i % 400;
			int garo = abs(199 - r);
			int sero = abs(149 - q);
			if (sqrt(garo*garo + sero * sero) < 30 && sqrt(garo*garo + sero * sero) > 27) {
				R[i] = 255;
				G[i] = 255;
				B[i] = 255;
			}
		}

		for (i = 198; i < WH; i += 400) {
			//중앙선 입니닷
			for (j = 0; j < 4; j++) {
				R[i + j] = 255;
				G[i + j] = 255;
				B[i + j] = 255;
			}

		}
		
		for (i = 0; i < WH; i++) {
			//움직이는 축구공을 만들겁니다. 노란색 공입니다.
			int q = i / 400;
			int r = i % 400;
			int garo = abs(25 - r);
			int sero = abs(149 - q);
			if (sqrt(garo*garo + sero * sero) < 9) {
				individual_idx = i + m;
				R[individual_idx] = 255;
				G[individual_idx] = 255;
				B[individual_idx] = 0;
			}
		}

		for (i = compositive_idx = 0; i < WH; i++, compositive_idx += 3) {
			RGB[compositive_idx] = R[i];
			RGB[compositive_idx + 1] = G[i];
			RGB[compositive_idx + 2] = B[i];
		}

		fff.write((const char*)RGB, WH * 3); //파일을 저장할 메모리 크기 wh
	}

	fff.close();
	return 0;*/


	int i, j, m, individual_idx, compositive_idx;
	ofstream fff;
	fff.open("videoSS.rgb", ios::out | ios::binary);
	unsigned char R[W*H], G[W*H], B[W*H], RGB[3 * WH];

	for (m = 0; m < 200; m += 5) {
		for (int i = 0; i < WH; i++) {
			R[i] = 255;
			G[i] = 255;
			B[i] = 0;
		}

		for (i = 80; i < 100; i++) {
			for (j = 100; j < 120; j++) {
				individual_idx = i * WIDTH + j + m;
				R[individual_idx] = 255;
				G[individual_idx] = 0;
				B[individual_idx] = 0;

			}
		}

		for (i = 0; i < 120000; i++, compositive_idx = 0; compositive_idx += 3) {
			RGB[compositive_idx] = R[i];
			RGB[compositive_idx + 1] = G[i];
			RGB[compositive_idx + 2] = B[i];
		}

		fff.write((const char*)RGB, 360000);

	}
}
