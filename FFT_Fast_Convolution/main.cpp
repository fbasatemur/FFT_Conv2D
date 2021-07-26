#include <iostream>
#include <iomanip>
#include "FFT_Conv2D.h"

using namespace std;

int main() {

	int m1x = 5, m1y = 5, m2x = 3, m2y = 3;
	float* m1 = new float[m1x * m1y]{ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24 };
	float* m2 = new float[m2x * m2y]{ 0,1,2,3,4,5,6,7,8 };

	int retW, retH;
	float* ret = FFTConv2D(m1, m1x, m1y, m2, m2x, m2y, retW, retH, "valid");

	for (int i = 0; i < retH; i++)
	{
		for (int j = 0; j < retW; j++)
			cout << std::fixed << std::setprecision(2) << ret[i * retW + j] << " ";

		cout << endl;
	}

	delete[] ret, delete[] m1, delete[] m2;
	return 0;
}