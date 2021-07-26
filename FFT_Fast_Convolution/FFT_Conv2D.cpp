#include "FFT_Conv2D.h"
#include <iostream>
#include <math.h>
#include "Process.h"


struct Output {
	float* pInput;
	float* pKernel;
	int outWidth;
	int outHeight;
};

void ComplexMult(const float* fft1Real, const float* fft1Imag, const float* fft2Real, const float* fft2Imag, float* outReal, float* outImag, int width, int height) {

	long imageSize = width * height;
	for (long i = 0; i < imageSize; i++)
	{
		outReal[i] = fft1Real[i] * fft2Real[i] - fft1Imag[i] * fft2Imag[i];
		outImag[i] = fft1Real[i] * fft2Imag[i] + fft1Imag[i] * fft2Real[i];
	}
}

void CircleShift(float* output, const float* input, int rows, int cols, int yshift, int xshift)
{
	for (int r = 0; r < rows; r++) {

		int newR = (r + yshift) % rows;
		if (newR < 0) newR = rows + newR;

		for (int c = 0; c < cols; c++) {

			int newC = (c + xshift) % cols;
			if (newC < 0) newC = cols + newC;

			output[newR * cols + newC] = input[r * cols + c];
		}
	}
}

float* ZeroExtend(const float* x, int sizeX, int sizeY, int exWidth, int exHeight) {

	float* extX = new float[exHeight * exWidth]();
	unsigned int byteSizeX = sizeX * sizeof(*x);

	for (int r = 0; r < sizeY; r++)
		memcpy((extX + r * exWidth), (x + r * sizeX), byteSizeX);

	return extX;
}

float* KernelShift(const float* x, int sizeX, int sizeY, int shiftX, int shiftY) {

	float* circExt = new float[sizeX * sizeY];
	CircleShift(circExt, x, sizeY, sizeX, -shiftY, -shiftX);
	return circExt;
}

Output* ProcessInputs(const float* input, int inputWidth, int inputHeight, const float* kernel, int kernelWidth, int kernelHeight, const char* shape) {

	Output* ret = new Output;

	if (!strcmp(shape, "full"))
	{
		int exWidth = inputWidth + kernelWidth - 1;
		int exHeight = inputHeight + kernelHeight - 1;
		ret->pInput = ZeroExtend(input, inputWidth, inputHeight, exWidth, exHeight);
		ret->pKernel = ZeroExtend(kernel, kernelWidth, kernelHeight, exWidth, exHeight);

		ret->outHeight = exHeight;
		ret->outWidth = exWidth;
	}
	else if (!strcmp(shape, "same"))
	{
		int exWidth = inputWidth, exHeight = inputHeight;
		int mmidX = floor(float(kernelWidth) / 2.0F);
		int mmidY = floor(float(kernelHeight) / 2.0F);
		exWidth += mmidX;
		exHeight += mmidY;
		ret->pInput = ZeroExtend(input, inputWidth, inputHeight, exWidth, exHeight);
		float* extKernel = ZeroExtend(kernel, kernelWidth, kernelHeight, exWidth, exHeight);
		ret->pKernel = KernelShift(extKernel, exWidth, exHeight, mmidX, mmidY);
		delete[] extKernel;

		ret->outHeight = exHeight;
		ret->outWidth = exWidth;
	}
	else if (!strcmp(shape, "valid"))
	{
		int mmidX = kernelWidth - 1;
		int mmidY = kernelHeight - 1;
		ret->pInput = new float[inputHeight * inputWidth];
		memcpy(ret->pInput, input, inputHeight * inputWidth * sizeof(float));
		float* extKernel = ZeroExtend(kernel, kernelWidth, kernelHeight, inputWidth, inputHeight);
		ret->pKernel = KernelShift(extKernel, inputWidth, inputHeight, mmidX, mmidY);
		delete[] extKernel;

		ret->outHeight = inputHeight;
		ret->outWidth = inputWidth;
	}
	else exit(-1);

	return ret;
}

/// <summary>
/// For shape options 'full', 'same' and 'valid' this function return the same result as "Time Domain Convolution".
/// FFTConv2D, uses multiplication in the frequency domain to compute the convolution.
/// It may be faster than "Time Domain Convolution" for masks above a certain size.
/// This should be checked experimentally for any given applicationand system.
/// </summary>
/// <param name="shape">
/// 'full'  - (default) returns the full 2-D convolution,
/// 'same'  - returns the central part of the convolution that is the same size as "input"(using zero padding),
/// 'valid' - returns only those parts of the convolution that are computed without the zero - padded edges.
/// </param>
/// <returns>
/// convolved result memory address returns
/// </returns>
float* FFTConv2D(const float* input, int inputWidth, int inputHeight, const float* kernel, int kernelWidth, int kernelHeight, int& outputWidth, int& outputHeight, const char* shape) {

	Output* outputs = ProcessInputs(input, inputWidth, inputHeight, kernel, kernelWidth, kernelHeight, shape);
	int outHeight = outputs->outHeight;
	int outWidth = outputs->outWidth;
	int outSize2D = outWidth * outHeight;

	float* m1_OR = new float[outSize2D];
	float* m1_OI = new float[outSize2D];
	float* m2_OR = new float[outSize2D];
	float* m2_OI = new float[outSize2D];

	FFT2D(outputs->pInput, m1_OR, m1_OI, outWidth, outHeight);
	FFT2D(outputs->pKernel, m2_OR, m2_OI, outWidth, outHeight);

	delete[] outputs->pInput;
	delete[] outputs->pKernel;
	delete[] outputs;

	float* multReal = new float[outSize2D];
	float* multImag = new float[outSize2D];

	ComplexMult(m1_OR, m1_OI, m2_OR, m2_OI, multReal, multImag, outWidth, outHeight);

	delete[] m1_OR, delete[] m1_OI;
	delete[] m2_OR, delete[] m2_OI;

	float* outReal = new float[outSize2D];
	float* outImag = new float[outSize2D];
	IFFT2D(outReal, outImag, multReal, multImag, outWidth, outHeight);

	delete[] multReal, delete[] multImag;
	delete[] outImag;

	float* ret = NULL;

	if (!strcmp(shape, "full"))
	{
		outputWidth = inputWidth + kernelWidth - 1;
		outputHeight = inputHeight + kernelHeight - 1;
		ret = outReal;
	}
	else if (!strcmp(shape, "same"))
	{
		outputWidth = inputWidth;
		outputHeight = inputHeight;
		ret = new float[outputWidth * outputHeight];
		unsigned int byteWidth = outputWidth * sizeof(*outReal);

		for (int r = 0; r < outputHeight; r++)
			memcpy((ret + r * outputWidth), (outReal + r * outWidth), byteWidth);

		delete[] outReal;
	}
	else if (!strcmp(shape, "valid"))
	{
		outputWidth = inputWidth - kernelWidth + 1;
		outputHeight = inputHeight - kernelHeight + 1;
		ret = new float[outputWidth * outputHeight];
		unsigned int byteWidth = outputWidth * sizeof(*outReal);

		for (int r = 0; r < outputHeight; r++)
			memcpy((ret + r * outputWidth), (outReal + r * outWidth), byteWidth);

		delete[] outReal;
	}
	else exit(-1);

	return ret;
}