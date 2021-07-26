
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
float* FFTConv2D(const float* input, int inputWidth, int inputHeight, const float* kernel, int kernelWidth, int kernelHeight, int& outputWidth, int& outputHeight, const char* shape = "full");