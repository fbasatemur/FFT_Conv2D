# Fast 2D Convolution
[<img src="https://img.shields.io/badge/FFTW-3.3.5-76B900?style=for-the-badge" style="vertical-align:top margin:6px 4px">](http://www.fftw.org/install/windows.html)

Implementation of 2D convolution using Fast Fourier Transformation (FFT). Convolution in time space is equal to multiplication in Frequency space. 

![formula](https://render.githubusercontent.com/render/math?math=\huge\color{Red}f(x,y)*h(x,y)=F(u,v).H(u,v))

This equation is expressed below.

<img src="https://www.programmersought.com/images/174/b9c351b333914cf23ff3aa43d2049a36.png" width = "700" height = "325">

## Run

The FFTW3 library in external_lib is compiled for Windows. For Ubuntu, FFTW 3.3.5 should be installed.  
Run steps for Windows: 
```shell
mkdir build && cd build
copy ..\external_lib\fftw-3.3.5-dll64\libfftw3f-3.dll .
cmake ..
make
```

***Or run with Clion***  
**Note**: A copy of libfftw3f-3.dll must be in the executable directory

## Defined Convulation Skils
:white_check_mark: ***full***: (default) returns the full 2-D convolution  
:white_check_mark: ***same***: returns the central part of the convolution that is the same size as "input"(using zero padding)  
:white_check_mark: ***valid***: returns only those parts of the convolution that are computed without the zero - padded edges.  

***full***             |  ***same***         | ***valid***
:-------------------------:|:-------------------------:|:-------------------------:
![](https://www.programmersought.com/images/370/4dff5e0f8d27089c57d46e2417ccbe62.png)  |  ![](https://www.programmersought.com/images/574/c30864b7815fea6b8479468f6651e906.png)  |  ![](https://www.programmersought.com/images/570/e7f858b23ba84e08db41c683711145aa.png)

## Result
This process may faster than Time Domain Convolution. When using large images (exp: 1024 x 1024 or 512 x 512 etc.) is about 2 times faster than Time Domain Convolution. This measurement is based on trial and error, you should observe speed up yourself.  
Example usage of FFTConv2D() function is specified in <a href="https://github.com/fbasatemur/FFT_Conv2D/blob/main/FFT_Fast_Convolution/main.cpp"> main</a>.
