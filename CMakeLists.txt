cmake_minimum_required(VERSION 3.15)
set(PROJECT_NAME FFT_Conv2D)
project(${PROJECT_NAME})

set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)

set(HEADERS
        include/FFT_Conv2D.h
        include/fftwf_.h)

set(SOURCES
        src/FFT_Conv2D.cpp
        main.cpp)

add_executable(${PROJECT_NAME} ${HEADERS} ${SOURCES})

message("-- Linking FFTW3")
target_include_directories(${PROJECT_NAME} PRIVATE "external_lib/fftw-3.3.5-dll64/")
target_link_directories(${PROJECT_NAME} PRIVATE "external_lib/fftw-3.3.5-dll64/")
target_link_libraries(${PROJECT_NAME} PRIVATE libfftw3f-3)
