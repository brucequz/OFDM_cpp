# OFDM_cpp


## Dependencies

1. fftw3 : Fast Fourier Transform library
2. MATLAB Data API

## Cmake commands
mkdir build
cd build
cmake ../
make &or& make --build .
./ofdm
## Note:

1. MATLAB does column-major reshape operations, while C++ performs row-major operation. Be cautious!
2. filter() in ofdm.cpp does convolution by brute force. Truncate the last $(CP_length_) elements to get the desired behavior.


