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

## Debug Info
I use the built-in lldb debugger for this project. 

To start debugging, use 'lldb ofdm' after cmake. Use "b ${file_name} : ${line_number}" to set a breakpoint. After setting breakpoints, use "r" to run the debugger.

* s: step in
* c: continue
* n: step over 
* p ${var_name}: print variable
* q: quit debugger
## Note:

1. MATLAB does column-major reshape operations, while C++ performs row-major operation. Be cautious!
2. filter() in ofdm.cpp does convolution by brute force. Truncate the last $(CP_length_) elements to get the desired behavior.
3. generateAWGN() now uses std::mt19937() random number generator. This is the key to ensure similar behavior with respect to MATLAB implementation. Using default rng would, for example, results in double symbol errors when running bpsk for 5000 iterations.
4. Why would not multiplying channel response by sqrt(rho) significantly impact 16qam but not bpsk or qpsk?

(10 x 16)
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  <---- an OFDM symbol
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
                                  .
                                  .
                                  .
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  * 
             *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *      

                      ^
                      |
                      |
                      |
                      one subcarrier