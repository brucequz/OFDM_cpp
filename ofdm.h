#ifndef OFDM
#define OFDM

#include <fftw3.h>

#include <map>
#include <string>
#include <unordered_map>

#include "constellation.h"

class Ofdm {
 public:
  // Constructor
  Ofdm(const std::map<std::string, int>& config);

  // Destructor
  ~Ofdm();

  // random number generator
  std::vector<int> generateRandomInt(const std::string& constl_type);

  std::vector<std::vector<std::complex<double>>> generateModulatedSignal(
      const std::vector<int>& integers, const std::string& constl_type);

  // Helper functions
  std::vector<std::complex<double>> flattenVector(
      const std::vector<std::vector<std::complex<double>>>& input);

  std::vector<int> convertIntToBits(const std::vector<int>& integers,
                                    const std::string& constl_type);

  std::vector<std::vector<std::complex<double>>> transpose2DComplexVector(
      const std::vector<std::vector<std::complex<double>>>& input);

  std::vector<std::vector<std::complex<double>>> reshape(
      const std::vector<std::complex<double>>& input, int num_rows,
      int num_cols);

  std::vector<std::vector<std::complex<double>>> columnMajorReshape(
      const std::vector<std::complex<double>>& input, int num_rows,
      int num_cols);

  // FFT and IFFT
  std::vector<std::vector<std::complex<double>>> fft(
      const std::vector<std::vector<std::complex<double>>>& input);

  std::vector<std::vector<std::complex<double>>> ifft(
      const std::vector<std::vector<std::complex<double>>>& input);

  // Cyclic Prefix
  std::vector<std::vector<std::complex<double>>> addCyclicPrefix(
      std::vector<std::vector<std::complex<double>>> vec2D);

  std::vector<std::vector<std::complex<double>>> removeCyclicPrefix(
      std::vector<std::complex<double>> input);

  // AWGN Noise
  std::vector<std::complex<double>> addAWGN(
      const std::vector<std::complex<double>>& signal);

  // Filter
  std::vector<std::complex<double>> filter(
      const std::vector<std::complex<double>>& signal,
      const std::vector<std::complex<double>>& filter_coeffs);

 private:
  Constellation* constl_;
  std::unordered_map<std::string, std::vector<std::complex<double>>>
      constellations_;  // modulation schemes
  int B_;               // number of OFDM symbols per transmitted frame
  int L_;               // number of subcarriers in each OFDM symbol
  int CP_length_;       // cyclic prefix length
  int Nh_;              // channel order
  std::vector<int> convertBits(int value, int num_bits);

  // You can also add private member functions here
};

#endif
