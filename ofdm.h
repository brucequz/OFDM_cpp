#ifndef OFDM
#define OFDM

#include <map>
#include <string>
#include <unordered_map>

#include "constellation.h"
#include <fftw3.h>

class Ofdm {
 public:
  // Constructor
  Ofdm(const std::map<std::string, int>& config);
  // Destructor
  ~Ofdm();
  // random number generator
  std::vector<int> generateRandomInt(const std::string& constl_type);
  std::vector<std::complex<double>> generateModulatedSignal(
      const std::vector<int>& integers, const std::string& constl_type);
  // Helper functions
  std::vector<std::vector<std::complex<double>>> reshapeVector(
      const std::vector<std::complex<double>>& input, int rows, int cols);
  std::vector<int> convertIntToBits(const std::vector<int>& integers,
                                    const std::string& constl_type);
  // FFT and IFFT

  std::vector<std::vector<std::complex<double>>> ifft(
      const std::vector<std::vector<std::complex<double>>>& input);

 private:
  Constellation* constl_;
  std::unordered_map<std::string, std::vector<std::complex<double>>>
      constellations_;  // modulation schemes
  int B_;               // number of OFDM symbols per transmitted frame
  int L_;               // number of subcarriers in each OFDM symbol
  int CP_length_;       // cyclic prefix length
  int Nh_;              // channel order
                        // You can add more private member variables here
  std::vector<int> convertBits(int value, int num_bits);

  // You can also add private member functions here
};

#endif
