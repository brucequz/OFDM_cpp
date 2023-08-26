#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "constellation.h"
#include "ofdm.h"

void print1DVector(const std::vector<int>& vec);

int main() {
  // complex number
  const std::complex<double> i(0.0, 1.0);
  std::cout << i.imag() << std::endl;
  std::cout << 5 * i.imag() << std::endl;

  // OFDM
  std::map<std::string, int> config;
  config["B"] = 10;         // number of OFDM symbol per frame
  config["L"] = 16;         // number of subcarriers per ofdm symbol
  config["CP_length"] = 4;  // cyclic prefix length
  config["Nh"] = 4;         // channel order

  Ofdm ofdm(config);

  // ------------------------------- Execution Cycle -----------------------------------
  // Symbol Error Rate for different modulation
  std::vector<double> Pe(SNR.size(), 0);
  int mc_N = 5000;  // maximum number of iterations to achieve sufficient errors
  // SNR
  double SNR_dB_start = 2.0;
  double SNR_dB_end = 20.0;
  double SNR_dB_step = 0.5;
  std::vector<double> SNR_dB = {};
  std::vector<double> SNR = {};
  for (double i = SNR_dB_start; i < SNR_dB_end; i += SNR_dB_step) {
    SNR_dB.push_back(i);
    SNR.push_back(pow(10.0, i / 10.0));
  }
  // Modulation Schemes
  std::vector<std::string> modulation_schemes = {"bpsk", "qpsk", "qam16"};

  for (std::string& mod : modulation_schemes) {
    for (int i = 0; i < SNR.size(); i++) {
      double rho = SNR[i];
      int err_cnt = 0;
      for (int mc_loop = 0; mc_loop < mc_N; mc_loop++) {
        std::vector<int> data_int = ofdm.generateRandomInt(mod);
        std::vector<complex<double>> symbol_flattened = ofdm.generateModulatedSignal(data_int, mod);
      }
    }
  }

  return 0;
}

void print1DVector(const std::vector<int>& vec) {
  for (const int& num : vec) {
    std::cout << num << " ";
  }
  std::cout << std::endl;
}