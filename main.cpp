#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "constellation.h"
#include "ofdm.h"

void print1DVector(const std::vector<int>& vec);
void output1DVector(std::ofstream& outputFile, const std::vector<int>& vec);
void printComplexVector2D(
    const std::vector<std::vector<std::complex<double>>>& vec2D);
void output2DComplexVector(
    std::ofstream& outputFile,
    const std::vector<std::vector<std::complex<double>>>& vec2D);

int main() {
  // output path
  std::ofstream outputFile("output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return 1;
  }
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

  // ------------------------------- Execution Cycle
  // ----------------------------------- Symbol Error Rate for different
  // modulation
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
  // Error Vector
  std::vector<double> Pe(SNR.size(), 0);
  // Modulation Schemes
  std::vector<std::string> modulation_schemes = {"bpsk", "qpsk", "qam16"};

  for (const std::string& mod : modulation_schemes) {
    for (int i = 0; i < SNR.size(); i++) {
      double rho = SNR[i];
      outputFile << "For SNR = " << rho << ": " << std::endl;
      int biterr_cnt = 0;
      for (int mc_loop = 0; mc_loop < mc_N; mc_loop++) {
        std::vector<int> data_int = ofdm.generateRandomInt(mod);
        outputFile << "Outputing integer vector" << std::endl;
        output1DVector(outputFile, data_int);
        std::vector<int> data_bits = ofdm.convertIntToBits(data_int, mod);
        outputFile << "Outputing bits vector" << std::endl;
        output1DVector(outputFile, data_bits);
        std::vector<std::vector<std::complex<double>>> data_sym =
            ofdm.generateModulatedSignal(data_int, mod);
        std::vector<std::complex<double>> data_flattened = ofdm.flattenVector(data_sym);
        // IFFT
        std::vector<std::vector<std::complex<double>>> data_ifft =
            ofdm.ifft(data_sym);
        ;
        outputFile << "Outputing ifft vector" << std::endl;
        output2DComplexVector(outputFile, data_ifft);
        break;
      }
      break;
    }
  }
  std::vector<std::vector<std::complex<double>>> input;
  std::vector<std::complex<double>> row_1 = {};
  std::vector<int> data_tx = {1,  -1, 1,  1,  -1, -1, 1, 1,
                              -1, -1, -1, -1, -1, -1, 1, -1};
  for (int i : data_tx) {
    row_1.push_back(static_cast<double>(i));
  }
  input.push_back(row_1);
  std::vector<std::vector<std::complex<double>>> result = ofdm.ifft(input);
  printComplexVector2D(result);
  outputFile.close();
  return 0;
}

void print1DVector(const std::vector<int>& vec) {
  for (const int& num : vec) {
    std::cout << num << " ";
  }
  std::cout << std::endl;
}

void output1DVector(std::ofstream& outputFile, const std::vector<int>& vec) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }

  for (int value : vec) {
    outputFile << value << " ";
  }
  outputFile << std::endl;
}

void printComplexVector2D(
    const std::vector<std::vector<std::complex<double>>>& vec2D) {
  for (const auto& row : vec2D) {
    for (const std::complex<double>& value : row) {
      std::cout << value << " ";
    }
    std::cout << std::endl;
  }
}

void output2DComplexVector(
    std::ofstream& outputFile,
    const std::vector<std::vector<std::complex<double>>>& vec2D) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }

  for (const auto& row : vec2D) {
    for (const std::complex<double>& value : row) {
      outputFile << "(" << value.real() << ", " << value.imag() << ") ";
    }
    outputFile << std::endl;
  }
}