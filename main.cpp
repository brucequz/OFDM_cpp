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
void output1DComplexVector(std::ofstream& outputFile, const std::vector<std::complex<double>>& vec);
void output2DComplexVector(
    std::ofstream& outputFile,
    const std::vector<std::vector<std::complex<double>>>& vec2D);
std::vector<std::complex<double>> normalize(
    const std::vector<std::complex<double>>& input);

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
  double SNR_dB_start = 0.0;
  double SNR_dB_end = 20.0;
  double SNR_dB_step = 2.0;
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
        std::vector<std::complex<double>> data_flattened =
            ofdm.flattenVector(data_sym);

        // IFFT
        std::vector<std::vector<std::complex<double>>> data_ifft =
            ofdm.ifft(data_sym);
        outputFile << "Outputing ifft vector" << std::endl;
        output2DComplexVector(outputFile, data_ifft);

        // cyclic prefix
        std::vector<std::vector<std::complex<double>>> data_cp =
            ofdm.addCyclicPrefix(data_ifft);
        outputFile << "Outputing cp vector" << std::endl;
        output2DComplexVector(outputFile, data_cp);

        // Reshape the BxN matrix to obtain the frame (1xTotal_length)
        // Total_length = (CP_length+L)*B
        std::vector<std::complex<double>> data_tx = ofdm.flattenVector(ofdm.transpose2DComplexVector(data_cp));

        // Define the channel response
        // deterministic frequency-selective channel
        std::complex<double> h1(3.0, 0.0);
        std::complex<double> h2(-1.0 *
                                std::exp(std::complex<double>(0.0, 0.13)));
        std::complex<double> h3(1.0 *
                                std::exp(std::complex<double>(0.0, -0.35)));
        std::complex<double> h4(0.0, 0.0);
        std::complex<double> h5(4.0 *
                                std::exp(std::complex<double>(0.0, 1.03)));

        // Create a vector of complex numbers
        std::vector<std::complex<double>> h = {h1, h2, h3, h4, h5};
        std::vector<std::complex<double>> h_normalized = normalize(h);

        std::vector<std::complex<double>> rec =
            ofdm.filter(data_tx, h_normalized);
        for (std::complex<double>& value : rec) {
          value *= std::sqrt(rho);
        }

        // Add Noise
        std::vector<std::complex<double>> received = ofdm.addAWGN(rec);
        outputFile << "Outputing received symbols (with noise)" << std::endl;
        output1DComplexVector(outputFile, received);

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

void output1DComplexVector(std::ofstream& outputFile, const std::vector<std::complex<double>>& vec) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }
  for (const auto& value : vec) {
    outputFile << "(" << value.real() << ", " << value.imag() << ") ";
  }
  outputFile << std::endl;
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

std::vector<std::complex<double>> normalize(
    const std::vector<std::complex<double>>& input) {
  double sum_of_abs_squares = 0.0;
  for (const std::complex<double>& value : input) {
    sum_of_abs_squares += std::norm(value);
  }

  double normalization_factor = std::sqrt(sum_of_abs_squares);

  std::vector<std::complex<double>> normalized_vector;
  for (const std::complex<double>& value : input) {
    normalized_vector.push_back(value / normalization_factor);
  }

  return normalized_vector;
}