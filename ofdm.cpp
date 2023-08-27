#include "ofdm.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>

Ofdm::Ofdm(const std::map<std::string, int>& config) {
  auto it = config.find("B");
  B_ = (it != config.end()) ? it->second : 10;

  it = config.find("L");
  L_ = (it != config.end()) ? it->second : 16;

  it = config.find("CP_length");
  CP_length_ = (it != config.end()) ? it->second : 4;

  it = config.find("Nh");
  Nh_ = (it != config.end()) ? it->second : 4;

  // constellation
  constl_ = new Constellation();
  constellations_["bpsk"] = constl_->bpsk();
  constellations_["qpsk"] = constl_->qpsk();
  constellations_["qam16"] = constl_->qam(16);
}

Ofdm::~Ofdm() {
  delete constl_;
  constl_ = nullptr;
}

std::vector<int> Ofdm::generateRandomInt(const std::string& constl_type) {
  // TODO: generate bits such that all modulation compare the number of
  // iterations
  //       when they make the same number of errors

  std::vector<int> randomIntegers;

  randomIntegers.reserve(B_ * L_);
  for (int i = 0; i < B_ * L_; ++i) {
    randomIntegers.push_back(std::rand() % constellations_[constl_type].size());
  }

  return randomIntegers;
}

std::vector<std::vector<std::complex<double>>> Ofdm::generateModulatedSignal(
    const std::vector<int>& integers, const std::string& constl_type) {
  std::vector<std::vector<std::complex<double>>> modulated_signal;
  std::vector<std::complex<double>> constellation =
      constellations_[constl_type];

  int total_elements = B_ * L_;

  if (integers.size() != total_elements) {
    std::cout << "Cannot generate modulated signal. Input vector size does not "
                 "match target size.\n";
    return {};
  }

  modulated_signal.resize(B_, std::vector<std::complex<double>>(L_));

  for (int i = 0; i < B_; ++i) {
    for (int j = 0; j < L_; ++j) {
      modulated_signal[i][j] = constellation[integers[i * L_ + j]];
    }
  }

  return modulated_signal;
}

std::vector<std::complex<double>> Ofdm::flattenVector(
    const std::vector<std::vector<std::complex<double>>>& input) {
  std::vector<std::complex<double>> flattened_vector;
  int rows = input.size();
  int cols = input[0].size();

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      flattened_vector.push_back(input[i][j]);
    }
  }

  return flattened_vector;
}

std::vector<int> Ofdm::convertBits(int value, int num_bits) {
  std::vector<int> bits;
  for (int i = num_bits - 1; i >= 0; --i) {
    int bit = (value >> i) & 1;
    bits.push_back(bit);
  }
  return bits;
}

std::vector<int> Ofdm::convertIntToBits(const std::vector<int>& integers,
                                        const std::string& constl_type) {
  std::vector<int> bits;

  int num_bits_per_symbol = 0;
  if (constl_type == "bpsk") {
    num_bits_per_symbol = 1;
  } else if (constl_type == "qpsk") {
    num_bits_per_symbol = 2;
  } else if (constl_type == "qam16") {
    num_bits_per_symbol = 4;
  }

  for (int value : integers) {
    std::vector<int> value_bits = convertBits(value, num_bits_per_symbol);
    bits.insert(bits.end(), value_bits.begin(), value_bits.end());
  }

  return bits;
}

std::vector<std::vector<std::complex<double>>> Ofdm::ifft(
    const std::vector<std::vector<std::complex<double>>>& input) {
  /*
    normalization factor sqrt(L) is incorporated in ifft
    " * sqrt(static_cast<double>(L_)); "
  */
  int num_rows = input.size();
  int num_cols = input[0].size();

  fftw_plan plan = fftw_plan_dft_1d(num_cols, nullptr, nullptr, FFTW_BACKWARD,
                                    FFTW_ESTIMATE);

  std::vector<std::vector<std::complex<double>>> output(
      num_rows, std::vector<std::complex<double>>(num_cols));

  for (int i = 0; i < num_rows; ++i) {
    std::vector<std::complex<double>> row_input(num_cols);
    std::vector<std::complex<double>> row_output(num_cols);

    for (int j = 0; j < num_cols; ++j) {
      row_input[j] = input[i][j];
    }

    fftw_complex* in = reinterpret_cast<fftw_complex*>(&row_input[0]);
    fftw_complex* out = reinterpret_cast<fftw_complex*>(&row_output[0]);

    fftw_execute_dft(plan, in, out);

    for (int j = 0; j < num_cols; ++j) {
      output[i][j] = row_output[j] / static_cast<double>(num_cols) *
                     sqrt(static_cast<double>(L_));
    }
  }

  fftw_destroy_plan(plan);

  return output;
}
