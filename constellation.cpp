#include "constellation.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

Constellation::Constellation() {};
 
std::vector<std::complex<double>> Constellation::bpsk() {
  std::vector<std::complex<double>> constl;
  constl.push_back(std::complex<double>(1.0, 0.0));   // symbol 0
  constl.push_back(std::complex<double>(-1.0, 0.0));  // symbol 1
  return constl;
}

std::vector<std::complex<double>> Constellation::qpsk() {
  std::vector<std::complex<double>> constl;
  constl.push_back(std::complex<double>(1.0, 1.0));    // symbol 00
  constl.push_back(std::complex<double>(1.0, -1.0));   // symbol 01
  constl.push_back(std::complex<double>(-1.0, 1.0));   // symbol 10
  constl.push_back(std::complex<double>(-1.0, -1.0));  // symbol 11
  return constl;
}

std::vector<std::complex<double>> Constellation::qam(int levels) {
  if (levels < 1 || levels % 2 != 0) {
    std::cerr << "Invalid QAM level" << std::endl;
    return {};
  }

  std::vector<std::complex<double>> constl;
  int sqrtLevels = static_cast<int>(std::sqrt(levels));
  double stepSize = 2.0 / (sqrtLevels - 1);

  for (int i = 0; i < sqrtLevels; ++i) {
    for (int j = 0; j < sqrtLevels; ++j) {
      constl.push_back(
          std::complex<double>(-1.0 + i * stepSize, -1.0 + j * stepSize));
    }
  }

  return constl;
}