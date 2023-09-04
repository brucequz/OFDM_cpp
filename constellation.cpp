#include "constellation.h"

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

Constellation::Constellation(){};

std::vector<std::complex<double>> Constellation::bpsk() {
  std::vector<std::complex<double>> constl;
  constl.push_back(std::complex<double>(-1.0, 0.0));  // symbol 0
  constl.push_back(std::complex<double>(1.0, 0.0));   // symbol 1
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
  double stepSize = std::sqrt(8.0 / 5.0);

  for (int i = 0; i < sqrtLevels; ++i) {
    for (int j = 0; j < sqrtLevels; ++j) {
      constl.push_back(std::complex<double>(-1.5 * stepSize + i * stepSize,
                                            -1.5 * stepSize + j * stepSize));
    }
  }

  return constl;
}

double Constellation::calculateAverageBitEnergy(
    const std::vector<std::complex<double>>& constellation) {
  // Get the number of symbols in the constellation
  int numSymbols = static_cast<int>(constellation.size());

  // Calculate the number of bits per symbol (log2)
  int bitsPerSymbol = static_cast<int>(std::log2(numSymbols));

  // Check if the number of symbols is not a power of 2
  if (numSymbols <= 0 || (numSymbols & (numSymbols - 1)) != 0) {
    std::cerr << "Error: Number of symbols should be a power of 2."
              << std::endl;
    return -1.0;  // Error
  }

  // Calculate the average symbol energy
  double averageSymbolEnergy = 0.0;

  for (const std::complex<double>& point : constellation) {
    // Calculate the magnitude (amplitude) of the complex point
    double magnitude = std::abs(point);

    // Add the squared magnitude to the average symbol energy
    averageSymbolEnergy += magnitude * magnitude;
  }

  // Normalize by dividing by the number of symbols
  averageSymbolEnergy /= numSymbols;

  // Calculate the average bit energy by dividing by the number of bits per
  // symbol
  double averageBitEnergy = averageSymbolEnergy / bitsPerSymbol;

  return averageBitEnergy;
}