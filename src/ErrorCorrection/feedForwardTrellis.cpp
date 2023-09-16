#include "feedForwardTrellis.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../helper.h"

FeedForwardTrellis::FeedForwardTrellis(int k, int n, int m,
                                       std::vector<int> poly)
    : k_(k),
      n_(n),
      numInputSymbols_(std::pow(2, k)),
      numOutputSymbols_(std::pow(2, n)),
      numStates_(std::pow(2, m)) {
  polynomials_ = poly;
  if (polynomials_.size() != n) {
    std::cerr << "INVALID POLYNOMIAL: mismatch between number of output symbols and polynomials" << std::endl;
  }
  int min_poly = 0;
  int max_poly = decToOct(static_cast<int>(std::pow(2.0, m+1)));
  for (int poly_oct : poly) {
    if (poly_oct <= min_poly || poly_oct >= max_poly) {
      std::cerr << "INVALID POLYNOMIAL: too large or small" << std::endl;
    }
  }
  nextStates_.resize(numStates_, std::vector<int>(numInputSymbols_));
  output_.resize(numStates_, std::vector<int>(numInputSymbols_));

  computeNextStates();
  computeOutput();
  // std::cout << "printing next states" << std::endl;
  // helper.print(nextStates_);

  // std::cout << "printing output" << std::endl;
  // helper.print(output_);
}

std::vector<int> FeedForwardTrellis::encode(const std::vector<int>& message) {
  // encode assuming zero states
  std::vector<int> encoded(message.size() * (n_ / k_), 0);
  int cur_state = 0;
  for (int i = 0; i < message.size(); ++i) {
    if (message[i] >= 2) {
      std::cerr << "INVALID INPUT MESSAGE: bit input greater than 1" << std::endl;
    }
    int output = output_[cur_state][message[i]];
    int j = n_;
    while (j != 0 && output > 0) {
      int remainder = output % 2;
      encoded[i * n_ + j - 1] = remainder;
      output /= 2;
      j--;
    }
    if (output != 0) {
      std::cerr << "INVALID OUTPUT MATRIX: TOO LARGE" << std::endl;
    }
    cur_state = nextStates_[cur_state][message[i]];
  }
  return encoded;
}

void FeedForwardTrellis::computeNextStates() {
  // compute the trellis based on polynomial
  // save the result in nextStates_
  // convert polynomial to decimal

  for (int input = 0; input < numInputSymbols_; ++input) {
    for (int state = 0; state < numStates_; ++state) {
      int power = std::log2(numStates_) - 1;
      nextStates_[state][input] =
          static_cast<int>(state / 2 + input * std::pow(2, power));
    }
  }
}

void FeedForwardTrellis::computeOutput() {
  std::vector<int> poly_in_dec;
  for (int octal : polynomials_) {
    poly_in_dec.push_back(octToDec(octal));
  }
  for (int input = 0; input < numInputSymbols_; ++input) {
    for (int state = 0; state < numStates_; ++state) {
      int cur_state = state + input * static_cast<int>(
                                          std::pow(2.0, std::log2(numStates_)));
      for (int poly_id = 0; poly_id < poly_in_dec.size(); ++poly_id) {
        int result = poly_in_dec[poly_id] & cur_state;
        int count = 0;
        while (result > 0) {
          if (result & 1) {
            count++;
          }
          result >>= 1;
        }
        output_[state][input] += static_cast<int>(
            (count % 2) * std::pow(2.0, poly_in_dec.size() - poly_id - 1));
      }
    }
  }
}

int FeedForwardTrellis::octToDec(int octal) {
  int base = 1;
  int decNum = 0;

  while (octal > 0) {
    int digit = octal % 10;
    decNum += digit * base;
    base *= 8;
    octal /= 10;
  }

  return decNum;
}

int FeedForwardTrellis::decToOct(int decimal) {
  int octal = 0;
  int base = 1; // Initialize the base for the rightmost digit

  while (decimal > 0) {
      int remainder = decimal % 8; // Get the remainder when dividing by 8
      octal += remainder * base;    // Add the remainder to the octal result
      decimal /= 8;                // Divide the decimal number by 8 to move to the next digit
      base *= 10;                  // Update the base for the next digit position
  }

  return octal;
}