#include "feedForwardTrellis.h"

#include <algorithm>
#include <cmath>
#include <iostream>

#include "../helper.h"

FeedForwardTrellis::FeedForwardTrellis(int k, int n, int m,
                                       std::vector<int> poly)
    : numInputSymbols_(std::pow(2, k)),
      numOutputSymbols_(std::pow(2, n)),
      numStates_(std::pow(2, m)) {
  polynomials_ = poly;
  nextStates_.resize(numStates_, std::vector<int>(numInputSymbols_));
  output_.resize(numStates_, std::vector<int>(numInputSymbols_));

  computeNextStates();
  computeOutput();
  // std::cout << "printing next states" << std::endl;
  // helper.print(nextStates_);

  // std::cout << "printing output" << std::endl;
  // helper.print(output_);
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
  for (int input = 0; input < (numOutputSymbols_ / numInputSymbols_); ++input) {
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