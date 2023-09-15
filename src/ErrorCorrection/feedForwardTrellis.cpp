#include "feedForwardTrellis.h"

#include <algorithm>
#include <cmath>
#include <iostream>

FeedForwardTrellis::FeedForwardTrellis(int k, int n, int m, std::vector<int> poly)
    : numInputSymbols_(std::pow(2, k)),
      numOutputSymbols_(std::pow(2, n)),
      numStates_(std::pow(2, m)) {
  polynomials_ = poly;
  nextStates_.resize(numStates_, std::vector<int>(numInputSymbols_));
  output_.resize(numStates_, std::vector<int>(numInputSymbols_));
}


void FeedForwardTrellis::computeNextStates() {
  // compute the trellis based on polynomial
  // save the result in nextStates_
  // convert polynomial to decimal
  std::vector<int> polyInDec;
  for (int octal : polynomials_) {
    polyInDec.push_back(octToDec(octal));
  }
  for (int input = 0; input < numInputSymbols_; ++input) {
    for (int state = 0; state < numStates_; ++state) {
      nextStates_[state][input] = polyInDec[input] ^ state;
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