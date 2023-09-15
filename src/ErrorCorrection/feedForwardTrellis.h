
#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H

#include <vector>

struct FeedForwardTrellis{
  FeedForwardTrellis(int k, int n, int m, std::vector<int> poly);
  ~FeedForwardTrellis();


  private:
    int numInputSymbols_;
    int numOutputSymbols_;
    int numStates_;
    std::vector<int> polynomials_;  // generator polynomial in octal
    std::vector<std::vector<int>> nextStates_;
    std::vector<std::vector<int>> output_;

    void computeNextStates();
    int octToDec(int octal);

};

#endif