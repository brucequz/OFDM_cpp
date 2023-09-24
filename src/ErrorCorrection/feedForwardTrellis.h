
#ifndef FEEDFORWARDTRELLIS_H
#define FEEDFORWARDTRELLIS_H

#include <vector>
#include "../helper.h"
#include "viterbiDecoder.h"
#include "viterbiCodec.h"


struct FeedForwardTrellis{
  FeedForwardTrellis(int k, int n, int m, std::vector<int> poly);
  // ~FeedForwardTrellis(); 

  std::vector<int> encode(const std::vector<int>& message);

  private:
    int k_;
    int n_;
    int numInputSymbols_;
    int numOutputSymbols_;
    int numStates_;
    int code_rate_;
    std::vector<int> polynomials_;  // generator polynomial in octal
    std::vector<std::vector<int>> nextStates_;
    std::vector<std::vector<int>> output_;

    void computeNextStates();
    void computeOutput();
    Helper helper;
    friend class ViterbiDecoder;
    friend class ViterbiCodec;
};

#endif