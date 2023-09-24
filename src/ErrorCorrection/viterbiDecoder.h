#ifndef VITERBIDECODER_H
#define VITERBIDECODER_H

#include <vector>

struct Cell;
class FeedForwardTrellis;

class ViterbiDecoder {
  public:
    ViterbiDecoder(const FeedForwardTrellis& FFT);
    void constructDecodeTrellis(const std::vector<int>& received);
    std::vector<int> decode(const std::vector<int>& received);
    std::vector<std::vector<Cell>> getTrellis() const {return trellis_states_;}

  private:
    int n_;
    int numStates_;
    const FeedForwardTrellis* trellis_;
    std::vector<std::vector<Cell>> trellis_states_;

};

struct Cell{
  int pathMetric = INT_MAX;
  int fatherState = -1;
  bool init = false;
  int subPathMetric = INT_MAX;
  int subFatherState = -1;
};



#endif