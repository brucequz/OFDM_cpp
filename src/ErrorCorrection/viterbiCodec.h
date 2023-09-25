#ifndef VITERBICODEC_H
#define VITERBICODEC_H

#include "feedForwardTrellis.h"
#include "viterbiDecoder.h"

struct Cell;
struct messageInformation;

class ViterbiCodec {
  public:
    ViterbiCodec(int k, int n, int m, std::vector<int> poly);
    ~ViterbiCodec();
    std::vector<int> encode(const std::vector<int>& message);
    std::vector<int> viterbiDecode(const std::vector<int>& coded);

    // temporary getters
    std::vector<std::vector<Cell>> getTrellis() {return trellis_states_;}
    
    


  private:
    int k_; // input message length
    int n_; // output message length
    int m_; // number of memory elements
    int numStates_;
    int list_size_;
    FeedForwardTrellis* trellis_ptr_;
    std::vector<std::vector<Cell>> trellis_states_;

    void constructTrellis(const std::vector<int>& coded);
};

struct messageInformation {
  messageInformation();
  
  std::vector<int> message;
  std::vector<int> path;
  std::pair<int, int> begin_end_states;
};

#endif