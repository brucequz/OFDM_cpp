#ifndef VITERBICODEC_H
#define VITERBICODEC_H

#include <vector>
#include <climits>

class FeedForwardTrellis;

struct Cell{
  int pathMetric = INT_MAX;
  int fatherState = -1;
  bool init = false;
  int subPathMetric = INT_MAX;
  int subFatherState = -1;
};

struct messageInformation {
  messageInformation() {};
  
  std::vector<int> message;
  std::vector<int> path;
  std::pair<int, int> begin_end_states;
};

class ViterbiCodec {
  public:
    ViterbiCodec(int k, int n, int m, std::vector<int> poly);
    ~ViterbiCodec();
    std::vector<int> encode(const std::vector<int>& message);
    messageInformation viterbiDecode(const std::vector<int>& coded);

    
    


  private:
    int k_; // input message length
    int n_; // output message length
    int m_; // number of memory elements
    double code_rate_;
    int numStates_;
    int list_size_;
    FeedForwardTrellis* trellis_ptr_;
    

    std::vector<std::vector<Cell>> constructTrellis(const std::vector<int>& coded);
};


#endif