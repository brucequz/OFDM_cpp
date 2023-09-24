#include "viterbiCodec.h"

#include "feedForwardTrellis.h"

namespace {

std::vector<int> convertIntToBits(int integer, const int& length) {
  if (integer < 0) {
    std::cerr << "CANNOT CONVERT: negative integer" << std::endl;
  } else if (std::ceil(std::log2(integer + 1)) > length) {
    std::cerr << "CANNOT CONVERT: integer too large" << std::endl;
  }
  std::vector<int> result(length, 0);
  int i = length - 1;
  while (integer > 0 && i >= 0) {
    int remainder = integer % 2;
    result[i] = remainder;
    integer /= 2;
    i--;
  }
  return result;
}

int hammingDistance(const std::vector<int> x, const std::vector<int>& y) {
  assert(x.size() == y.size());
  int distance = 0;
  for (int i = 0; i < x.size(); ++i) {
    distance += (x[i] != y[i]);
  }
  return distance;
}

}  // namespace

ViterbiCodec::ViterbiCodec(int k, int n, int m, std::vector<int> poly)
      : k_(k),
      n_(n),
      m_(m) {
  trellis_ptr_ = new FeedForwardTrellis(k, n, m, poly);
  numStates_ = std::pow(2, m);
}

ViterbiCodec::~ViterbiCodec() {
  delete trellis_ptr_;
  trellis_ptr_ = nullptr;
}

std::vector<int> ViterbiCodec::encode(const std::vector<int>& message) {
  return trellis_ptr_->encode(message);
}

std::vector<int> ViterbiCodec::viterbiDecode(const std::vector<int>& coded) {
  constructTrellis(coded);
  std::vector<int> result;

  return result;
}

void ViterbiCodec::constructTrellis(const std::vector<int>& coded) {
  /*
   *  Construct a reverse trellis states with cells containing relevant
   * information for viterbi decoding. 
   *  Cell 
   * {
   *    - branchMetric:
   *    - pathMetric:
   *    - init: whether or not a cell is activated, aka. reachable in the
   * decoding process.
   * }
   */
  int message_length = coded.size();
  int total_length = message_length / n_;
  trellis_states_.resize(numStates_,
                         std::vector<Cell>(total_length + 1));
  trellis_states_[0][0].pathMetric = 0;
  trellis_states_[0][0].init = true;

  for (int cur_stage = 0; cur_stage < total_length; ++cur_stage) {
    for (int cur_state = 0; cur_state < numStates_; ++cur_state) {
      if (!trellis_states_[cur_state][cur_stage].init) {
        continue;
      }
      int cur_path_metric = trellis_states_[cur_state][cur_stage].pathMetric;
      auto begin = coded.begin() + cur_stage * n_;
      auto end = begin + n_;
      std::vector<int> target_message(begin, end);
      // activate the next states
      for (int i = 0; i < trellis_ptr_->nextStates_[cur_state].size(); ++i) {
        int next_state = trellis_ptr_->nextStates_[cur_state][i];
        trellis_states_[next_state][cur_stage + 1].init = true;

        int possible_output = trellis_ptr_->output_[cur_state][i];
        std::vector<int> expected_output = convertIntToBits(possible_output, n_);

        int branch_metric = hammingDistance(expected_output, target_message);
        int temp_path_metric = cur_path_metric + branch_metric;

        Cell* target_cell = &trellis_states_[next_state][cur_stage + 1];

        if (!target_cell->init) {
          // if the next state is not initialized, we temporarily store the path
          // metric
          target_cell->init = true;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else if (target_cell->pathMetric > temp_path_metric) {
          // the current path metric is better
          target_cell->subPathMetric = target_cell->pathMetric;
          target_cell->subFatherState = target_cell->fatherState;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else {
          // the current path metric is worse
          target_cell->subPathMetric = temp_path_metric;
          target_cell->subFatherState = cur_state;
        }
      }
    }
  }
}