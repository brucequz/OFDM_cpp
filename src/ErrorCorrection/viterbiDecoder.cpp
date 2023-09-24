#include "viterbiDecoder.h"

#include <cassert>

#include "feedForwardTrellis.h"

namespace {

int hammingDistance(const std::vector<int> x, const std::vector<int>& y) {
  assert(x.size() == y.size());
  int distance = 0;
  for (int i = 0; i < x.size(); ++i) {
    distance += (x[i] != y[i]);
  }
  return distance;
}

}  // namespace

ViterbiDecoder::ViterbiDecoder(const FeedForwardTrellis& FFT) : trellis_(&FFT) {
  n_ = trellis_->n_;
  numStates_ = trellis_->numStates_;
}

void ViterbiDecoder::constructDecodeTrellis(const std::vector<int>& received) {
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
  int message_length = received.size();
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
      auto begin = received.begin() + cur_stage * n_;
      auto end = begin + n_;
      std::vector<int> target_message(begin, end);
      // activate the next states
      for (int i = 0; i < trellis_->nextStates_[cur_state].size(); ++i) {
        int next_state = trellis_->nextStates_[cur_state][i];
        trellis_states_[next_state][cur_stage + 1].init = true;

        int possible_output = trellis_->output_[cur_state][i];
        std::vector<int> expected_output =
            trellis_->helper.convertIntToBits(possible_output, n_);

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

std::vector<int> ViterbiDecoder::decode(const std::vector<int>& received) {
  // construct trace back trellis
  constructDecodeTrellis(received);
  
  std::vector<int> result;

  return result;
}
