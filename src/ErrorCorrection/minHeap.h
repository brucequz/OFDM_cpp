#ifndef MINHEAP_H
#define MINHEAP_H

#include <deque>  // Include deque for the MinHeap implementation
#include <vector>

struct DetourNode {
  DetourNode(){};
  int path_metric;
  bool operator>(const DetourNode& other) const {
    return this->path_metric > other.path_metric;
  }
  bool operator<(const DetourNode& other) const {
    return this->path_metric < other.path_metric;
  }
};

class MinHeap {
 public:
  MinHeap();
  void insert(DetourNode node);
  DetourNode pop();
  int size() { return heap_.size(); }

 private:
  std::deque<DetourNode> heap_;
  void heapify(int index);
  int leftChildIndex(int index) { return 2 * index + 1; }
  int rightChildIndex(int index) { return 2 * index + 2; }
  int parentIndex(int index) { return (index - 1) / 2; }
};

#endif
