#include "minHeap.h"

MinHeap::MinHeap() {}

void MinHeap::insert(DetourNode node) {
  /*
   * append it to the end, and then swap with parent if necessary.
   */
  heap_.push_back(node); 
  int index = heap_.size() - 1;
  while (index > 0 && heap_[parentIndex(index)] > node) {
    std::swap(heap_[index], heap_[parentIndex(index)]);
    index = parentIndex(index);
  }
}

DetourNode MinHeap::pop() {
  DetourNode min_node = heap_.front();
  heap_.pop_front();
  heapify(0);
  return min_node;
}

void MinHeap::heapify(int index) {
  int left_index = leftChildIndex(index);
  int right_index = rightChildIndex(index);

  int min_index = index;
  if (left_index < heap_.size() &&
      heap_[leftChildIndex(index)] < heap_[index]) {
    min_index = left_index;
  }
  if (right_index < heap_.size() &&
      heap_[rightChildIndex(index)] < heap_[index]) {
    min_index = right_index;
  }
  if (min_index != index) {
    std::swap(heap_[min_index], heap_[index]);
    heapify(min_index);
  }
}