#include <gtest/gtest.h>
#include "../src/ErrorCorrection/minHeap.h"

TEST(MinHeapTest, HeapOperation) {
  MinHeap heap;
  DetourNode node1;
  node1.path_metric = 5;
  DetourNode node2;
  node2.path_metric = 3;
  DetourNode node3;
  node3.path_metric = 1;
  DetourNode node4;
  node4.path_metric = 10;
  DetourNode node5;
  node5.path_metric = 7;

  heap.insert(node1);
  heap.insert(node2);
  heap.insert(node3);
  heap.insert(node4);
  heap.insert(node5);

  EXPECT_EQ(heap.size(), 5);
  EXPECT_EQ(heap.pop().path_metric, 1);
  EXPECT_EQ(heap.size(), 4);
  EXPECT_EQ(heap.pop().path_metric, 3);
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}