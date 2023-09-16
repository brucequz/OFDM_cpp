#include <gtest/gtest.h>
#include <vector>
#include "../src/ErrorCorrection/feedForwardTrellis.h"
#include "../src/helper.h"

// Define test cases
TEST(TrellisTest, DecodeNoError) {
  std::vector<int> poly = {13, 17};
  FeedForwardTrellis trellis(1, 2, 3, poly);
  std::vector<int> input_1 = {1, 0, 1, 0, 1, 1, 1};
  std::vector<int> encoded = trellis.encode(input_1);
  std::vector<int> expect_encoded = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1};
  std::vector<int> input_2 = {0, 0, 1, 0, 1, 1, 0};
  std::vector<int> encoded_2 = trellis.encode(input_2);
  std::vector<int> expect_encoded_2 = {0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0};
  
  EXPECT_EQ(encoded,expect_encoded);
  EXPECT_EQ(encoded_2,expect_encoded_2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}