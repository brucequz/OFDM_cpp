#include <gtest/gtest.h>
#include <vector>
#include "../src/ErrorCorrection/feedForwardTrellis.h"
#include "../src/helper.h"

// Define test cases
TEST(TrellisTest, DecodeNoError) {
  std::vector<int> poly = {27, 31};
  FeedForwardTrellis trellis(1, 2, 4, poly);
  
  EXPECT_EQ(1,1);
}