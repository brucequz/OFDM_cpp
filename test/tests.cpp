#include <gtest/gtest.h>
#include "../src/ErrorCorrection/hammingCode.h"  // Include the header for your HammingCode class

// Define test cases
TEST(HammingCodeTest, DecodeNoError) {
    // Create an instance of your HammingCode class
    HammingCode hammingCode;

    // Define a received codeword with no errors
    Eigen::VectorXi receivedCodeword(7);
    receivedCodeword << 1, 0, 1, 0, 1, 1, 0;

    // Perform decoding
    Eigen::VectorXi decodedMessage = hammingCode.decodeHammingCode(receivedCodeword);

    // Check if the decoded message matches the expected result
    Eigen::VectorXi expectedMessage(4);
    expectedMessage << 1, 0, 1, 0;
    Eigen::VectorXi Message(4);
    Message << 0, 0, 1, 0;
    EXPECT_EQ(decodedMessage, expectedMessage);
    EXPECT_NE(decodedMessage, Message);
}

TEST(HammingCodeTest, DecodeWithError) {
    // Create an instance of your HammingCode class
    HammingCode hammingCode;

    // Define a received codeword with a single-bit error
    Eigen::VectorXi receivedCodeword(7);
    receivedCodeword << 1, 0, 1, 0, 1, 1, 1;

    // Perform decoding
    Eigen::VectorXi decodedMessage = hammingCode.decodeHammingCode(receivedCodeword);

    // Check if the decoded message matches the expected result after error correction
    Eigen::VectorXi expectedMessage(4);
    expectedMessage << 1, 0, 1, 0;
    Eigen::VectorXi Message(4);
    Message << 1, 0, 0, 0;

    EXPECT_NE(decodedMessage, expectedMessage);
    EXPECT_EQ(decodedMessage, Message);
}

// Add more test cases as needed

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
