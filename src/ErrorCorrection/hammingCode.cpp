#include "hammingCode.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// Constructor
HammingCode::HammingCode() {
    // TODO: Initialize any class members if needed
    generatorMatrix << 1, 0, 0, 0, 1, 1, 1,
                       0, 1, 0, 0, 1, 1, 0,
                       0, 0, 1, 0, 1, 0, 1,
                       0, 0, 0, 1, 0, 1, 1;
    parityCheckMatrix << 1, 1, 1, 0, 1, 0, 0,
                         1, 1, 0, 1, 0, 1, 0,
                         1, 0, 1, 1, 0, 0, 1;
    // TODO: initialize k_ and v_
    k_ = generatorMatrix.rows();
    v_ = generatorMatrix.cols();
}

// Destructor
HammingCode::~HammingCode() {
    // TODO: Clean up any allocated resources if needed
}

// Function to generate Hamming code from input data

Eigen::VectorXi HammingCode::encodeHammingCode(const Eigen::VectorXi& data) {
  // Calculate the number of chunks (blocks) required to encode the input
  int numChunks = (data.size() + k_ - 1) / k_;

  // Calculate the size of the padded chunk
  int paddedChunkSize = k_;

  // Initialize the encoded data vector with appropriate size
  Eigen::VectorXi encodedData(numChunks * v_);

  for (int chunkIdx = 0; chunkIdx < numChunks; ++chunkIdx) {
      // Calculate the start and end indices for the current chunk
      int startIdx = chunkIdx * k_;
      int endIdx = std::min<int>((chunkIdx + 1) * k_, data.size());

      // Extract the current chunk of input data
      Eigen::VectorXi chunk = data.segment(startIdx, endIdx - startIdx);

      // Pad the current chunk with zeros if it's shorter than k_
      if (chunk.size() < paddedChunkSize) {
          Eigen::VectorXi paddedChunk(paddedChunkSize);
          paddedChunk.head(chunk.size()) = chunk;
          paddedChunk.tail(paddedChunkSize - chunk.size()).setZero();
          chunk = paddedChunk;
      }

      // Encode the current chunk and append to the output vector
      Eigen::VectorXi encodedChunk = (chunk.transpose() * generatorMatrix).cast<int>();
      for (int i = 0; i < encodedChunk.size(); ++i) {
          encodedData[chunkIdx * v_ + i] = encodedChunk[i] % 2;
      }
  }

  return encodedData;
}

// Function to decode received Hamming code
Eigen::VectorXi HammingCode::decodeHammingCode(Eigen::VectorXi received) {
  Eigen::VectorXi syndrome = parityCheckMatrix * received;
  for (int i = 0; i < syndrome.size(); ++i) {
    syndrome[i] %= 2;
  }
  Eigen::VectorXi decodedMessage;
  // std::cout << "syndrome is: " << syndrome << std::endl;
  // std::cout << "syndrome is zero: " << syndrome.isZero() << std::endl;

  // Check for errors (non-zero syndrome)
  if (syndrome.isZero()) {
      // Extract the original message bits (first 4 bits of the received codeword)
      decodedMessage = received.head(4);
  } else {
      // Error correction: Find the bit position of the error and correct it
      int errorPosition = 0;
      for (int i = 0; i < 7; ++i) {
          if (syndrome == parityCheckMatrix.col(i)) {
              errorPosition = i;
              break;
          }
      }
      // std::cout << "error position: " << errorPosition << std::endl;

      if (errorPosition >= 0) {

          // Correct the error by flipping the corresponding bit
          // std:: cout << "original received: " << received.transpose() << std::endl;
          received(errorPosition) ^= 1;
          // std:: cout << "corrected received: " << received.transpose() << std::endl;

          // Extract the corrected message bits (first 4 bits of the received codeword)
          decodedMessage = received.head(4);
      } else {
          // std::cout << "Error position not found. Unable to correct." << std::endl;
      }
  }
  // std::cout << "decoded again: " << decodedMessage.transpose() << std::endl;
  return decodedMessage;
}
