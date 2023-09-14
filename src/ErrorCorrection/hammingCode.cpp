#include "hammingCode.h"
#include <iostream>

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
}

// Destructor
HammingCode::~HammingCode() {
    // TODO: Clean up any allocated resources if needed
}

// Function to generate Hamming code from input data

Eigen::VectorXi HammingCode::encodeHammingCode(const Eigen::VectorXi& data) {
    // Ensure that the input data has a size of 4 (information bits)
    assert(data.size() == 4);

    Eigen::VectorXi encodedData = (data.transpose() * generatorMatrix).cast<int>();

    for (int i = 0; i < encodedData.size(); ++i) {
        encodedData[i] %= 2;
    }

    return encodedData;
}

// Function to decode received Hamming code
Eigen::VectorXi HammingCode::decodeHammingCode(Eigen::VectorXi received) {
  Eigen::VectorXi syndrome = parityCheckMatrix * received;
  Eigen::VectorXi decodedMessage;

  // Check for errors (non-zero syndrome)
  if (syndrome.isZero()) {
      // Extract the original message bits (first 4 bits of the received codeword)
      decodedMessage = received.head(4);
  } else {
      for (int i = 0; i < syndrome.size(); ++i) {
          syndrome[i] %= 2;
      }
      std::cout << "syndrome is: " << syndrome << std::endl;
      // Error correction: Find the bit position of the error and correct it
      int errorPosition = 0;
      for (int i = 0; i < 7; ++i) {
          if (syndrome == parityCheckMatrix.col(i)) {
              errorPosition = i;
              break;
          }
      }
      std::cout << "error position: " << errorPosition << std::endl;

      if (errorPosition >= 0) {

          // Correct the error by flipping the corresponding bit
          std:: cout << "original received: " << received.transpose() << std::endl;
          received(errorPosition) ^= 1;
          std:: cout << "corrected received: " << received.transpose() << std::endl;

          // Extract the corrected message bits (first 4 bits of the received codeword)
          decodedMessage = received.head(4);
      } else {
          std::cout << "Error position not found. Unable to correct." << std::endl;
      }
  }
  std::cout << "decoded again: " << decodedMessage.transpose() << std::endl;
  return decodedMessage;
}
