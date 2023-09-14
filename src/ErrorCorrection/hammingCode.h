#ifndef HAMMINGCODE_H
#define HAMMINGCODE_H

#include "../../include/eigen/Dense"

class HammingCode {
public:
    // Constructor
    HammingCode();

    // Destructor
    ~HammingCode();

    // Function to generate Hamming code from input data
    Eigen::VectorXi encodeHammingCode(const Eigen::VectorXi& data);

    // Function to decode received Hamming code
    Eigen::VectorXi decodeHammingCode(Eigen::VectorXi received);

private:
    // TODO: Add any private members or helper functions if needed
    Eigen::Matrix<int, 4, 7> generatorMatrix;
    Eigen::Matrix<int, 3, 7> parityCheckMatrix;
};

#endif  // HAMMINGCODE_H
