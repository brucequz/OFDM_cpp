#ifndef CONVOLUTIONALCODE_H
#define CONVOLUTIONALCODE_H

#include <vector>

class ConvolutionalCode {
public:
    // Constructor
    ConvolutionalCode(int constraintLength, const std::vector<int>& generatorPolynomials);

    // Destructor
    ~ConvolutionalCode();

    // Encode a block of input bits
    std::vector<int> encode(const std::vector<int>& inputBits);

    // Decode a block of received bits using Viterbi decoding
    std::vector<int> decode(const std::vector<int>& receivedBits);

private:
    // Private members and functions for internal use
    int constraintLength_;
    std::vector<int> generatorPolynomials_;
};

#endif // CONVOLUTIONALCODE_H
