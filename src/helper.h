#ifndef HELPER_H
#define HELPER_H

#include <vector>
#include <iostream>

class Helper{
  public: 
    Helper(){};
    ~Helper(){};
    template <typename T>
    void print(const std::vector<T>& vec) {
      for (const T& element : vec) {
          std::cout << element << " ";
      }
      std::cout << std::endl;
    }

    template <typename T>
    void print(const std::vector<std::vector<T>>& matrix) {
      for (const std::vector<T>& row : matrix) {
          for (const T& element : row) {
              std::cout << element << " ";
          }
          std::cout << std::endl;
    }
    }

    std::vector<int> convertIntToBits(int integer, int length) const {
      if (integer < 0) {
        std::cerr << "CANNOT CONVERT: negative integer" << std::endl;
      } else if (std::ceil(std::log2(integer + 1)) > length) {
        std::cerr << "CANNOT CONVERT: integer too large" << std::endl;
      }
      std::vector<int> result(length, 0);
      int i = length - 1;
      while (integer > 0 && i >= 0) {
        int remainder = integer % 2;
        result[i] = remainder;
        integer /= 2;
        i--;
      }
      return result;
    }

    int countSetBits(int num) const {
      int count = 0;
      while (num > 0) {
          count += num & 1;
          num >>= 1;
      }
      return count;
    }

    int countDifferentBits(const std::vector<int>& vector1, const std::vector<int>& vector2) const {
      // Ensure both vectors have the same size (handle size mismatch as needed)
      if (vector1.size() != vector2.size()) {
          std::cerr << "Vector size mismatch!" << std::endl;
          return -1; // Return an error code or handle the error as appropriate
      }

      int different_bit_count = 0;

      for (size_t i = 0; i < vector1.size(); ++i) {
          // XOR the corresponding elements and count the set bits
          int difference = vector1[i] ^ vector2[i];
          different_bit_count += countSetBits(difference);
      }

      return different_bit_count;
    }
};

#endif