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
};





#endif