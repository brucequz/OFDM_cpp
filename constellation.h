#ifndef CONSTELLATION_H
#define CONSTELLATION_H

#include <complex>
#include <vector>

class Constellation {
 public:
  Constellation();
  std::vector<std::complex<double>> bpsk();
  std::vector<std::complex<double>> qpsk();
  std::vector<std::complex<double>> qam(int levels);

};

#endif  // CONSTELLATION_H
