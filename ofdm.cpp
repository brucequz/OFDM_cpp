#include "ofdm.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <Eigen/Dense>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

Ofdm::Ofdm(const std::map<std::string, int>& config) : generator(12345) {
  auto it = config.find("B");
  B_ = (it != config.end()) ? it->second : 10;

  it = config.find("L");
  L_ = (it != config.end()) ? it->second : 16;

  it = config.find("CP_length");
  CP_length_ = (it != config.end()) ? it->second : 4;

  it = config.find("Nh");
  Nh_ = (it != config.end()) ? it->second : 4;

  // Pilot Config
  it = config.find("Pilot");
  pilot_type_ = (it != config.end()) ? it->second : 0;

  if (pilot_type_ == 1) {
    pilot_tx_ = generateBlockPilotSymbol(L_);
  }

  // constellation
  constl_ = new Constellation();
  constellations_["bpsk"] = constl_->bpsk();
  std::cout << "Average Bit Energy for bpsk: "
            << constl_->calculateAverageBitEnergy(constellations_["bpsk"])
            << std::endl;
  constellations_["qpsk"] = constl_->qpsk();
  std::cout << "Average Bit Energy for qpsk: "
            << constl_->calculateAverageBitEnergy(constellations_["qpsk"])
            << std::endl;
  constellations_["qam16"] = constl_->qam(16);
  std::cout << "Average Bit Energy for 16qam: "
            << constl_->calculateAverageBitEnergy(constellations_["qam16"])
            << std::endl;
}

Ofdm::~Ofdm() {
  delete constl_;
  constl_ = nullptr;
}

std::vector<int> Ofdm::generateRandomInt(const std::string& constl_type) {
  // TODO: generate bits such that all modulation compare the number of
  // iterations
  //       when they make the same number of errors

  std::vector<int> randomIntegers;
  if (pilot_type_ == 0) {
    randomIntegers.reserve(B_ * L_);
    for (int i = 0; i < B_ * L_; ++i) {
      randomIntegers.push_back(std::rand() %
                               constellations_[constl_type].size());
    }
  } else if (pilot_type_ == 1) {
    randomIntegers.reserve((B_ - 1) * L_);
    for (int i = 0; i < (B_ - 1) * L_; ++i) {
      randomIntegers.push_back(std::rand() %
                               constellations_[constl_type].size());
    }
  }

  return randomIntegers;
}

std::vector<std::vector<std::complex<double>>> Ofdm::generateModulatedSignal(
    const std::vector<int>& integers, const std::string& constl_type) {
  std::vector<std::vector<std::complex<double>>> modulated_signal;
  std::vector<std::complex<double>> constellation =
      constellations_[constl_type];
  
  if (pilot_type_ == 0) {
    int total_elements = B_ * L_;
    if (integers.size() != total_elements) {
      std::cout << "Cannot generate modulated signal. Input vector size does not "
                  "match target size.\n";
      return {};
    }
    modulated_signal.resize(B_, std::vector<std::complex<double>>(L_));
    for (int i = 0; i < L_; ++i) {
      for (int j = 0; j < B_; ++j) {
        modulated_signal[j][i] = constellation[integers[i * B_ + j]];
      }
    }
    return modulated_signal;
  } else if (pilot_type_ == 1) {
    int total_elements = (B_ - 1) * L_;
    if (integers.size() != total_elements) {
      std::cout << "Cannot generate modulated signal. Input vector size does not "
                  "match target size.\n";
      return {};
    }
    modulated_signal.resize(B_, std::vector<std::complex<double>>(L_));
    for (int i = 0; i < L_; ++i) {
      modulated_signal[0][i] = pilot_tx_[i];
    }
    for (int i = 0; i < L_; ++i) {
      for (int j = 1; j < B_; ++j) {
        modulated_signal[j][i] = constellation[integers[i * (B_-1) + j - 1]];
      }
    }
    return modulated_signal;
  }
  return {};
}

std::vector<int> Ofdm::convertBits(int value, int num_bits) {
  std::vector<int> bits;
  for (int i = num_bits - 1; i >= 0; --i) {
    int bit = (value >> i) & 1;
    bits.push_back(bit);
  }
  return bits;
}

std::vector<int> Ofdm::convertIntToBits(const std::vector<int>& integers,
                                        const std::string& constl_type) {
  std::vector<int> bits;

  int num_bits_per_symbol = 0;
  if (constl_type == "bpsk") {
    num_bits_per_symbol = 1;
  } else if (constl_type == "qpsk") {
    num_bits_per_symbol = 2;
  } else if (constl_type == "qam16") {
    num_bits_per_symbol = 4;
  }

  for (int value : integers) {
    std::vector<int> value_bits = convertBits(value, num_bits_per_symbol);
    bits.insert(bits.end(), value_bits.begin(), value_bits.end());
  }

  return bits;
}

std::vector<std::vector<double>> Ofdm::transpose2DVector(
    const std::vector<std::vector<double>>& input) {
  int rows = input.size();
  int cols = input[0].size();

  std::vector<std::vector<double>> transposed(cols, std::vector<double>(rows));

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      transposed[j][i] = input[i][j];
    }
  }

  return transposed;
}

std::vector<std::vector<std::complex<double>>> Ofdm::transpose2DComplexVector(
    const std::vector<std::vector<std::complex<double>>>& input) {
  int rows = input.size();
  int cols = input[0].size();

  std::vector<std::vector<std::complex<double>>> transposed(
      cols, std::vector<std::complex<double>>(rows));

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      transposed[j][i] = input[i][j];
    }
  }

  return transposed;
}

std::vector<std::vector<std::complex<double>>> Ofdm::reshape(
    const std::vector<std::complex<double>>& input, int num_rows,
    int num_cols) {
  int total_elements = num_rows * num_cols;

  if (input.size() != total_elements) {
    std::cout
        << "Cannot reshape. Input vector size does not match target size.\n";
    return {};
  }

  std::vector<std::vector<std::complex<double>>> reshaped(
      num_rows, std::vector<std::complex<double>>(num_cols));

  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      reshaped[i][j] = input[i * num_cols + j];
    }
  }

  return reshaped;
}

std::vector<std::vector<std::complex<double>>> Ofdm::columnMajorReshape(
    const std::vector<std::complex<double>>& input, int num_rows,
    int num_cols) {
  int total_elements = num_rows * num_cols;

  if (input.size() != total_elements) {
    std::cout << "Cannot column-major reshape. Input vector size does not "
                 "match target size.\n";
    return {};
  }

  std::vector<std::vector<std::complex<double>>> reshaped(
      num_rows, std::vector<std::complex<double>>(num_cols));

  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      reshaped[i][j] = input[i + j * num_rows];
    }
  }

  return reshaped;
}

int Ofdm::symbolErrorCount(const std::vector<int>& vector1,
                           const std::vector<int>& vector2) {
  // Check if the vectors have the same size; if not, return -1 to indicate an
  // error.
  if (vector1.size() != vector2.size()) {
    return -1;
  }

  int errorCount = 0;

  // Iterate through the elements of both vectors and compare them.
  for (size_t i = 0; i < vector1.size(); ++i) {
    if (vector1[i] != vector2[i]) {
      // If elements are different, increment the error count.
      errorCount++;
    }
  }

  return errorCount;
}

double Ofdm::totalTxSymbols(const int& iterations) {
  if (pilot_type_ == 0) {
    return static_cast<double>(iterations * B_ * L_);
  } else if (pilot_type_ == 1) {
    return static_cast<double>(iterations * (B_-1) * L_);
  }
  return 0.0;
}

double Ofdm::totalTxBits(const int& iterations, const std::string& mod) {
  if (pilot_type_ == 0) {
    return static_cast<double>(log2(constellations_[mod].size()) * iterations * B_ * L_);
  } else if (pilot_type_ == 1) {
    return static_cast<double>(log2(constellations_[mod].size()) * iterations * (B_-1) * L_);
  }

  return 0.0;
}

std::vector<std::vector<std::complex<double>>> Ofdm::fft(
    const std::vector<std::vector<std::complex<double>>>& input) {
  /*
    normalization factor sqrt(L) is incorporated in ifft
    " / sqrt(static_cast<double>(L_)); "
  */
  int numRows = input.size();
  int numCols = input[0].size();

  // Allocate space for the FFTW input and output arrays
  fftw_complex* in = reinterpret_cast<fftw_complex*>(
      fftw_malloc(sizeof(fftw_complex) * numCols));
  fftw_complex* out = reinterpret_cast<fftw_complex*>(
      fftw_malloc(sizeof(fftw_complex) * numCols));

  // Create an FFTW plan for the 1D FFT
  fftw_plan plan =
      fftw_plan_dft_1d(numCols, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

  // Initialize the output vector
  std::vector<std::vector<std::complex<double>>> output(
      numRows, std::vector<std::complex<double>>(numCols));

  // Perform FFT along rows
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      in[j][0] = input[i][j].real();
      in[j][1] = input[i][j].imag();
    }
    fftw_execute(plan);
    for (int j = 0; j < numCols; ++j) {
      output[i][j] = std::complex<double>(out[j][0], out[j][1]) /
                     std::sqrt(static_cast<double>(L_));
    }
  }

  // Clean up FFTW resources
  fftw_destroy_plan(plan);
  fftw_free(in);
  fftw_free(out);

  return output;
}

std::vector<std::complex<double>> Ofdm::fft(
    const std::vector<std::complex<double>>& input, int n) {
  // Validate input size
  if (n <= 0) {
    throw std::invalid_argument("Output size n must be positive");
  }

  int input_size = input.size();

  // Prepare result vector
  std::vector<std::complex<double>> result(n);

  // Create a plan for the forward FFT
  fftw_plan plan_forward =
      fftw_plan_dft_1d(n,
                       reinterpret_cast<fftw_complex*>(
                           const_cast<std::complex<double>*>(result.data())),
                       reinterpret_cast<fftw_complex*>(
                           const_cast<std::complex<double>*>(result.data())),
                       FFTW_FORWARD, FFTW_ESTIMATE);

  // Copy input, pad with zeros or truncate as needed
  if (input_size < n) {
    std::copy(input.begin(), input.end(), result.begin());
    std::fill(result.begin() + input_size, result.end(),
              std::complex<double>(0.0, 0.0));
  } else {
    std::copy(input.begin(), input.begin() + n, result.begin());
  }

  // Perform the forward FFT
  fftw_execute(plan_forward);

  // Clean up the plan
  fftw_destroy_plan(plan_forward);

  return result;
}

std::vector<std::vector<std::complex<double>>> Ofdm::ifft(
    const std::vector<std::vector<std::complex<double>>>& input) {
  /*
    normalization factor sqrt(L) is incorporated in ifft
    " * sqrt(static_cast<double>(L_)); "
  */
  int num_rows = input.size();
  int num_cols = input[0].size();

  fftw_plan plan = fftw_plan_dft_1d(num_cols, nullptr, nullptr, FFTW_BACKWARD,
                                    FFTW_ESTIMATE);

  std::vector<std::vector<std::complex<double>>> output(
      num_rows, std::vector<std::complex<double>>(num_cols));

  for (int i = 0; i < num_rows; ++i) {
    std::vector<std::complex<double>> row_input(num_cols);
    std::vector<std::complex<double>> row_output(num_cols);

    for (int j = 0; j < num_cols; ++j) {
      row_input[j] = input[i][j];
    }

    fftw_complex* in = reinterpret_cast<fftw_complex*>(&row_input[0]);
    fftw_complex* out = reinterpret_cast<fftw_complex*>(&row_output[0]);

    fftw_execute_dft(plan, in, out);

    for (int j = 0; j < num_cols; ++j) {
      output[i][j] = row_output[j] / static_cast<double>(num_cols) *
                     sqrt(static_cast<double>(L_));
    }
  }

  fftw_destroy_plan(plan);

  return output;
}

std::vector<std::vector<std::complex<double>>> Ofdm::addCyclicPrefix(
    std::vector<std::vector<std::complex<double>>> vec2D) {
  int numRows = vec2D.size();
  int numCols = vec2D[0].size();

  if (CP_length_ >= numCols) {
    std::cerr << "Number of columns to copy is greater than or equal to the "
                 "total number of columns."
              << std::endl;
    return {};
  }

  for (int i = 0; i < numRows; ++i) {
    std::vector<std::complex<double>> newCols(vec2D[i].end() - CP_length_,
                                              vec2D[i].end());
    vec2D[i].insert(vec2D[i].begin(), newCols.begin(), newCols.end());
  }
  return vec2D;
}

std::vector<std::vector<std::complex<double>>> Ofdm::removeCyclicPrefix(
    std::vector<std::complex<double>> input) {
  std::vector<std::vector<std::complex<double>>> rec_reshaped =
      transpose2DComplexVector(columnMajorReshape(input, CP_length_ + L_, B_));

  std::vector<std::vector<std::complex<double>>> rec_sans_cp;
  for (const auto& row : rec_reshaped) {
    int num_cols_original = row.size();
    int num_cols_remaining = num_cols_original - CP_length_;

    if (num_cols_remaining > 0) {
      std::vector<std::complex<double>> truncated_row(row.begin() + CP_length_,
                                                      row.end());
      rec_sans_cp.push_back(truncated_row);
    }
  }
  return rec_sans_cp;
}

std::vector<std::complex<double>> Ofdm::addAWGN(
    const std::vector<std::complex<double>>& signal) {
  std::default_random_engine generator;
  std::normal_distribution<double> distribution(0.0, 1.0);

  std::vector<std::complex<double>> noisy_signal;
  noisy_signal.reserve(signal.size());

  for (const std::complex<double>& value : signal) {
    double noise_real = distribution(generator) * std::sqrt(0.5 / L_);
    double noise_imag = distribution(generator) * std::sqrt(0.5 / L_);

    std::complex<double> noisy_value(value.real() + noise_real,
                                     value.imag() + noise_imag);
    noisy_signal.push_back(noisy_value);
  }

  return noisy_signal;
}

std::vector<std::complex<double>> Ofdm::generateNoise(const double& mean,
                                                      const double& stddev,
                                                      size_t size) {
  std::normal_distribution<double> distribution(mean, stddev);

  std::vector<std::complex<double>> noise;
  noise.reserve(size);

  for (int i = 0; i < size; i++) {
    double noise_real = distribution(generator);
    double noise_imag = distribution(generator);

    std::complex<double> noise_value(noise_real, noise_imag);
    noise.push_back(noise_value);
  }

  return noise;
}

std::vector<std::complex<double>> Ofdm::filter(
    const std::vector<std::complex<double>>& signal,
    const std::vector<std::complex<double>>& filter_coeffs) {
  int signal_size = signal.size();
  int filter_size = filter_coeffs.size();

  std::vector<std::complex<double>> result(signal_size + filter_size - 1,
                                           std::complex<double>(0.0, 0.0));

  for (int i = 0; i < signal_size; ++i) {
    for (int j = 0; j < filter_size; ++j) {
      result[i + j] += signal[i] * filter_coeffs[j];
    }
  }

  for (int i = 0; i < CP_length_; ++i) {
    result.pop_back();
  }

  return result;
}

std::vector<std::vector<int>> Ofdm::decode(
    std::vector<std::vector<std::complex<double>>> received,
    std::vector<std::complex<double>> channel_response,
    const std::string& constl_type, const std::string& compensation) {
  
  if (pilot_type_ == 1) {
    // remove pilot symbol before decoding
    received.erase(received.begin());
  }

  if (compensation == "zero") {
    return decodeNoCompensation(received, constl_type);
  } else if (compensation == "full") {
    return decodeFullCompensation(received, channel_response, constl_type);     
  }
  
  return decodeNoCompensation(received, constl_type);
}

std::vector<std::vector<int>> Ofdm::decodeNoCompensation(
    std::vector<std::vector<std::complex<double>>>& received,
    const std::string& constl_type) {
  /**
   *  @brief  Decode the input B x L matrix of received symbols
   *
   *  @param  received A BxL matrix of received symbols
   *  @param  constl_type a string indicating constellation type
   *
   *  @result returns a 2d vector of decoded symbols
   */

  std::vector<std::vector<int>> dec_sym;

  std::vector<std::complex<double>> constellation =
      constellations_[constl_type];
  
  for (std::vector<std::complex<double>> rec_symbol : received) {
    std::vector<std::vector<double>> det;
    for (int i = 0; i < constellation.size(); i++) {
      std::vector<double> det_ind;
      for (int j = 0; j < rec_symbol.size(); j++) {
        det_ind.push_back(
            squareEuclideanDistance(rec_symbol[j], constellation[i]));
      }
      det.push_back(det_ind);
    }
    std::vector<std::vector<double>> det_transpose = transpose2DVector(det);
    dec_sym.push_back(findMinInd(det_transpose));
  }

  return dec_sym;
}

std::vector<std::vector<int>> Ofdm::decodeFullCompensation(
    std::vector<std::vector<std::complex<double>>>& received,
    std::vector<std::complex<double>> channel_response,
    const std::string& constl_type) {
  /**
   *  @brief  Decode the input B x L matrix of received symbols
   *
   *  @param  received A BxL matrix of received symbols
   *  @param  channel_response channel impulse response
   *  @param  constl_type a string indicating constellation type
   *
   *  @result returns a 2d vector of decoded symbols
   */

  std::vector<std::vector<int>> dec_sym;

  std::vector<std::complex<double>> constellation =
      constellations_[constl_type];

  for (std::vector<std::complex<double>> rec_symbol : received) {
    std::vector<std::vector<double>> det;
    for (int i = 0; i < constellation.size(); i++) {
      std::vector<double> det_ind;
      for (int j = 0; j < channel_response.size(); j++) {
        det_ind.push_back(squareEuclideanDistance(
            complexDivision(rec_symbol[j], channel_response[j]),
            constellation[i]));
      }
      det.push_back(det_ind);
    }
    std::vector<std::vector<double>> det_transpose = transpose2DVector(det);
    dec_sym.push_back(findMinInd(det_transpose));
  }

  return dec_sym;
}

std::vector<std::complex<double>> Ofdm::generateBlockPilotSymbol(
    const int& pilot_length) {
  // Define the known pilot symbols for 16-QAM (4 possible symbols)
  std::vector<std::complex<double>> pilotSymbols = {
      {1.0, 1.0}, {-1.0, 1.0}, {1.0, -1.0}, {-1.0, -1.0}};

  // Create a block pilot symbol by repeating the known symbols for each
  // subcarrier
  std::vector<std::complex<double>> blockPilotSymbol;
  for (int i = 0; i < pilot_length; ++i) {
    blockPilotSymbol.push_back(pilotSymbols[i % 4]);
  }

  return blockPilotSymbol;
}

std::vector<std::complex<double>> Ofdm::estimateChannelML(
    const std::vector<std::complex<double>>& pilot_rx) {
  /*
  Returns a estimation of the Channel impulse response (CIR).

  */
  std::vector<std::complex<double>> h_hat;
  for (int i = 0; i < pilot_rx.size(); ++i) {
    h_hat.push_back(complexDivision(pilot_rx[i], pilot_tx_[i]));
  }
  return h_hat;
}

std::vector<std::complex<double>> Ofdm::estimateChannelLS(
    const std::vector<std::complex<double>>& pilot_rx) {
  /* Least squares estimation of channel frequency response H_hat
  */
  Eigen::VectorXcd h_hat_eigen;
  // Convert std::vector to Eigen::complex vector
  Eigen::VectorXcd diagonalVector = Eigen::Map<const Eigen::VectorXcd>(pilot_tx_.data(), pilot_tx_.size());
  // diagonalize the Eigen::vector
  Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic> diagonalMatrix(diagonalVector);
  // Take the inverse 
  Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic> inverseDiagonalMatrix = diagonalMatrix.inverse();

  Eigen::VectorXcd eigenComplexVectorRx = Eigen::Map<const Eigen::VectorXcd>(pilot_rx.data(), pilot_rx.size());

  // Multiply inverseDiagonalMatrix with pilot_rx_
  h_hat_eigen = inverseDiagonalMatrix * eigenComplexVectorRx;

  // Convert h_hat_eigen to a std::vector of 1D
  std::vector<std::complex<double>> H_hat(h_hat_eigen.data(), h_hat_eigen.data() + h_hat_eigen.size());

  return H_hat;
}

std::vector<std::complex<double>> Ofdm::estimateChannelMMSE(
      const std::vector<std::complex<double>>& pilot_rx) {
  /* Minimum Mean Squared Error estimation of channel frequency response H_hat
  */
  
  // data read-in
  std::vector<std::complex<double>> H; 
  for (size_t i = 0; i < pilot_rx.size(); ++i) {
    H.push_back(complexDivision(pilot_rx[i], pilot_tx_[i]));
  }
  std::vector<std::complex<double>> h_std = ifft({H})[0];

  Eigen::VectorXcd h = Eigen::Map<const Eigen::VectorXcd>(h_std.data(), h_std.size());
  Eigen::VectorXcd x = Eigen::Map<const Eigen::VectorXcd>(pilot_tx_.data(), pilot_tx_.size());
  Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> X(x);
  Eigen::DiagonalMatrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> X_H(x.conjugate());
  Eigen::VectorXcd y = Eigen::Map<const Eigen::VectorXcd>(pilot_rx.data(), pilot_rx.size());
  
  // divide to get h
  // Eigen::VectorXcd h = y.array() / x.array();
  std::complex<double> h_mean = h.mean();

  Eigen::MatrixXcd R_hh(h.size(), h.size());
  for (int i = 0; i < h.size(); ++i) {
      for (int j = 0; j < h.size(); ++j) {
          R_hh(i, j) = (h[i] - h_mean) * std::conj(h[j] - h_mean);
      }
  }
  
  Eigen::MatrixXcd F(L_, L_);
  for (int i = 0; i < L_; ++i) {
    for (int j = 0; j < L_; ++j) {
        std::complex<double> exponent(0.0, -2.0*M_PI*i*j/static_cast<double>(L_));
        F(i, j) = 1/static_cast<double>(L_) * std::exp(exponent);
    }
  }
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> F_H = F.adjoint();

  Eigen::MatrixXcd R_hY = R_hh * F_H * X_H;

  Eigen::MatrixXcd R_YY(h.size(), h.size());
  R_YY = X * F * R_hh * F_H * X_H + 0.5/L_ * Eigen::MatrixXd::Identity(x.size(), x.size());

  Eigen::VectorXcd H_est =  F * R_hY * R_YY.inverse() * y;

  std::vector<std::complex<double>> H_est_vector(H_est.data(), H_est.data() + H_est.size());
  return H_est_vector;
}

double Ofdm::squareEuclideanDistance(const std::complex<double>& a,
                                     const std::complex<double>& b) {
  double realDiff = std::real(a) - std::real(b);
  double imagDiff = std::imag(a) - std::imag(b);
  return realDiff * realDiff + imagDiff * imagDiff;
}

std::complex<double> Ofdm::complexDivision(const std::complex<double>& a,
                                           const std::complex<double>& b) {
  double denominator = b.real() * b.real() + b.imag() * b.imag();
  std::complex<double> result = a * std::conj(b) / denominator;
  return result;
}

std::vector<int> Ofdm::findMinInd(
    const std::vector<std::vector<double>>& matrix) {
  std::vector<int> minIndices;

  for (const auto& row : matrix) {
    if (row.empty()) {
      // Handle empty row, return -1 or some other sentinel value
      minIndices.push_back(-1);
    } else {
      double minVal = row[0];
      int minIndex = 0;

      for (int i = 1; i < row.size(); ++i) {
        if (row[i] < minVal) {
          minVal = row[i];
          minIndex = i;
        }
      }

      minIndices.push_back(minIndex);
    }
  }

  return minIndices;
}

std::complex<double> Ofdm::calculateMean(
    const std::vector<std::complex<double>>& input) {
  if (input.empty()) {
    // Handle the case of an empty vector
    return std::complex<double>(0.0, 0.0);
  }

  std::complex<double> sum(0.0, 0.0);

  for (const std::complex<double>& element : input) {
    sum += element;
  }

  return sum / static_cast<double>(input.size());
}

std::complex<double> Ofdm::calculateStandardDeviation(
    const std::vector<std::complex<double>>& input) {
  if (input.empty()) {
    // Handle the case of an empty vector
    return std::complex<double>(0.0, 0.0);
  }

  // Step 1: Calculate the mean
  std::complex<double> mean = calculateMean(input);

  // Step 2: Calculate the squared differences
  std::complex<double> squaredDifferencesSum(0.0, 0.0);
  for (const std::complex<double>& element : input) {
    std::complex<double> diff = element - mean;
    squaredDifferencesSum +=
        std::norm(diff);  // Square and sum real and imaginary parts
  }

  // Step 3: Calculate the mean of squared differences
  double meanSquaredDifferences =
      std::abs(squaredDifferencesSum) / static_cast<double>(input.size());

  // Step 4: Take the square root to get the standard deviation
  return std::sqrt(meanSquaredDifferences);
}