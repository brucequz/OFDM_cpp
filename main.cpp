#include <mat.h>

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "constellation.h"
#include "ofdm.h"

void print1DVector(const std::vector<int>& vec);
void output1DVector(std::ofstream& outputFile, const std::vector<int>& vec);
void printComplexVector2D(
    const std::vector<std::vector<std::complex<double>>>& vec2D);
void output1DComplexVector(std::ofstream& outputFile,
                           const std::vector<std::complex<double>>& vec);
void output2DComplexVector(
    std::ofstream& outputFile,
    const std::vector<std::vector<std::complex<double>>>& vec2D);
std::vector<std::complex<double>> normalize(
    const std::vector<std::complex<double>>& input);
std::vector<int> readVectorFromMATFile(const std::string& filePath,
                                       const std::string& variableName);
std::vector<std::vector<std::complex<double>>> readComplexMatFile(
    const std::string& filename, const std::string& variableName);
std::vector<std::complex<double>> addComplexVectors(
    const std::vector<std::complex<double>>& vector1,
    const std::vector<std::complex<double>>& vector2);

int main() {
  // output path
  std::ofstream outputFile("output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return 1;
  }
  // complex number
  const std::complex<double> i(0.0, 1.0);
  std::cout << i.imag() << std::endl;
  std::cout << 5 * i.imag() << std::endl;

  // OFDM
  std::map<std::string, int> config;
  config["B"] = 10;         // number of OFDM symbol per frame
  config["L"] = 16;         // number of subcarriers per ofdm symbol
  config["CP_length"] = 4;  // cyclic prefix length
  config["Nh"] = 4;         // channel order

  Ofdm ofdm(config);

  // ------------------------------- Execution Cycle
  // ----------------------------------- Symbol Error Rate for different
  // modulation
  int mc_N = 500;  // maximum number of iterations to achieve sufficient errors

  // SNR
  double SNR_dB_start = 0.0;
  double SNR_dB_end = 20.0;
  double SNR_dB_step = 2.0;
  std::vector<double> SNR_dB = {};
  std::vector<double> SNR = {};
  for (double i = SNR_dB_start; i <= SNR_dB_end; i += SNR_dB_step) {
    SNR_dB.push_back(i);
    SNR.push_back(pow(10.0, i / 10.0));
  }

  // Error Vector
  std::vector<double> Pe(SNR.size(), 0);

  // Modulation Schemes
  std::vector<std::string> modulation_schemes = {"bpsk", "qpsk", "qam16"};

  for (const std::string& mod : modulation_schemes) {
    outputFile << "For modulation: " << mod << ": --------------------------------------"
               << std::endl;
    for (int i = 0; i < SNR.size(); i++) {
      double rho = SNR[i];
      outputFile << "For SNR = " << rho << ": " << std::endl;
      int biterr_cnt = 0;
      int symerr_cnt = 0;
      for (int mc_loop = 0; mc_loop < mc_N; mc_loop++) {
        std::vector<int> data_int = ofdm.generateRandomInt(mod);

        // std::string filename = "../dat_ind.mat";
        // data_int = readVectorFromMATFile(filename, "dat_ind_reshape");
        // outputFile << "Outputing integer vector" << std::endl;
        // output1DVector(outputFile, data_int);

        std::vector<int> data_bits = ofdm.convertIntToBits(data_int, mod);
        // outputFile << "Outputing bits vector" << std::endl;
        // output1DVector(outputFile, data_bits);

        std::vector<std::vector<std::complex<double>>> data_sym =
            ofdm.generateModulatedSignal(data_int, mod);
        std::vector<std::complex<double>> data_flattened =
            ofdm.columnMajorFlatten(data_sym);
        // outputFile << "Outputing modulated signal vector" << std::endl;
        // output1DComplexVector(outputFile, data_flattened);

        // IFFT
        std::vector<std::vector<std::complex<double>>> data_ifft =
            ofdm.ifft(data_sym);
        // outputFile << "Outputing ifft vector" << std::endl;
        // output2DComplexVector(outputFile, data_ifft);

        // cyclic prefix
        std::vector<std::vector<std::complex<double>>> data_cp =
            ofdm.addCyclicPrefix(data_ifft);
        // outputFile << "Outputing cp vector" << std::endl;
        // output2DComplexVector(outputFile, data_cp);

        // Reshape the BxN matrix to obtain the frame (1xTotal_length)
        // Total_length = (CP_length+L)*B
        std::vector<std::complex<double>> data_tx =
            ofdm.columnMajorFlatten(ofdm.transpose2DComplexVector(data_cp));

        // Define the channel response
        // deterministic frequency-selective channel
        std::complex<double> h1(3.0, 0.0);
        std::complex<double> h2(-1.0 *
                                std::exp(std::complex<double>(0.0, 0.13)));
        std::complex<double> h3(1.0 *
                                std::exp(std::complex<double>(0.0, -0.35)));
        std::complex<double> h4(0.0, 0.0);
        std::complex<double> h5(4.0 *
                                std::exp(std::complex<double>(0.0, 1.03)));

        // Create a vector of complex numbers
        std::vector<std::complex<double>> h = {h1, h2, h3, h4, h5};
        std::vector<std::complex<double>> h_normalized = normalize(h);

        std::vector<std::complex<double>> rec =
            ofdm.filter(data_tx, h_normalized);
        for (std::complex<double>& value : rec) {
          value *= std::sqrt(rho);
        }
        // outputFile << "Outputing received symbols (no noise)" << std::endl;
        // output1DComplexVector(outputFile, rec);

        // Generate Random Noise
        std::vector<std::complex<double>> noise =
            ofdm.generateNoise(0.0, std::sqrt(0.5 / config["L"]), rec.size());

        std::vector<std::complex<double>> received =
            addComplexVectors(rec, noise);

        // remove cyclic prefix
        std::vector<std::vector<std::complex<double>>> rec_sans_cp =
            ofdm.removeCyclicPrefix(received);

        // FFT
        std::vector<std::vector<std::complex<double>>> rec_f =
            ofdm.fft(rec_sans_cp);
        //  outputFile << "Outputing fft result" << std::endl;
        // output2DComplexVector(outputFile, rec_f);

        // channel FFT
        std::vector<std::complex<double>> H_f =
            ofdm.fft(normalize(h), config["L"]);
        for (auto& val : H_f) {
          val *= std::sqrt(rho);
        }

        // Decoding
        std::vector<std::vector<int>> dec_sym = ofdm.decode(rec_f, H_f, mod);
        // outputFile << "Outputing decoding result" << std::endl;
        // for (auto& row : dec_sym) {
        //   output1DVector(outputFile, row);
        // }

        // Compare and count error
        std::vector<int> recsym_flatten = ofdm.columnMajorFlatten(dec_sym);
        // outputFile << "Outputing flattened decoding result" << std::endl;
        // output1DVector(outputFile, recsym_flatten);

        // Bit Error Calculation
        std::vector<int> rec_bits = ofdm.convertIntToBits(recsym_flatten, mod);

        symerr_cnt += ofdm.symbolErrorCount(data_int, recsym_flatten);
        biterr_cnt += ofdm.symbolErrorCount(data_bits, rec_bits);
        // outputFile << "In iteration " << mc_loop << ": , there are " <<
        // ofdm.symbolErrorCount(data_int, recsym_flatten) << " errors" <<
        // std::endl;
      }
      outputFile << "Symbol error after 5000 iterations: " << symerr_cnt
                 << std::endl;
      outputFile << "Bit error after 5000 iterations: " << biterr_cnt
                 << std::endl;
    }
  }
  /*
  // IFFT Test
  std::vector<std::vector<std::complex<double>>> input;
  std::vector<std::complex<double>> row_1 = {};
  std::vector<int> data_tx = {1,  -1, 1,  1,  -1, -1, 1, 1,
                              -1, -1, -1, -1, -1, -1, 1, -1};
  for (int i : data_tx) {
    row_1.push_back(static_cast<double>(i));
  }
  input.push_back(row_1);
  std::vector<std::vector<std::complex<double>>> result = ofdm.ifft(input);
  printComplexVector2D(result);

  // MATLAB DATA API Test
  std::string filename = "../rec.mat";
  std::vector<std::vector<std::complex<double>>> complexData =
      readComplexMatFile(filename, "rec");
  outputFile << "FFT Test" << std::endl;
  outputFile << "Outputing rec data (without noise)" << std::endl;
  output2DComplexVector(outputFile, complexData);

  std::vector<std::vector<std::complex<double>>> rec_sans_cp =
      ofdm.removeCyclicPrefix(complexData[0]);
  outputFile << "Outputing remove cyclic prefix result" << std::endl;
  output2DComplexVector(outputFile, rec_sans_cp);
  std::cout << "Remove cp size: " << rec_sans_cp.size() << " x "
            << rec_sans_cp[0].size() << std::endl;

  // FFT Test
  std::vector<std::vector<std::complex<double>>> rec_f = ofdm.fft(rec_sans_cp);
  outputFile << "Outputing fft result" << std::endl;
  output2DComplexVector(outputFile, rec_f);
  */

  outputFile.close();
  return 0;
}

void print1DVector(const std::vector<int>& vec) {
  for (const int& num : vec) {
    std::cout << num << " ";
  }
  std::cout << std::endl;
}

void output1DVector(std::ofstream& outputFile, const std::vector<int>& vec) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }

  for (int value : vec) {
    outputFile << value << " ";
  }
  outputFile << std::endl;
}

void printComplexVector2D(
    const std::vector<std::vector<std::complex<double>>>& vec2D) {
  for (const auto& row : vec2D) {
    for (const std::complex<double>& value : row) {
      std::cout << value << " ";
    }
    std::cout << std::endl;
  }
}

void output1DComplexVector(std::ofstream& outputFile,
                           const std::vector<std::complex<double>>& vec) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }
  for (const auto& value : vec) {
    outputFile << "(" << value.real() << ", " << value.imag() << ") ";
  }
  outputFile << std::endl;
}

void output2DComplexVector(
    std::ofstream& outputFile,
    const std::vector<std::vector<std::complex<double>>>& vec2D) {
  if (!outputFile.is_open()) {
    std::cerr << "Output file is not open for writing." << std::endl;
    return;
  }

  for (const auto& row : vec2D) {
    for (const std::complex<double>& value : row) {
      outputFile << "(" << value.real() << ", " << value.imag() << ") ";
    }
    outputFile << std::endl;
  }
}

std::vector<std::complex<double>> normalize(
    const std::vector<std::complex<double>>& input) {
  double sum_of_abs_squares = 0.0;
  for (const std::complex<double>& value : input) {
    sum_of_abs_squares += std::norm(value);
  }

  double normalization_factor = std::sqrt(sum_of_abs_squares);

  std::vector<std::complex<double>> normalized_vector;
  for (const std::complex<double>& value : input) {
    normalized_vector.push_back(value / normalization_factor);
  }

  return normalized_vector;
}

std::vector<int> readVectorFromMATFile(const std::string& filePath,
                                       const std::string& variableName) {
  std::vector<int> dataVector;

  // Load the MATLAB file
  MATFile* matFile = matOpen(filePath.c_str(), "r");

  if (matFile == nullptr) {
    std::cerr << "Error: Could not open MATLAB file." << std::endl;
    return dataVector;
  }

  // Read the 1D vector from the MATLAB file
  mxArray* matlabData = matGetVariable(matFile, variableName.c_str());

  if (matlabData == nullptr) {
    std::cerr << "Error: Could not read data from MATLAB file." << std::endl;
    mxDestroyArray(matlabData);
    matClose(matFile);
    return dataVector;
  }

  if (!mxIsDouble(matlabData) || mxGetNumberOfDimensions(matlabData) > 2) {
    std::cerr << "Error: Invalid data format in MATLAB file." << std::endl;
    mxDestroyArray(matlabData);
    matClose(matFile);
    return dataVector;
  }

  double* data = mxGetPr(matlabData);
  int dataSize = static_cast<int>(mxGetNumberOfElements(matlabData));

  // Convert MATLAB data to a C++ vector of ints
  dataVector.reserve(dataSize);
  for (int i = 0; i < dataSize; ++i) {
    dataVector.push_back(static_cast<int>(data[i]) - 1);
  }

  // Clean up
  mxDestroyArray(matlabData);
  matClose(matFile);

  return dataVector;
}

std::vector<std::vector<std::complex<double>>> readComplexMatFile(
    const std::string& filename, const std::string& variableName) {
  // Open the .mat file
  MATFile* matFile = matOpen(filename.c_str(), "r");
  if (!matFile) {
    std::cerr << "Failed to open the .mat file." << std::endl;
    return {};
  }

  // Get the complex variable from the .mat file
  mxArray* mxArrayComplex = matGetVariable(matFile, variableName.c_str());
  if (!mxArrayComplex) {
    std::cerr << "Failed to get the variable from the .mat file." << std::endl;
    matClose(matFile);
    return {};
  }

  // Extract complex data from mxArray and store it in a 2D vector
  mwSize numRows = mxGetM(mxArrayComplex);
  mwSize numCols = mxGetN(mxArrayComplex);

  const std::complex<double>* complexData =
      reinterpret_cast<const std::complex<double>*>(
          mxGetComplexDoubles(mxArrayComplex));

  std::vector<std::vector<std::complex<double>>> complexVector2D(
      numRows, std::vector<std::complex<double>>(numCols));
  for (mwSize i = 0; i < numRows; ++i) {
    for (mwSize j = 0; j < numCols; ++j) {
      complexVector2D[i][j] = complexData[i + j * numRows];
    }
  }

  // Clean up and close the .mat file
  mxDestroyArray(mxArrayComplex);
  matClose(matFile);

  return complexVector2D;
}

std::vector<std::complex<double>> addComplexVectors(
    const std::vector<std::complex<double>>& vector1,
    const std::vector<std::complex<double>>& vector2) {
  // Ensure both input vectors have the same size
  if (vector1.size() != vector2.size()) {
    throw std::invalid_argument("Vector sizes must match for addition.");
  }

  std::vector<std::complex<double>> result;
  result.reserve(vector1.size());

  for (size_t i = 0; i < vector1.size(); ++i) {
    result.push_back(vector1[i] + vector2[i]);
  }

  return result;
}