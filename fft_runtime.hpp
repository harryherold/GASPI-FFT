/*
 * fft_runtime.hpp
 *
 *  Created on: 21.06.2012
 *      Author: Christian Herold
 */

#ifndef FFT_RUNTIME_HPP_
#define FFT_RUNTIME_HPP_
#include <complex.h>
#include <fftw3.h>
#include "utils.hpp"
#include "fft_computation.hpp"
#include "rdma_manager.hpp"

class FftRuntime {

public:

  explicit FftRuntime(unsigned long vectorlength, unsigned int splitCount, gaspi_segment_id_t seg);
  ~FftRuntime();
  void startRuntime();
  void distributeVectors();
  void receiveVector();
  void initialOffsets();
  void validateFFT();
  int getActualMergeNodeID(int expOf2);
  int calcReverseBitOrder(int number);
  unsigned long getStartPosInGroup(int exponent);
  double generateFakeData(size_t idx);

private:
  FftComputation *compute;
  RdmaManager * rdma;
  int * nodes;
  gaspi_rank_t  rank;
  int levelCount;
  gaspi_rank_t nodecount;
  gaspi_rank_t master_rank;
  unsigned int splitCount;
  unsigned long totalVectorLength;
  static const unsigned int intMax = 1073741824;

};

#endif /* FFT_RUNTIME_HPP_ */
