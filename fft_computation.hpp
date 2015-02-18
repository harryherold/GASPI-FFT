/*
 * fft_computation.hpp
 *
 *  Created on: 21.06.2012
 *      Author: Christian Herold
 */
#ifndef FFT_COMPUTATION_HPP
#define FFT_COMPUTATION_HPP
#include <complex.h>
#include <fftw3.h>
#include "rdma_manager.hpp"

class FftComputation {

public:
  explicit FftComputation(unsigned long length);
  ~FftComputation();

  void radix2FFT(unsigned int level);

  void calculateTwiddles(unsigned long kmin, unsigned long mergelength);

  void calculateFftw();

  void printFftw();

  unsigned long getVectorLength();

private:
  fftw_complex * srcVector;
  fftw_complex * finalVector;
  unsigned long vectorlength;
  unsigned long totalLength;
  fftw_plan fftwPlan;
  RdmaManager * rdma;
  fftw_complex * twiddles;

};
#endif
