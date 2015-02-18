/*
 * fft_computation.cpp
 *
 *  Created on: 21.06.2012
 *      Author: Christian Herold
 */
#include <iostream>
#include <math.h>
#include <assert.h>
#include "fft_computation.hpp"

//------------------------------------------------------------------------------
FftComputation::FftComputation(unsigned long length)
{
  this->vectorlength = length;
  rdma = RdmaManager::getInstance();
  twiddles = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) * rdma->getBufferLength() );
  assert(twiddles);

  finalVector = (fftw_complex *) (((char *) rdma->getRdmaPointer())
      + rdma->getCalcBufferOffset1());

  srcVector   = (fftw_complex *) (((char *) rdma->getRdmaPointer())
      + rdma->getInitialOffset1());
}

//------------------------------------------------------------------------------
FftComputation::~FftComputation()
{
  if (twiddles != NULL) {
    fftw_free(twiddles);
    twiddles = NULL;
  }
}

//------------------------------------------------------------------------------
void FftComputation::calculateFftw()
{
  finalVector = (fftw_complex *) (((char *) rdma->getRdmaPointer())
      + rdma->getCalcBufferOffset1());

  srcVector   = (fftw_complex *) (((char *) rdma->getRdmaPointer())
      + rdma->getInitialOffset1());

  fftwPlan = fftw_plan_dft_1d(vectorlength, srcVector, finalVector,
        FFTW_FORWARD, FFTW_ESTIMATE);

  assert(fftwPlan);

  fftw_execute(fftwPlan);

  fftw_destroy_plan(fftwPlan);

}

//------------------------------------------------------------------------------
void FftComputation::printFftw()
{
  gaspi_printf("Result of 1 d fftw\n");
  finalVector = (fftw_complex *) (((char *) rdma->getRdmaPointer())
      + rdma->getCalcBufferOffset1());
  for (unsigned long i = 0; i < vectorlength; i++)
    gaspi_printf("%lf + %lf i\n", creal(finalVector[i]), cimag(finalVector[i]));
}

//------------------------------------------------------------------------------
unsigned long FftComputation::getVectorLength()
{
  return vectorlength;
}

//------------------------------------------------------------------------------
void FftComputation::calculateTwiddles(unsigned long kmin,
                                       unsigned long mergelength)
{
  for (unsigned int i = 0; i < rdma->getBufferLength(); i++, kmin++) {
    twiddles[i] = cexp(I * M_PI * 2 * kmin * -1 / mergelength);
  }
}

//------------------------------------------------------------------------------
void FftComputation::radix2FFT(unsigned int level)
{
  unsigned long bufferlength = rdma->getBufferLength();
  fftw_complex minuend = 0.0, subrathend = 0.0;

  for (unsigned long i = 0; i < bufferlength; i++)
  {

    minuend = rdma->getVectorElement(i, level);

    subrathend = rdma->getVectorElement(i + bufferlength, level);

    (*(*rdma)[i]) = minuend + (subrathend * twiddles[i]);

    (*(*rdma)[i + bufferlength]) = minuend - (subrathend * twiddles[i]);
  }
}
