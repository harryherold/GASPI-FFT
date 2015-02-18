/*
 * fft_runtime.cpp
 *
 *  Created on: 21.06.2012
 *      Author: Christian Herold
 */
#include <iostream>

#include <GASPI.h>
#include <unistd.h>

#include <complex.h>
#include <cmath>
#include <assert.h>
#include "fft_runtime.hpp"

FftRuntime::FftRuntime(unsigned long vectorlength,
					   unsigned int splitcount,
					   gaspi_segment_id_t seg )
:master_rank( 0 )
{
  gaspi_proc_rank( &rank );
  gaspi_proc_num( &nodecount );
  splitCount = splitcount;
  /*
   * get the number of communication partner
   * and store these in the rdma (sorted by iteration)
   */
  levelCount = log2(nodecount);
  /*
   * Initial the portion of the RDMA per Node
   */
  totalVectorLength = vectorlength;

  rdma = RdmaManager::getInstance();
  rdma->initial( seg );
  initialOffsets();

  nodes = new int[levelCount];
  assert(nodes);

  for (int i = levelCount - 1; i >= 0; i--) {
    nodes[i] = getActualMergeNodeID(i);
  }

  rdma->initialNodeEntries(nodes, levelCount);
  /*
   * send data to the worker nodes
   */
  compute = new FftComputation(rdma->getBufferLength() * splitCount);
  assert(compute);
  gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK);
  if ( rank == master_rank )
  {
    rdma->distributeVectors(splitCount, totalVectorLength);
  }
  else
  {
    rdma->waitOnNotifies( master_rank + 10 , 1 );
  }
}
//------------------------------------------------------------------------------
void FftRuntime::initialOffsets()
{
  unsigned long bufferLengthPerNode = (totalVectorLength
      / (nodecount * splitCount));

  unsigned long buffersize = bufferLengthPerNode * sizeof(fftw_complex);

  rdma->setLengthperBuffer(bufferLengthPerNode);

  rdma->setCalcOffsets((unsigned long) 0, buffersize);

  rdma->setRecvBuffersOffset( 2 * buffersize );

  rdma->setNodeCount(nodecount);

  unsigned int recvBufferCount = levelCount;

  unsigned long initalSize = buffersize * splitCount;

  unsigned long initialOffset1 = (2 * buffersize) + recvBufferCount * buffersize;

  unsigned long initialOffset2 = initialOffset1 + initalSize;

  rdma->setInitialOffsets(initialOffset1, initialOffset2);

}

//------------------------------------------------------------------------------
int FftRuntime::getActualMergeNodeID(int expOf2)
{
  int neighbour = -1;
  int border = (int) pow(2.0, expOf2);
  for (int i = border; i < nodecount; i += ((int) pow(2.0, expOf2 + 1))) {
    if (rank < i) {
      neighbour = rank + border;
      break;
    } else if (rank < (i + border)) {
      neighbour = rank - border;
      break;
    }
  }
  if (neighbour == -1) {
    std::cerr << "ERROR # getActualMergeNodeID #Can't find a merge neighbor"
        << std::endl;
    return -1;
  }
  return neighbour;
}

//------------------------------------------------------------------------------
int FftRuntime::calcReverseBitOrder(int number)
{
  int num = number, i, h;
  int digits = levelCount;
  for (h = i = 0; i < digits; i++) {
    h = h << 1;
    h += (num & 1);
    num >>= 1;
  }
  return h;
}

//------------------------------------------------------------------------------
unsigned long FftRuntime::getStartPosInGroup(int exponent)
{
  int groupsize = nodecount / pow(2.0, exponent);
  int pos = calcReverseBitOrder(rank);
  unsigned long kmin = pos - (((int) (pos / groupsize)) * groupsize);
  return kmin * rdma->getBufferLength();
}

//------------------------------------------------------------------------------
void FftRuntime::startRuntime()
{

  int levelcounter = levelCount - 1;

  compute->calculateFftw();
  //----------------------------------------------------------------------------

  while (0 <= levelcounter)
  {

    int mergeNode = getActualMergeNodeID(levelcounter);

    if (rank < mergeNode) {
      rdma->sendbuffer = calc_buffer2;
      rdma->writeVectorToNode(levelcounter);
      unsigned long kmin = getStartPosInGroup(levelcounter);
      unsigned long mergeLength = totalVectorLength / pow(2.0, levelcounter);
      compute->calculateTwiddles(kmin, mergeLength);
    } else {
      rdma->sendbuffer = calc_buffer1;
      rdma->writeVectorToNode(levelcounter);
      unsigned long kmin = getStartPosInGroup(levelcounter);
      unsigned long mergeLength = totalVectorLength / pow(2.0, levelcounter);
      compute->calculateTwiddles(kmin, mergeLength);
    }

    gaspi_printf("Wait on Notify %d\n",levelcounter);
    rdma->waitOnNotifies( levelcounter , 1 );
    gaspi_wait( 0 , GASPI_BLOCK );
    compute->radix2FFT(levelcounter);

    levelcounter--;
  }

  gaspi_printf("Main Computation finished\n");

  if (rank == 0)
  {
    rdma->copyCalcBufferToResultBuffer(totalVectorLength);
    rdma->waitOnNotifies( ((int) log2(nodecount)) + 1 , nodecount - 1 );
  }
  else
  {
    rdma->writeResultToMaster(calcReverseBitOrder(rank), totalVectorLength);
    gaspi_wait( 0 , GASPI_BLOCK );
  }
}

//------------------------------------------------------------------------------
void FftRuntime::validateFFT()
{
  fftw_complex * pResult = (fftw_complex *) ((char *) rdma->getRdmaPointer()
      + rdma->getInitialOffset1());

  fftw_complex * in = (fftw_complex *) fftw_malloc(
      sizeof(fftw_complex) * totalVectorLength);
  fftw_complex * out = (fftw_complex *) fftw_malloc(
      sizeof(fftw_complex) * totalVectorLength);
  fftw_plan plan = fftw_plan_dft_1d(totalVectorLength, in, out, FFTW_FORWARD,
      FFTW_ESTIMATE);
  for (unsigned long i = 0; i < totalVectorLength; i++) {
    in[i] = rdma->generateFakeData(i,totalVectorLength);
  }
  fftw_execute(plan);

  gaspi_printf("Result of 1d FFT\n");
  double diff = 0.0;
  double max = 0.0;
  for (unsigned long i = 0; i < totalVectorLength; i++) {
    double tmp_diff_real = std::abs(creal(out[i]))
        - std::abs(creal(pResult[i]));
    double tmp_max_real =
        std::abs(creal(out[i])) > std::abs(creal(pResult[i])) ?
            creal(out[i]) : creal(pResult[i]);

    diff = std::abs(tmp_diff_real) > std::abs(diff) ? tmp_diff_real : diff;
    max = std::abs(tmp_max_real) > std::abs(max) ? tmp_max_real : max;

    double tmp_diff_imag = std::abs(cimag(out[i]))
        - std::abs(cimag(pResult[i]));
    double tmp_max_imag =
        std::abs(cimag(out[i])) > std::abs(cimag(pResult[i])) ?
            cimag(out[i]) : cimag(pResult[i]);

    diff = std::abs(tmp_diff_imag) > std::abs(diff) ? tmp_diff_imag : diff;
    max = std::abs(tmp_max_imag) > std::abs(max) ? tmp_max_imag : max;

  }
  std::cout << "Relativer Fehler " <<  (double) (std::abs(diff) / std::abs(max)) << "\n";
  std::cout << "Abweichnung max. " << std::abs(diff) << "\n";
  std::cout << "max factor " << std::abs(max) << "\n";

  fftw_free(in);
  fftw_free(out);
  fftw_destroy_plan(plan);
}

//------------------------------------------------------------------------------
FftRuntime::~FftRuntime()
{
  rdma->destroyInstance();
  delete compute;
  delete nodes;
}
