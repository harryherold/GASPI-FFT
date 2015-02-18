/*
 * local_rdma.hpp
 *
 *  Created on: 15.06.2012
 *      Author: Christian Herold
 */

#ifndef LOCAL_RDMA_HPP_
#define LOCAL_RDMA_HPP_
#include <complex.h>
#include <fftw3.h>



#include <vector>
#include "utils.hpp"

class RdmaManager {

public:
  struct RemoteNodeEntry
  {
    unsigned int    nodeid;
    unsigned long   recvBuffer_Offset;
  };
  std::vector<RemoteNodeEntry> allNodes;
  use_calcBuffer_t             calcBuffer;
  send_t                       sendbuffer;

  static RdmaManager* getInstance( void );
  void                initial( gaspi_segment_id_t seg );
  void                destroyInstance();
  void                checkDmaQueue(gaspi_queue_id_t queue);

  void                setCalcOffsets(unsigned long calcoffset_1, unsigned long calcoffset_2);
  void                setRecvBuffersOffset(unsigned long offset);
  void                setLengthperBuffer(unsigned long length);
  void                setCalcBuffer(use_calcBuffer_t buffer);
  void                setNodeCount(unsigned int nodeCount);
  void                setInitialOffsets(unsigned long offset1, unsigned long offset2);

  void*               getRdmaPointer();
  fftw_complex*       getStartAddress();
  unsigned long       getBufferLength();
  unsigned long       getCalcBufferOffset1();
  unsigned long       getCalcBufferOffset2();
  unsigned long       getInitialOffset1();
  unsigned long       getInitialOffset2();
  fftw_complex        getVectorElement(size_t idx, unsigned int level);

  unsigned long       getRecvBuffersOffset( void );

  void                waitOnNotifies( gaspi_notification_id_t   id_begin,
                                      gaspi_notification_id_t   id_count );
  void                distributeVectors(int splitCount, unsigned long totalVectorLength);
  bool                writeVectorToNode(int level);
  bool                writeResultToMaster(int reverseBitOrderOfRank, unsigned long totalVectorLength);
  void                copyCalcBufferToResultBuffer(unsigned long totalVectorLength);
  void                printNodeEntries();
  void                initialNodeEntries(int * nodes, unsigned int nodeCount);
  double              generateFakeData(size_t idx , unsigned long totalVectorLength);

  fftw_complex *      operator [](size_t idx);

private:
  static const unsigned int intMax = 1073741824;
  gaspi_timeout_t     timeout;

  unsigned long       recvBuffersOffset;
  unsigned long       notifyOffset;
  unsigned long       calcOffset_1;
  unsigned long       calcOffset_2;
  unsigned long       initialOffset_1;
  unsigned long       initialOffset_2;
  unsigned long       bufferlength;
  static RdmaManager* singleton;

  void*                pRdmaSegment;
  gaspi_segment_id_t   used_segment;
  unsigned int         nodecount;
  gaspi_rank_t         rank;
  gaspi_notification_t flag_value;

  fftw_complex          getLocalElement(size_t idx);
  fftw_complex          getRemoteElement(size_t idx, unsigned int level);

  RdmaManager(){}
  ~RdmaManager(){}

};
#endif /* LOCAL_RDMA_HPP_ */
