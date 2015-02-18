/*
 * local_rdma.cpp
 *
 *  Created on: 15.06.2012
 *      Author: Christian Herold
 */
#include <iostream>
#include <GASPI.h>
#include <assert.h>
#include <cstring>
#include <cmath>
#include "rdma_manager.hpp"

RdmaManager * RdmaManager::singleton = NULL;

RdmaManager * RdmaManager::getInstance()
{
  if (singleton == NULL) 
  {
    singleton = new RdmaManager();
    assert(singleton);
  }

  return singleton;
}
//------------------------------------------------------------------------------
void RdmaManager::initial( gaspi_segment_id_t seg )
{
  pRdmaSegment = NULL;
  used_segment = seg;
  gaspi_segment_ptr( used_segment , &pRdmaSegment );
  gaspi_proc_rank( &rank );
  timeout = GASPI_BLOCK;
  flag_value = 42;
}
//------------------------------------------------------------------------------
void RdmaManager::destroyInstance()
{
  if (singleton != NULL) 
  {
    delete singleton;
    singleton = NULL;
  }
}

//------------------------------------------------------------------------------
void * RdmaManager::getRdmaPointer()
{
  gaspi_segment_ptr( used_segment , &pRdmaSegment );
  return pRdmaSegment;
}

//------------------------------------------------------------------------------
fftw_complex * RdmaManager::getStartAddress()
{
  gaspi_pointer_t pSegment = NULL;
  gaspi_segment_ptr( used_segment , &pSegment );
  return ((fftw_complex *) pRdmaSegment);
}

//------------------------------------------------------------------------------
void RdmaManager::setCalcOffsets(unsigned long calcoffset_1,
    unsigned long calcoffset_2)
{
  calcOffset_1 = calcoffset_1;
  calcOffset_2 = calcoffset_2;
}

//------------------------------------------------------------------------------
fftw_complex * RdmaManager::operator [](size_t idx)
{
  char *ptr = (char *) pRdmaSegment;
  if ((idx < bufferlength) && (idx >= 0)) 
  {
    ptr += calcOffset_1 + (sizeof(fftw_complex) * idx);
    fftw_complex * tmp = (fftw_complex *) ptr;
    return tmp;
  }
  else if (((idx >= bufferlength) && (idx < (2 * bufferlength))))
  {
    ptr += calcOffset_2 + (sizeof(fftw_complex) * (idx - bufferlength));
    fftw_complex * tmp = (fftw_complex *) ptr;
    return tmp;
  }
  else
  {
    std::cerr << "ERROR # RdmaManager::operator [] # index out of bound"
        << std::endl;
    return NULL;
  }
}
//------------------------------------------------------------------------------
unsigned long RdmaManager::getRecvBuffersOffset( void )
{
  return recvBuffersOffset;
}
//------------------------------------------------------------------------------
void RdmaManager::setRecvBuffersOffset(unsigned long offset)
{
  recvBuffersOffset = offset;
}

//------------------------------------------------------------------------------
fftw_complex RdmaManager::getLocalElement(size_t idx)
{
  char *ptr = (char *) pRdmaSegment;
  if (idx > bufferlength)
  {
    std::cerr << "ERROR # RdmaManager::getLocalElement # Index is to high"
        << std::endl;
    return 0.0;
  }
  else if (idx < 0)
  {
    std::cerr << "ERROR # RdmaManager::getLocalElement # Index lesser than 0"
        << std::endl;
    return 0.0;
  }
  
  if (sendbuffer == calc_buffer1)
    ptr += calcOffset_2;
  else if (sendbuffer == calc_buffer2)
    ptr += calcOffset_1;
  else
  {
    std::cerr << "unexpected value of sendbuffer in getLocalElement"
        << std::endl;
    return 0.0;
  }
  return (((fftw_complex *) ptr)[idx]);
}

//------------------------------------------------------------------------------
fftw_complex RdmaManager::getRemoteElement(size_t idx, unsigned int level)
{
  char *ptr = (char *) pRdmaSegment;
  if (idx > bufferlength) 
  {
    std::cerr << "ERROR # RdmaManager::getRemoteElement # Index is to high"
        << std::endl;
    return 0.0;
  }
  else if (idx < 0)
  {
    std::cerr << "ERROR # RdmaManager::getRemoteElement # Index lesser than 0"
        << std::endl;
    return 0.0;
  }
  ptr += allNodes[level].recvBuffer_Offset;

  return (((fftw_complex *) ptr)[idx]);
}

//------------------------------------------------------------------------------
fftw_complex RdmaManager::getVectorElement(size_t idx, unsigned int level)
{
  if ((idx < bufferlength) && (sendbuffer == calc_buffer2))
    return getLocalElement(idx);
  else if ((idx >= bufferlength) && (sendbuffer != calc_buffer2))
    return getLocalElement(idx - bufferlength);
  else if ((idx >= bufferlength) && (sendbuffer == calc_buffer2))
    return getRemoteElement(idx - bufferlength, level);
  else if (((idx < bufferlength) && (sendbuffer != calc_buffer2)))
    return getRemoteElement(idx, level);
  else
    std::cerr
        << "unexpected combination of index and sendbuffer in getVectorElement() "
        << std::endl;
  return 0.0;
}

//------------------------------------------------------------------------------
void RdmaManager::checkDmaQueue( gaspi_queue_id_t queue )
{
  gaspi_number_t queue_size = 0;
  gaspi_queue_size( queue , &queue_size);

  gaspi_number_t queue_size_max = 0;
  gaspi_queue_size_max ( &queue_size_max );

  if(  queue_size >= queue_size_max )
  {
    if( GASPI_ERROR == gaspi_wait( queue , GASPI_BLOCK ) )
    {
      gaspi_printf("wait failed\n");
    }
  }
}

//------------------------------------------------------------------------------
/*
 * @Override
 * needed levelcount and nodelist
 * i.e.
 * nodelist[2] = 2^2
 * nodelist[1] = 2^1
 * nodelist[0] = 2^0
 */
/*
 * Initialisierung allNodes[node]->recvBuffer_Offset ändern
 * Jede Node schreibt den passenden Offset zu seinem Recv-Buffer
 * in die Nachbarhosts , damit hat jede Node nur einen Offsetwert
 * von jeder Node
 */

//------------------------------------------------------------------------------
bool RdmaManager::writeVectorToNode(int level)
{
  unsigned int sendOffset = 0;
  gaspi_return_t ret;

  if (sendbuffer == calc_buffer2)
  {
    sendOffset = calcOffset_2;
  }
  else if (sendbuffer == calc_buffer1)
  {
    sendOffset = calcOffset_1;
  }

  unsigned long remoteOffset = allNodes[level].recvBuffer_Offset;
  gaspi_queue_id_t queue0  = 0;
  unsigned int nodeid = allNodes[level].nodeid;


  unsigned long send_size =  bufferlength * sizeof(fftw_complex);
  int maxSends = 1;

  if (send_size > intMax)
  {
    maxSends = send_size / intMax;
    send_size = intMax;
  }
  for(int i = 0; i < maxSends ; i++)
  {
    checkDmaQueue(queue0);
    ret = gaspi_write( used_segment,
                       sendOffset + (i * send_size),
                       nodeid,
                       used_segment,
                       remoteOffset + (i * send_size),
                       send_size,
                       queue0,
                       GASPI_BLOCK);

    if (ret != GASPI_SUCCESS)
    {
      std::cerr << "ERROR # writeVectorToNode() # write Dma failed" << std::endl;
      exit(2);
    }
  }
  ret = gaspi_notify( used_segment,
		  	  	  	  nodeid,
					  level,
					  42,
					  queue0,
					  GASPI_BLOCK );
  gaspi_printf("write notify %d to %d\n", level, nodeid);
  if (ret != GASPI_SUCCESS)
  {
	  std::cerr << "ERROR # writeVectorToNode() # write Dma failed" << std::endl;
      exit(2);
  }

  return true;
}

//------------------------------------------------------------------------------
void RdmaManager::copyCalcBufferToResultBuffer(unsigned long totalVectorLength)
{

  fftw_complex * pEvenSrc = (fftw_complex *) ((char *) pRdmaSegment
      + calcOffset_1);

  fftw_complex * pEvenDest = (fftw_complex *) ((char *) pRdmaSegment
      + initialOffset_1);

  memcpy(pEvenDest, pEvenSrc, bufferlength * sizeof(fftw_complex));

  fftw_complex * pOddSrc = (fftw_complex *) ((char *) pRdmaSegment
      + calcOffset_2);

  fftw_complex * finalVector = (fftw_complex *) (((char *) pRdmaSegment)
      + getCalcBufferOffset1());

  unsigned long oddEvenDispl = (totalVectorLength / 2) * sizeof(fftw_complex);



  fftw_complex * pOddDest = (fftw_complex *) ((char *) pRdmaSegment
      + initialOffset_1 + oddEvenDispl);



  memcpy(pOddDest, pOddSrc, bufferlength * sizeof(fftw_complex));

}

//------------------------------------------------------------------------------
bool RdmaManager::writeResultToMaster(  int reverseBitOrderOfRank,
                                        unsigned long totalVectorLength )
{

  unsigned long evenOffset = initialOffset_1
      + (reverseBitOrderOfRank * bufferlength * sizeof(fftw_complex));


  unsigned long oddOffset = evenOffset
      + ((totalVectorLength / 2) * sizeof(fftw_complex));

  unsigned long send_size = bufferlength * sizeof(fftw_complex);

  int maxSends = 1;

  if (send_size > intMax)
  {
    maxSends = send_size / intMax;
    send_size = intMax;
  }

  for(int i = 0; i < maxSends ; i++)
  {
    checkDmaQueue(0);
    gaspi_return_t ret = gaspi_write( used_segment,
                                      calcOffset_1 + (i * send_size),
                                      0,
                                      used_segment,
                                      evenOffset + (i * send_size),
                                      send_size,
                                      0,
                                      GASPI_BLOCK);

    if (ret != GASPI_SUCCESS ) {
      std::cerr << "ERROR # writeResultToMaster() # write Dma failed"
          << std::endl;
      return false;
    }

    checkDmaQueue(0);
    ret = gaspi_write_notify( used_segment,
                              calcOffset_2 + (i * send_size),
                              0,
                              used_segment,
                              oddOffset + (i * send_size),
                              send_size,
                              ((int) log2(nodecount)) + rank,
                              42,
                              0,
                              GASPI_BLOCK);
    if (ret != GASPI_SUCCESS) {
      std::cerr << "ERROR # writeResultToMaster() # write Dma failed"
          << std::endl;
      return false;
    }

  }
  return true;
}

//------------------------------------------------------------------------------
void RdmaManager::initialNodeEntries(int * nodes, unsigned int nodeCount)
{
  for (unsigned int i = 0; i < nodeCount; i++) {
    allNodes.push_back(RemoteNodeEntry());
    allNodes[i].nodeid = nodes[i];
    allNodes[i].recvBuffer_Offset = recvBuffersOffset
        + ( bufferlength * sizeof(fftw_complex) ) * i;
  }
}

//------------------------------------------------------------------------------
void RdmaManager::printNodeEntries()
{
  for (unsigned long i = 0; i < allNodes.size(); i++) {
    gaspi_printf("Level %ld id %d bufferoffset1: %ld\n", i, allNodes[i].nodeid,
        allNodes[i].recvBuffer_Offset);
  }
}

//------------------------------------------------------------------------------
unsigned long RdmaManager::getBufferLength()
{
  return bufferlength;
}

//------------------------------------------------------------------------------
unsigned long RdmaManager::getCalcBufferOffset1()
{
  return calcOffset_1;
}

//------------------------------------------------------------------------------
unsigned long RdmaManager::getCalcBufferOffset2()
{
  return calcOffset_2;
}

//------------------------------------------------------------------------------
unsigned long RdmaManager::getInitialOffset1()
{
  return initialOffset_1;
}

//------------------------------------------------------------------------------
unsigned long RdmaManager::getInitialOffset2()
{
  return initialOffset_2;
}

//------------------------------------------------------------------------------
void RdmaManager::setLengthperBuffer(unsigned long length)
{
  bufferlength = length;
}

//------------------------------------------------------------------------------
void RdmaManager::setNodeCount(unsigned int nodeCount)
{
  nodecount = nodeCount;
}

//------------------------------------------------------------------------------
void RdmaManager::setInitialOffsets(unsigned long offset1,
    unsigned long offset2)
{
  initialOffset_1 = offset1;
  initialOffset_2 = offset2;
}
//------------------------------------------------------------------------------
double RdmaManager::generateFakeData(size_t idx , unsigned long totalVectorLength)
{
  unsigned long sig = totalVectorLength / 8;
  unsigned long i = 1;
  for (; i <= 8; i++) {
    if (idx < (i * sig)) {
      if (i % 2 != 0)
        return 1.0;
      else
        return -1.0;
    }
  }
  std::cerr << "Initial-Error" << std::endl;

  return 0.0;
}
//------------------------------------------------------------------------------
void RdmaManager::waitOnNotifies( gaspi_notification_id_t   id_begin,
                                  gaspi_notification_id_t   id_count )
{
  gaspi_notification_t    tmp;
  gaspi_return_t          retval;
  gaspi_notification_id_t first_id;

  for( int i = id_begin ; i < ( id_begin + id_count ) ; i++ )
  {
    retval = gaspi_notify_waitsome( used_segment,
                                    id_begin,
                                    id_count,
                                    &first_id,
                                    GASPI_BLOCK );

    if( retval == GASPI_ERROR )
    {
      gaspi_printf("Wait-Error in waitOnInitialVector\n");
    }
    gaspi_notify_reset( used_segment, first_id , &tmp );
  }
}
//------------------------------------------------------------------------------
void RdmaManager::distributeVectors(int splitCount, unsigned long totalVectorLength)
{
  gaspi_segment_id_t GaspiQueue = 0;
  unsigned long initial_offsets[2] = { initialOffset_1 , initialOffset_2 };
  unsigned long send_size = bufferlength * splitCount * sizeof(fftw_complex);
  int maxSends = 1;
  gaspi_return_t retval;

  if (send_size > intMax)
  {
    maxSends = (send_size / intMax);
    send_size = intMax;
  }

  for (int node = 1; node < nodecount; node++)
  {
    GaspiQueue = node % 2;

    fftw_complex * pInitialBuffer = (fftw_complex *) ((char *) pRdmaSegment
      + initial_offsets[GaspiQueue]);

    for (unsigned long i = 0; i < (bufferlength * splitCount); i++)
    {
      pInitialBuffer[i] = generateFakeData(node + (i * nodecount),totalVectorLength);
    }

    for (int i = 0; i < maxSends; i++)
    {
      checkDmaQueue(GaspiQueue);
      retval = gaspi_write(  used_segment,
                             initial_offsets[node % 2] + (i * send_size),
                             node,
                             used_segment,
                             initialOffset_1 + (i * send_size),
                             send_size,
                             GaspiQueue,
                             GASPI_BLOCK );

      if (retval == GASPI_ERROR)
      {
        std::cerr << "write_notify failed in function distributeVectors()"
        << std::endl;
        return;
      }
    }

    retval = gaspi_notify( used_segment, node, rank + 10, flag_value, GaspiQueue, GASPI_BLOCK);

    if (retval == GASPI_ERROR)
    {
    	std::cerr << "write_notify failed in function distributeVectors()"
        << std::endl;
        return;
    }
  }

  fftw_complex * pInitialBuffer_1 = (fftw_complex *) ((char *) pRdmaSegment
      + initialOffset_1);

  gaspi_wait( 0 , GASPI_BLOCK );

  for (unsigned long i = 0; i < (bufferlength * splitCount); i++)
  {
    pInitialBuffer_1[i] = generateFakeData((i * nodecount),totalVectorLength);
  }
}
