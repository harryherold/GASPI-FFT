/*
 * utils.hpp
 *
 *  Created on: 15.06.2012
 *      Author: Christian Herold
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <GASPI.h>

#include <iostream>
#include <stdlib.h>

typedef enum use_CalcBuffer_t {
  calc_buffer_1, calc_buffer_2
} use_calcBuffer_t;

//------------------------------------------------------------------------------

typedef enum use_RecvBuffer_t {
  recv_buffer_1, recv_buffer_2
} use_recvBuffer_t;

//------------------------------------------------------------------------------

typedef enum Send_t {
  calc_buffer1, calc_buffer2
} send_t;

gaspi_rank_t
gaspi_bcast_binominal(  gaspi_segment_id_t  seg_id,
                        unsigned long       offset,
                        unsigned long       bytesize,
                        gaspi_rank_t        root );

#endif /* UTILS_HPP_ */
