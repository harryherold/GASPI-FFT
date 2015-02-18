#include <GASPI.h>

#include <iostream>
#include <signal.h>
#include <cstdlib>
#include <cmath>
#include <cstdlib>
#include <assert.h>
#include "utils.hpp"
#include "fft_runtime.hpp"
#include <sstream>
#include <sys/time.h>


#define FFTW_COMPLEX 16
unsigned int cycle = 1;
//--------------------------------------------------------------------------------------------
unsigned long calcMemoryReservation(unsigned long vectorlength, gaspi_rank_t rankcount)
{
    unsigned long PartnerCount = log2( rankcount );

    unsigned long memPerBuffer = ((vectorlength / (unsigned long) rankcount) / 2)
                                 * FFTW_COMPLEX;

    unsigned long resultMem = vectorlength * FFTW_COMPLEX;

    unsigned long exponent = log2( resultMem + memPerBuffer * 2
     + (PartnerCount * (memPerBuffer + sizeof(int))) + sizeof(int)) + 1;

    return (unsigned long) 1 << exponent;
}

//--------------------------------------------------------------------------------------------
bool checkArguments( int argc ,char **argv, bool & val , unsigned long & length )
{
  val = false;
  if(argc < 3)
  {
    std::cout << "Not enough arguments given\n";
    return false;
  } 
  else if( argc == 4 ) 
  {
    if(argv[3][0] == 'v')
    {
      val = true;
    }
    else
    {
      std::cout << "Wrong Mode given\n";
      return false;
    }
  }
  int mem_val = std::atoi(argv[1]);
  unsigned long multiplicator = 0;

  char unit = argv[2][0];

  if(unit == 'G')
  {
    multiplicator = (unsigned long) std::pow(2.0,30);
  } 
  else if(unit == 'M') 
  {
    multiplicator = (unsigned long) std::pow(2.0,20);
  } 
  else if(unit =='U')
  {
    length = mem_val;
    return true;
  } 
  else
  {
    std::cout << "Wrong Unit given\n";
    return false;
  }
  length = (multiplicator * mem_val) / 16;
  return true;
}
//--------------------------------------------------------------------------------------------
int main(int argc, char ** argv)
{
  bool               validation   = false;
  unsigned long      masterLength = 0;

  gaspi_segment_id_t used_segment = 3;
  gaspi_segment_id_t coll_segment = 1;
  gaspi_return_t     ret          = GASPI_SUCCESS;
  gaspi_rank_t       rank;
  gaspi_rank_t       rankcount;

  if ( !checkArguments(argc,argv,validation,masterLength) ) 
  {
    std::cout << "Not enough arguments" << std::endl;
    std::cout << "How to use :" << std::endl;
    std::cout << "./gpi_run.sh -n 16 ./bin/main <size> <memory unit> [ v ]\n";
    std::cout << "memory unit :                        G...Gigabyte\n";
    std::cout << "                                     M...Megabyte\n";
    std::cout << "Options:\n";
    std::cout << "v                    enable the correctness check\n";
    std::cout << "Example one gigabyte with correctness check:\n";
    std::cout << "./gpi_run.sh -n 16 ./bin/main 1 G v\n\n";
    std::cout << "Example eight megabyte without correctness check:\n";
    std::cout << "./gpi_run.sh -n 16 ./bin/main 8 M \n";
    gaspi_proc_term( GASPI_BLOCK );
    return 0;
  }

  if ( gaspi_proc_init( GASPI_BLOCK ) != GASPI_SUCCESS )
  {
    std::cerr << "GASPI is down" << std::endl;
    gaspi_proc_term( GASPI_BLOCK );
    return -1;
  }

  gaspi_proc_num( &rankcount );
  gaspi_proc_rank( &rank );

  struct timeval startTV_incl, endTV_incl, startTV_excl, endTV_excl;

  if( rank == 0 )
  {
    gettimeofday(&startTV_incl, NULL);
  }

  unsigned long initialLength = 0;
  gaspi_pointer_t pRdma;

  ret = gaspi_segment_create( coll_segment,
                              sizeof(unsigned long) * rankcount,
                              GASPI_GROUP_ALL,
                              GASPI_BLOCK,
                              GASPI_MEM_INITIALIZED );

  if( ret != GASPI_SUCCESS )
  {
    gaspi_printf("gaspi_segment_create for coll failed\n");
    return ret;
  }

  gaspi_segment_ptr( coll_segment , &pRdma );

  if(rank == 0 )
  {
    *( (unsigned long *) pRdma) = masterLength;
  }

  gaspi_bcast_binominal( coll_segment, 0UL, sizeof(unsigned long), 0 );

  initialLength = *((unsigned long *) ( pRdma ));

  gaspi_segment_delete( coll_segment );

  gaspi_size_t seg_size = calcMemoryReservation( initialLength, rankcount );

  ret = gaspi_segment_create( used_segment,
                              seg_size,
                              GASPI_GROUP_ALL,
                              GASPI_BLOCK,
                              GASPI_MEM_INITIALIZED );
  if( ret != GASPI_SUCCESS )
  {
    gaspi_printf("gaspi_segment_create for used segment failed\n");
    return ret;
  }
  while( cycle > 0 )
  {
      if( rank == 0 )
        gettimeofday( &startTV_excl, 0 );

      FftRuntime f2(initialLength, 2, used_segment);
      f2.startRuntime();

      gaspi_printf("All done\n");

     if( rank == 0  && validation )
        f2.validateFFT();

      gaspi_barrier( GASPI_GROUP_ALL , GASPI_BLOCK );
      if( rank == 0 )
      {
        gettimeofday( &endTV_excl, 0 );
        gaspi_printf("excl. execution time in secs  : %lu\n",endTV_excl.tv_sec  - startTV_excl.tv_sec);
        gaspi_printf("excl. execution time in usecs : %lu\n",endTV_excl.tv_usec - startTV_excl.tv_usec);
      }
      cycle--;
  }
  if( gaspi_segment_delete( used_segment ) != GASPI_SUCCESS )
  {
    gaspi_printf("Segment-deletion failed\n");
  }
  if( rank == 0 )
  {
    gettimeofday(&endTV_incl, NULL);
    gaspi_printf("incl Execution time in secs  : %lu\n",endTV_incl.tv_sec - startTV_incl.tv_sec);
    gaspi_printf("incl execution time in usecs : %lu\n",endTV_incl.tv_usec - startTV_incl.tv_usec);
  }
  gaspi_proc_term( GASPI_BLOCK );
  return 0;
}
