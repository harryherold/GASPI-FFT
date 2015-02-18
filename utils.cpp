#include "utils.hpp"

unsigned int
npot( unsigned int v )
{
  v--;
  v |= v >>  1;
  v |= v >>  2;
  v |= v >>  4;
  v |= v >>  8;
  v |= v >> 16;
  return v + 1;
}

int
calculate_comm_partners( int*  parent,
                         int** children,
                         int   me,
                         int   root )
{
  gaspi_rank_t size;
  gaspi_proc_num( &size );
  unsigned int size_pot = npot( size );

  unsigned int d;
  unsigned int number_of_children = 0;

  /* Be your own parent, ie. the root, by default */
  *parent = me;

  me -= root;
  if( me < 0 )
    me += size;

  /* Calculate the number of children for me */
  for ( d = 1; d ; d <<= 1 ) {
    /* Actually break condition */
    if ( d > size_pot ) {
      break;
    }

    /* Check if we are actually a child of someone */
    if ( me & d ) {

      /* Yes, set the parent to our real one, and stop */
      *parent = me ^ d;

      (*parent) += root;
      if( (*parent) >= size )
        (*parent) -= size;
      break;
    }

    /* Only count real children, of the virtual hypercube */
    if ( ( me ^ d ) < size ) {
      number_of_children++;
    }
  }

  /* Put the ranks of all children into a list and return */
  *children = (int *) malloc( sizeof( **children ) * number_of_children );
  unsigned int child = number_of_children;

  d >>= 1;
  while ( d ) {
    if ( ( me ^ d ) < size ) {
      ( *children )[ --child ] = me ^ d;
      ( *children )[ child ] += root;
      if( ( *children )[ child ] >= size )
        ( *children )[ child ] -= size;
    }
    d >>= 1;
  }

  return number_of_children;
}
gaspi_rank_t
gaspi_bcast_binominal(  gaspi_segment_id_t  seg_id,
                        unsigned long       offset,
                        unsigned long       bytesize,
                        gaspi_rank_t        root )
{
  int                     children_count;
  int                     child;
  int                     parent;
  int*                    children = NULL;

  gaspi_notification_id_t first_id;
  gaspi_return_t          retval;
  gaspi_rank_t            rank;
  gaspi_timeout_t         timeout = 2000;
  gaspi_notification_id_t notify_id = 0;
  gaspi_rank_t            rankcount;
  short                   queue = 0;

  gaspi_proc_rank( &rank );
  gaspi_proc_num( &rankcount );

  children_count = calculate_comm_partners( &parent , &children , rank , root );

  /*
   * parents + children wait for upper parents data
   */
  if( rank != parent )
  {
    retval = gaspi_notify_waitsome( seg_id,
                                    notify_id ,
                                    1 ,
                                    (gaspi_notification_id_t * const) &first_id,
                                    timeout );

    gaspi_notification_t val = 0;
    gaspi_notify_reset( seg_id, first_id , &val );
  }
  /*
   * write to all childs
   */
  for ( child = 0; child < children_count ; child++ )
  {
    retval = gaspi_write_notify(  seg_id,
                                  offset,
                                  children[child],
                                  seg_id,
                                  offset,
                                  bytesize,
                                  notify_id,
                                  43,
                                  queue,
                                  GASPI_BLOCK );
    if( retval != 0 )
    {
      std::cerr << "write data to %d failed " << children[child] << "\n";
    }
    gaspi_wait( queue , GASPI_BLOCK);
  }
  gaspi_barrier( GASPI_GROUP_ALL , GASPI_BLOCK );
  free( children );
  return retval;
}
