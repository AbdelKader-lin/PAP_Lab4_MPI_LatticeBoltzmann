/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement non-blocking 1D communication scheme
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
// NEW:
//     - >>> Non-blocking communications <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex3(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation then ex1
	lbm_comm_init_ex1(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex3( lbm_comm_t * comm , lbm_mesh_t * mesh ) {
	
	// Get info
	int rank ;
	int comm_size ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank ) ;
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size ) ;


	// Number of doubles in one col
	int col_size = comm->height * DIRECTIONS ;

	double * left_ghost = lbm_mesh_get_cell( mesh , 0 , 0 ) ;
	double * left_real  = lbm_mesh_get_cell( mesh , 1 , 0 ) ;
	double * right_real = lbm_mesh_get_cell( mesh , comm->width - 2 , 0 ) ;
	double * right_ghost= lbm_mesh_get_cell( mesh , comm->width - 1 , 0 ) ;

	MPI_Request RequArray[ 4 ] ;
	int i = 0 ;

	// left -> right : receive from left then send to right.
	if ( rank > 0 ) {
		MPI_Irecv( left_ghost, col_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD , &RequArray[ i++ ] ) ;
	}
	if ( rank < comm_size - 1 ) {
		MPI_Isend( right_real, col_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD , &RequArray[ i++ ] ) ;
	}

	// right -> left : receive from right then send to left.
	if ( rank < comm_size - 1 ) {
		MPI_Irecv( right_ghost, col_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD , &RequArray[ i++ ] ) ;
	}
	if ( rank > 0 ) {
		MPI_Isend( left_real, col_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD , &RequArray[ i++ ] ) ;
	}

	MPI_Waitall( i , RequArray , MPI_STATUSES_IGNORE ) ; // i because it is the number of requests

}
