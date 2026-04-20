/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement odd/even 1D blocking communication scheme 
//       along X axis.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
// NEW:
//     - >>> Odd/even communication ordering <<<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex2( lbm_comm_t * comm , int total_width , int total_height ) {
	//we use the same implementation then ex1
	lbm_comm_init_ex1( comm , total_width , total_height ) ;
}

/****************************************************/
void lbm_comm_ghost_exchange_ex2( lbm_comm_t * comm , lbm_mesh_t * mesh ) {
	//
	// TODO: Implement the 1D communication with blocking MPI functions using
	//       odd/even communications.
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	
	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

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


	if ( rank % 2 != 0 ){ // First send, then receive

		if ( rank < comm_size - 1 ) {
			MPI_Send( right_real, col_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD ) ;
		}
		if ( rank > 0 ) {
			MPI_Recv( left_ghost, col_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
		}

		if ( rank > 0 ) {
			MPI_Send( left_real, col_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD ) ;
		}
		if ( rank < comm_size - 1 ) {
			MPI_Recv( right_ghost, col_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
		}

	} else { // R then S

		if ( rank > 0 ) {
			MPI_Recv( left_ghost, col_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
		}
		if ( rank < comm_size - 1 ) {
			MPI_Send( right_real, col_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD ) ;
		}

		if ( rank < comm_size - 1 ) {
			MPI_Recv( right_ghost, col_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
		}
		if ( rank > 0 ) {
			MPI_Send( left_real, col_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD ) ;
		}

	}

}
