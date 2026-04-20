/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
//
// GOAL: Implement a 1D communication scheme along
//       X axis with blocking communications.
//
// SUMMARY:
//     - 1D splitting along X
//     - Blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/****************************************************/
void lbm_comm_init_ex1( lbm_comm_t * comm , int total_width , int total_height ) {
	//
	// TODO: calculate the splitting parameters for the current task.
	//
	// HINT: You can look in exercise_0.c to get an example for the sequential case.
	//

	// Get info
	int rank ;
	int comm_size ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank ) ;
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size ) ;

	// TODO: calculate the number of tasks along X axis and Y axis.
	// 1D Splitting : We only split along the X axis
	comm->nb_x = comm_size ;
	comm->nb_y = 1 ;

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank ;
	comm->rank_y = 0 ; // Because we only split along X, so all processes are position 0.

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = ( total_width ) / ( comm_size ) + 2 ; // '+2' for the ghost cells
	comm->height = total_height + 2 ;

	// TODO : calculate the absolute position in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	comm->x = rank * ( total_width ) / ( comm_size ) ;
	comm->y = 0 ;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex1 ( lbm_comm_t * comm , lbm_mesh_t * mesh ) {
	//
	// TODO: Implement the 1D communication with blocking MPI functions (MPI_Send & MPI_Recv)
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[DIRECTIONS] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
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

	// left -> right : receive from left then send to right.
	if ( rank > 0 ) {
		MPI_Recv( left_ghost, col_size, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
	}
	if ( rank < comm_size - 1 ) {
		MPI_Send( right_real, col_size, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD ) ;
	}

	// right -> left : receive from right then send to left.
	if ( rank < comm_size - 1 ) {
		MPI_Recv( right_ghost, col_size, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE ) ;
	}
	if ( rank > 0 ) {
		MPI_Send( left_real, col_size, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD ) ;
	}

}
