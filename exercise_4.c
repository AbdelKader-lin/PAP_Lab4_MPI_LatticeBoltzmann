/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication scheme with
//       8 neighbors using manual copy for non
//       contiguous side and blocking communications
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
//     - Manual copy for non continguous cells
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

#include <string.h>
#include <stdlib.h>

/****************************************************/

/*
According to Lab3 :
use MPI_Dims_create to define the dimensions of the grid,
and MPI_Cart_create to create the grid
*/
void lbm_comm_init_ex4( lbm_comm_t * comm , int total_width , int total_height ) {
	//
	// TODO: calculate the splitting parameters for the current task.
	//

	// Get info
	int rank ;
	int comm_size ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank ) ;
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size ) ;


	// TODO: calculate the number of tasks along X axis and Y axis.
	// int MPI_Dims_create(int nnodes, int ndims, int dims[])
	int dims[ 2 ] = { 0 , 0 } ;
	MPI_Dims_create( comm_size , 2 , dims ) ; // We def the dimensions of the grid
	comm->nb_x = dims[ 0 ] ;
	comm->nb_y = dims[ 1 ] ;

	// TODO: calculate the current task position in the splitting
	comm->rank_x = rank % dims[ 0 ] ;
	comm->rank_y = rank / dims[ 0 ] ;

	// TODO : calculate the local sub-domain size (do not forget the 
	//        ghost cells). Use total_width & total_height as starting 
	//        point.
	comm->width = ( total_width / dims[ 0 ] ) + 2 ;
	comm->height = ( total_height / dims[ 1 ] ) + 2 ;

	// TODO : calculate the absolute position  (in cell number) in the global mesh.
	//        without accounting the ghost cells
	//        (used to setup the obstable & initial conditions).
	// x , y : rank position
	comm->x = comm->rank_x * ( total_width / dims[ 0 ] ) ;
	comm->y = comm->rank_y * ( total_height / dims[ 1 ] ) ;

	//OPTIONAL : if you want to avoid allocating temporary copy buffer
	//           for every step :
	//comm->buffer_recv_down, comm->buffer_recv_up, comm->buffer_send_down, comm->buffer_send_up

	// comm->width * 9 * sizeof( double ) because :
	// nb of cells per row * 9 ( number of values per cell ) * sizeof( double )
	comm->buffer_recv_down = malloc ( comm->width * 9 * sizeof( double ) ) ;
	comm->buffer_recv_up = malloc ( comm->width * 9 * sizeof( double ) ) ;

	comm->buffer_send_down = malloc ( comm->width * 9 * sizeof( double ) ) ;
	comm->buffer_send_up= malloc ( comm->width * 9 * sizeof( double ) ) ;

	//if debug print comm
	//lbm_comm_print(comm);
}

/****************************************************/
void lbm_comm_release_ex4( lbm_comm_t * comm ) {
	//free allocated ressources
	free( comm->buffer_recv_down ) ;
	free( comm->buffer_recv_up ) ; 
	free( comm->buffer_send_down ) ;
	free( comm->buffer_send_up ) ;
}



int get_rank( lbm_comm_t * comm , int x , int y ){

	// Fct to get the rank of the process

	if ( x < 0 || y < 0 || x >= comm->nb_x || y >= comm->nb_y ) {
		return MPI_PROC_NULL ;
	}
	
	/*
	 *		x=0   x=1   x=2   x=3
	 *	y=0 [  P0  |  P1  |  P2  |  P3  ]
	 *	y=1 [  P4  |  P5  |  P6  |  P7  ]
	 */
	// r = x + y * nombre de cols
	return ( x + y * comm->nb_x ) ;

}

/****************************************************/
void lbm_comm_ghost_exchange_ex4( lbm_comm_t * comm , lbm_mesh_t * mesh ) {
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - manual copy in temp buffer for non contiguous side 
	//
	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	// TIP: create a function to get the target rank from x,y task coordinate. 
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/write communication
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications

	// ###############################

	// left / write communications

	// Number of doubles in one col
	int col_size = comm->height * DIRECTIONS ;

	double * left_ghost = lbm_mesh_get_cell( mesh , 0 , 0 ) ;
	double * left_real  = lbm_mesh_get_cell( mesh , 1 , 0 ) ;
	double * right_real = lbm_mesh_get_cell( mesh , comm->width - 2 , 0 ) ;
	double * right_ghost= lbm_mesh_get_cell( mesh , comm->width - 1 , 0 ) ;


	/*
	 *  Now we compute the ranks for the 4 cells around. 
	 *	The reason we do it now and not in 1D, is the boundaries.
	 * 	
	 * 
	 *		x=0   x=1   x=2   x=3
	 *	y=0 [  P0  |  P1  |  P2  |  P3  ]
	 *	y=1 [  P4  |  P5  |  P6  |  P7  ]
	 */
	int lrank = get_rank( comm , comm->rank_x - 1 , comm->rank_y ) ;
	int r_rank = get_rank( comm , comm->rank_x + 1 , comm->rank_y ) ;
	int uprank = get_rank( comm , comm->rank_x , comm->rank_y - 1 ) ;
	int drank = get_rank( comm , comm->rank_x , comm->rank_y + 1 ) ;


	MPI_Sendrecv( left_real , col_size , MPI_DOUBLE , lrank , 0 , 
		right_ghost , col_size , MPI_DOUBLE , r_rank , 0 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( right_real , col_size , MPI_DOUBLE , r_rank , 1 ,
        left_ghost , col_size , MPI_DOUBLE , lrank , 1 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;


	// up / down communications
	int row_size = comm->width * DIRECTIONS ;
	/*
	To send non-contiguous data, in this exercise you will have to use a temporary buffer into which
	you will manually copy the non-contiguous data before sending. The data will also be received
	in a temporary buffer and will be then copied into the correct ghost cells.
	*/

	/*double * up_ghost = lbm_mesh_get_cell( mesh , 0 , 0 ) ;
	double * up_real  = lbm_mesh_get_cell( mesh , 1 , 0 ) ;
	double * down_real = lbm_mesh_get_cell( mesh , comm->width - 2 , 0 ) ;
	double * down_ghost= lbm_mesh_get_cell( mesh , comm->width - 1 , 0 ) ;
	*/

	// To send 

	// memcpy( *to , *from , numBytes ) ;
	for ( int i = 0 ; i < comm->width ; i++ ) {
		memcpy( &comm->buffer_send_up[ i * DIRECTIONS ] , lbm_mesh_get_cell( mesh , i , 1 ) , DIRECTIONS * sizeof( double ) ) ;
	}

	for ( int i = 0 ; i < comm->width ; i++ ) {
		memcpy( &comm->buffer_send_down[ i * DIRECTIONS ] , lbm_mesh_get_cell( mesh , i , comm->height - 2 ) , DIRECTIONS * sizeof( double ) ) ;
	}


	MPI_Sendrecv( comm->buffer_send_up , row_size , MPI_DOUBLE , uprank , 2 , 
		comm->buffer_recv_down , row_size , MPI_DOUBLE , drank , 2 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( comm->buffer_send_down , row_size , MPI_DOUBLE , drank , 3 ,
        comm->buffer_recv_up , row_size , MPI_DOUBLE , uprank , 3 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	// To receive

	// memcpy( *to , *from , numBytes ) ;
	for ( int i = 0 ; i < comm->width ; i++ ) {
		memcpy( lbm_mesh_get_cell( mesh , i , 0 ) , &comm->buffer_recv_up[ i * DIRECTIONS ] , DIRECTIONS * sizeof( double ) ) ;
	}

	for ( int i = 0 ; i < comm->width ; i++ ) {
		memcpy( lbm_mesh_get_cell( mesh , i , comm->height - 1 ) , &comm->buffer_recv_down[ i * DIRECTIONS ] , DIRECTIONS * sizeof( double ) ) ;
	}


	// ###############################
	// ###############################
	// ###############################

	// Now, the diagonales

	int up_l_rank = get_rank( comm , comm->rank_x - 1 , comm->rank_y - 1 ) ;
	int up_r_rank = get_rank( comm , comm->rank_x + 1 , comm->rank_y - 1 ) ;
	int down_l_rank = get_rank( comm , comm->rank_x - 1 , comm->rank_y + 1 ) ;
	int down_r_drank = get_rank( comm , comm->rank_x + 1 , comm->rank_y + 1 ) ;


	double * up_left_ghost = lbm_mesh_get_cell( mesh , 0 , 0 ) ;
	double * up_left_real  = lbm_mesh_get_cell( mesh , 1 , 1 ) ;
	double * up_right_real = lbm_mesh_get_cell( mesh , comm->width - 2 , 1 ) ;
	double * up_right_ghost= lbm_mesh_get_cell( mesh , comm->width - 1 , 0 ) ;

	double * down_left_ghost = lbm_mesh_get_cell( mesh , 0 , comm->height - 1 ) ;
	double * down_left_real  = lbm_mesh_get_cell( mesh , 1 , comm->height - 2 ) ;
	double * down_right_real = lbm_mesh_get_cell( mesh , comm->width - 2 , comm->height - 2 ) ;
	double * down_right_ghost= lbm_mesh_get_cell( mesh , comm->width - 1 , comm->height - 1 ) ;

	
	MPI_Sendrecv( up_left_real , DIRECTIONS , MPI_DOUBLE , up_l_rank , 4 , 
		down_right_ghost , DIRECTIONS , MPI_DOUBLE , down_r_drank , 4 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( up_right_real , DIRECTIONS , MPI_DOUBLE , up_r_rank , 5 ,
        down_left_ghost , DIRECTIONS , MPI_DOUBLE , down_l_rank , 5 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( down_left_real , DIRECTIONS , MPI_DOUBLE , down_l_rank , 6 , 
		up_right_ghost , DIRECTIONS , MPI_DOUBLE ,up_r_rank , 6 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( down_right_real , DIRECTIONS , MPI_DOUBLE , down_r_drank , 7 ,
        up_left_ghost , DIRECTIONS , MPI_DOUBLE , up_l_rank , 7 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

}
