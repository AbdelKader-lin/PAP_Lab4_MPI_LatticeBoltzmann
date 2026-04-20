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
//      8 neighbors using MPI types for non contiguous
//      side.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - Blocking communications
// NEW:
//     - >>> MPI type for non contiguous cells <<<
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"

/*
count
    number of blocks (nonnegative integer) 
blocklength
    number of elements in each block (nonnegative integer) 
stride
    number of elements between start of each block (integer) 
oldtype
    old datatype (handle) 
*/
/****************************************************/
void lbm_comm_init_ex5( lbm_comm_t * comm , int total_width , int total_height ) {
	//we use the same implementation than ex5 execpt for type creation
	lbm_comm_init_ex4( comm , total_width , total_height ) ;

	//TODO: create MPI type for non contiguous side in comm->type.
	MPI_Type_vector( comm->width , DIRECTIONS , comm->height * DIRECTIONS , MPI_DOUBLE , &comm->type ) ; //  Creates a MPI vector type.
	MPI_Type_commit( &comm->type ) ;

}




/****************************************************/
void lbm_comm_release_ex5( lbm_comm_t * comm ) {
	//we use the same implementation than ex5 except for type destroy
	lbm_comm_release_ex4( comm );

	//TODO: release MPI type created in init.
	MPI_Type_free( &comm->type ) ;


}

/****************************************************/
void lbm_comm_ghost_exchange_ex5( lbm_comm_t * comm , lbm_mesh_t * mesh ) {
	//
	// TODO: Implement the 2D communication with :
	//         - blocking MPI functions
	//         - use MPI type for non contiguous side 
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


	MPI_Sendrecv( lbm_mesh_get_cell( mesh , 0 , 1 ) , 1 , comm->type , uprank , 2 , 
		lbm_mesh_get_cell( mesh , 0 , comm->height - 1 ) , 1 , comm->type , drank , 2 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;

	MPI_Sendrecv( lbm_mesh_get_cell( mesh , 0 , comm->height - 2 ) , 1 , comm->type , drank , 3 ,
        lbm_mesh_get_cell( mesh , 0 , 0 ) , 1 , comm->type , uprank , 3 ,
        MPI_COMM_WORLD , MPI_STATUS_IGNORE ) ;


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
