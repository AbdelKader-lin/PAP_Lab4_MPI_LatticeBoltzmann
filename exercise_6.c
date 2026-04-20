/*****************************************************
    AUTHOR  : Sébastien Valat
    MAIL    : sebastien.valat@univ-grenoble-alpes.fr
    LICENSE : BSD
    YEAR    : 2021
    COURSE  : Parallel Algorithms and Programming
*****************************************************/

//////////////////////////////////////////////////////
//
// Goal: Implement 2D grid communication with non-blocking
//       messages.
//
// SUMMARY:
//     - 2D splitting along X and Y
//     - 8 neighbors communications
//     - MPI type for non contiguous cells
// NEW:
//     - Non-blocking communications
//
//////////////////////////////////////////////////////

/****************************************************/
#include "src/lbm_struct.h"
#include "src/exercises.h"


static int get_rank_2d(lbm_comm_t *comm, int x, int y)
{
	// return MPI_PROC_NULL if outside domain
	if (x < 0 || x >= comm->nb_x || y < 0 || y >= comm->nb_y)
		return MPI_PROC_NULL;

	int coords[2] = {x, y};
	int rank;
	MPI_Cart_rank(comm->communicator, coords, &rank);
	return rank;
}

/****************************************************/
void lbm_comm_init_ex6(lbm_comm_t * comm, int total_width, int total_height)
{
	//we use the same implementation than ex5
	lbm_comm_init_ex5(comm, total_width, total_height);
}

/****************************************************/
void lbm_comm_release_ex6(lbm_comm_t * comm)
{
	//we use the same implementation than ext 5
	lbm_comm_release_ex5(comm);
}

/****************************************************/
void lbm_comm_ghost_exchange_ex6(lbm_comm_t * comm, lbm_mesh_t * mesh)
{
	//
	// TODO: Implement the 2D communication with :
	//         - non-blocking MPI functions
	//         - use MPI type for non contiguous side 
	//
	// TIP: The previous trick require to make two batch of non-blocking communications.
    MPI_Request reqs[16];
    int r = 0;

	// To be used:
	//    - DIRECTIONS: the number of doubles composing a cell
	//    - double[9] lbm_mesh_get_cell(mesh, x, y): function to get the address of a particular cell.
	//    - comm->width : The with of the local sub-domain (containing the ghost cells)
	//    - comm->height : The height of the local sub-domain (containing the ghost cells)
	//
	int left, right, top, bottom;
    //direct neighbors
    MPI_Cart_shift(comm->communicator, 0, 1, &left, &right);   // X direction
    MPI_Cart_shift(comm->communicator, 1, 1, &top, &bottom);   // Y direction

    int w = comm->width;
    int h = comm->height;

	// TIP: create a function to get the target rank from x,y task coordinate.
	// TIP: You can use MPI_PROC_NULL on borders.
	// TIP: send the corner values 2 times, with the up/down/left/right communication // right, right ? :)))
	//      and with the diagonal communication in a second time, this avoid
	//      special cases for border tasks.

	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 1), h - 2, comm->type,left, 0, comm->communicator, &reqs[r++]);
    MPI_Irecv(lbm_mesh_get_cell(mesh, w - 1, 1), h - 2, comm->type, right, 0, comm->communicator, &reqs[r++]);
    MPI_Irecv(lbm_mesh_get_cell(mesh, 1, 0), (w - 2) * DIRECTIONS, MPI_DOUBLE,top, 1, comm->communicator, &reqs[r++]);
    MPI_Irecv(lbm_mesh_get_cell(mesh, 1, h - 1), (w - 2) * DIRECTIONS, MPI_DOUBLE,bottom, 1, comm->communicator, &reqs[r++]);




	//example to access cell
	//double * cell = lbm_mesh_get_cell(mesh, local_x, local_y);
	//double * cell = lbm_mesh_get_cell(mesh, comm->width - 1, 0);

	//TODO:
	//   - implement left/write communications
	//   - implement top/bottom communication (non contiguous)
	//   - implement diagonal communications


	// send to left right top bottom
	// ()
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 1), h - 2, comm->type,left, 0, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, w - 2, 1), h - 2, comm->type,right, 0, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 1), (w - 2) * DIRECTIONS, MPI_DOUBLE,top, 1, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, h - 2), (w - 2) * DIRECTIONS, MPI_DOUBLE,bottom, 1, comm->communicator, &reqs[r++]);

	MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);


	// for diagonal : 
	r = 0;
	int coords[2];
	MPI_Cart_coords(comm->communicator, comm->rank_x, 2, coords);

	int diag_rank;
// for up left
	diag_rank = get_rank_2d(comm, coords[0] - 1, coords[1] - 1);

	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, 0), DIRECTIONS, MPI_DOUBLE,diag_rank, 2, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, 1), DIRECTIONS, MPI_DOUBLE, diag_rank, 2, comm->communicator, &reqs[r++]);

// for up right
	diag_rank = get_rank_2d(comm, coords[0] + 1, coords[1] - 1);

	MPI_Irecv(lbm_mesh_get_cell(mesh, w - 1, 0), DIRECTIONS, MPI_DOUBLE,diag_rank, 3, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, w - 2, 1), DIRECTIONS, MPI_DOUBLE,diag_rank, 3, comm->communicator, &reqs[r++]);


	// for down left
	diag_rank = get_rank_2d(comm, coords[0] - 1, coords[1] + 1);	

	MPI_Irecv(lbm_mesh_get_cell(mesh, 0, h - 1), DIRECTIONS, MPI_DOUBLE,diag_rank, 4, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, 1, h - 2), DIRECTIONS, MPI_DOUBLE,diag_rank, 4, comm->communicator, &reqs[r++]);


	// for down right
	diag_rank = get_rank_2d(comm, coords[0] + 1, coords[1] + 1);
	MPI_Irecv(lbm_mesh_get_cell(mesh, w - 1, h - 1), DIRECTIONS, MPI_DOUBLE,diag_rank, 5, comm->communicator, &reqs[r++]);
	MPI_Isend(lbm_mesh_get_cell(mesh, w - 2, h - 2), DIRECTIONS, MPI_DOUBLE,diag_rank, 5, comm->communicator, &reqs[r++]);


	// wait for all communications to end : 
	MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);






}
