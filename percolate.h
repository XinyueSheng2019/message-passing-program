#include <mpi.h>

/*System size L*/

#define L 288

struct topology
{
    //this struct includes the dimensions, periods and coords of this topology, and the rank for each process and the number of processes.

    int rank,p;
    int dims[2];
    int periods[2];
    int coord[2]; 
};


struct neighbor
{
    //Each process' neighbors are put into a struct for the sake of convenience.

    int right_nbr, left_nbr, up_nbr, down_nbr;
};

struct mapinfo
{
    //this struct includes all related variables for the grid.

    int seed, nhole, step, maxstep, oldval, newval, nchange, printfreq, sumchange;
    int itop, ibot, percolate_success;
    int extension;
    int NL,M,N;
    double rho;
    double sumdigit, smalldigit;
    double random_number;
    
};

struct dynamicmaps
{
    //this struct includes all arralloc maps.

    int **map;
    int **finalmap;
    int **smallmap;
    int **demomap;
    int **old;
    int **new;
};

struct communication
{
    //this struct includes the communication related information.

    MPI_Request request1,request2,request3,request4,request5;
    MPI_Status status;
};


void set_parameters(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr,  MPI_Comm* ring_comm);
int  init_map(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, int* argc, char** argv);
void deliver_smallmap_to_old(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps);
void divide_map(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps);
void percolate_process(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm *ring_comm);
void extension_swap(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm* ring_comm, int* subarray_down_send_exten, int* subarray_down_recv_exten, int* subarray_up_send, int* subarray_up_recv);
void swap_halos(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm *ring_comm, int* subarray_right_send, int* subarray_right_recv, int* subarray_left_send, int* subarray_left_recv, int* subarray_up_send, int* subarray_up_recv, int* subarray_down_send, int* subarray_down_recv, int* subarray_down_send_exten, int* subarray_down_recv_exten);
void compare_squares(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps);
void send_smallmap(struct topology* topo, struct dynamicmaps* maps, struct communication *comm, MPI_Datatype* smallmap_type,MPI_Comm* ring_comm);
void combine_smallmap(struct topology* topo, struct dynamicmaps* maps, struct communication *comm, struct mapinfo* info, MPI_Datatype* smallmap_type, MPI_Comm* ring_comm);
void drop_extension(struct dynamicmaps* maps);
double get_sum_of_smallmap_digit(struct topology* topo,struct mapinfo* info, struct dynamicmaps* maps);
void judge_percolate(struct mapinfo* info, struct dynamicmaps* maps);
void percwritedynamic(char *percfile, int **map, int l, int ncluster);
void percwrite(char *percfile, int map[L][L], int ncluster);
void *arralloc(size_t size, int ndim, ...);
void rinit(int ijkl); //Random numbers
float uni(void);
