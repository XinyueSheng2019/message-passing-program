#include <stdio.h>
#include <stdlib.h>
#include "percolate.h"
#include <mpi.h>

void set_parameters(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, MPI_Comm* ring_comm)
{
  topo->periods[0] = 1;
  topo->periods[1] = 0;
  topo->dims[0] = 0;
  topo->dims[1] = 0;

  MPI_Comm_rank(MPI_COMM_WORLD,&topo->rank);
  MPI_Comm_size(MPI_COMM_WORLD,&topo->p);
  MPI_Dims_create(topo->p,2,topo->dims);
  MPI_Cart_create(MPI_COMM_WORLD,2,topo->dims,topo->periods,0,ring_comm);
  MPI_Cart_shift(*ring_comm,0,1,&nbr->up_nbr,&nbr->down_nbr);
  MPI_Cart_shift(*ring_comm,1,1,&nbr->left_nbr,&nbr->right_nbr);
  MPI_Cart_coords(*ring_comm, topo->rank, 2, topo->coord); 
  
  info->extension = 0;
  info->NL = L;

  if (L%topo->dims[0]!=0||L%topo->dims[1]!=0) 
  {
    info->extension = topo->p-(L%topo->p);
    info->NL = L + info->extension;
  }

  info->M = info->NL/topo->dims[0];
  info->N = info->NL/topo->dims[1];
 
  maps->map = arralloc(sizeof(int),2,info->NL,info->NL);
  maps->finalmap = arralloc(sizeof(int),2,L,L);
  maps->demomap = arralloc(sizeof(int),2,info->M,info->N);
  maps->smallmap = arralloc(sizeof(int),2,info->M,info->N);
  maps->old = arralloc(sizeof(int),2,info->M+2,info->N+2);
  maps->new= arralloc(sizeof(int),2,info->M+2,info->N+2);
   
 
  
  
  
    
}