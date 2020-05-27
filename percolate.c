#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "percolate.h"

int main(int argc, char *argv[])
{ 
  MPI_Init(NULL,NULL);
  
  int i, j;
  struct topology *topo, topo1; topo = &topo1;
  struct dynamicmaps *maps, maps1; maps = &maps1;
  struct mapinfo *info, mapinfo1; info = &mapinfo1;
  struct neighbor *nbr, nbr1;nbr = &nbr1;
  struct communication *comm, comm1; comm = &comm1;
  MPI_Comm ring_comm;

  set_parameters(topo, info, maps, nbr, &ring_comm);
  init_map(topo, info, maps, &argc, argv);
  MPI_Bcast(&(maps->map[0][0]),info->NL*info->NL,MPI_INT,0,ring_comm); //Broadcast map to each processor
  divide_map(topo, info, maps);
  deliver_smallmap_to_old(topo, info, maps);
  percolate_process(topo, info, maps, nbr, comm, &ring_comm);

 //this vector is used to send smallmap to rank 0 processor
  MPI_Datatype smallmap_type;
  MPI_Type_vector(info->M*info->N,1,1,MPI_INT,&smallmap_type);
  MPI_Type_commit(&smallmap_type);
  send_smallmap(topo, maps, comm, &smallmap_type, &ring_comm);

  if (topo->rank == 0)
  {
    combine_smallmap(topo, maps, comm, info, &smallmap_type,&ring_comm);
    free(maps->smallmap);
    drop_extension(maps);
    judge_percolate(info, maps);
    percwritedynamic("map.pgm", maps->finalmap, L, 3);
    free(maps->map);
  }
   
  MPI_Finalize();
  return 0;
}

void divide_map(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps)
{
   //divide map into corresponding smallmap
  int m,n,i,j;
  m = 0;
  for (i = topo->coord[0]*info->M;i<(topo->coord[0]*info->M+info->M);i++)
  {
    n = 0;
    for(j = topo->coord[1]*info->N;j<(topo->coord[1]*info->N+info->N);j++)
    {  
      maps->smallmap[m][n] = maps->map[i][j];
      n++;
    }  
    m++;
  }
}

void send_smallmap(struct topology* topo, struct dynamicmaps* maps, struct communication* comm, MPI_Datatype* smallmap_type, MPI_Comm* ring_comm)
{
  //send the small map to rank 0 process
  int i;
  for(i = 0;i<topo->p;i++)
  {
    if (topo->rank == i)
    {
       MPI_Isend(&maps->smallmap[0][0],1,*smallmap_type,0,0,*ring_comm, &comm->request5); 
    }
  }
}
void combine_smallmap(struct topology* topo, struct dynamicmaps* maps, struct communication *comm, struct mapinfo* info, MPI_Datatype* smallmap_type, MPI_Comm* ring_comm)
{
  //receive the smallmap from each processor and combine them into map(with extension)
  int m,n,i,j,t,k,s;
  s = 0;
  for(i = 0;i<topo->dims[0];i++)
  {
    for(j = 0;j<topo->dims[1];j++)
    {
      MPI_Recv(&maps->demomap[0][0],1,*smallmap_type,s,0,*ring_comm,&comm->status); 
      m = 0;
      for(t = i*info->M;t<i*info->M+info->M;t++){
        n = 0;
        for(k = j*info->N;k<j*info->N+info->N;k++){
          maps->map[t][k] = maps->demomap[m][n];
          n++;
        }
        m++;
      }
      s++;  
    }
  }
  


}

void drop_extension(struct dynamicmaps* maps)
{
  //drop the extension (NL-L) part of the map
  
  int i,j;
  for(i = 0;i<L;i++)
  {
    for(j = 0;j<L;j++)
    {
      maps->finalmap[i][j] = maps->map[i][j];
    }
  }
  
}

void judge_percolate(struct mapinfo* info, struct dynamicmaps* maps)
{
  //judge whether this map is percolated
  
  info->percolate_success = 0;

  for (info->itop=0; info->itop < L; info->itop++)
  {
    if (maps->finalmap[info->itop][L-1] > 0)
    {
      for (info->ibot=0; info->ibot < L; info->ibot++)
        {
          if (maps->finalmap[info->itop][L-1] == maps->finalmap[info->ibot][0])
          {
            info->percolate_success = 1;
          }
        }
    }
  }

  if (info->percolate_success != 0)
    {
      printf("percolate: cluster DOES percolate\n");
    }
  else
    {
      printf("percolate: cluster DOES NOT percolate\n");
    }

}