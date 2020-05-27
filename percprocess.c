#include <stdio.h>
#include <stdlib.h>
#include "percolate.h"
#include <mpi.h>


void percolate_process(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm *ring_comm)
{
  //this function includes the whole process of percolating map.

  int i,j;
  double start_time, end_time;

  info->maxstep = 16*L;
  info->printfreq = 100;
  info->step = 1;
  info->nchange = 1;

  //declare send and receive arrays
  int subarray_right_send[info->M+2];
  int subarray_right_recv[info->M+2];
  int subarray_left_send[info->M+2];
  int subarray_left_recv[info->M+2];
  int subarray_up_send[info->N+2];
  int subarray_up_recv[info->N+2];
  int subarray_down_send[info->N+2];
  int subarray_down_recv[info->N+2];
  int subarray_down_send_exten[info->N+2];
  int subarray_down_recv_exten[info->N+2];

  //record this while loop time
  MPI_Barrier(*ring_comm);
  start_time = MPI_Wtime();

  while (info->step <= info->maxstep)
    {
      int i,j;
       info->nchange = 0;
       info->smalldigit = 0.00;
       info->sumdigit = 0.00;
       swap_halos(topo, info, maps, nbr, comm, ring_comm, subarray_right_send, subarray_right_recv, subarray_left_send, subarray_left_recv, subarray_up_send, subarray_up_recv, subarray_down_send, subarray_down_recv, subarray_down_send_exten, subarray_down_recv_exten);
       
       //calculate the sum of digit in each processor
       info->smalldigit = get_sum_of_smallmap_digit(topo, info, maps);
      
       MPI_Allreduce(&info->nchange,&info->sumchange,1,MPI_INT, MPI_SUM,*ring_comm);
       MPI_Allreduce(&info->smalldigit,&info->sumdigit,1,MPI_DOUBLE, MPI_SUM,*ring_comm);
  
    //calculate the number of changes, average of the map for every 100 step in rank 0 processor.
    if (topo->rank == 0)
    {  
      
      if (info->step % info->printfreq == 0)
      {
        printf("percolate: number of changes on step %d is %d\n", info->step, info->sumchange);
        printf("the average of this map is %.2lf\n",(double)info->sumdigit/(L*L)); 
      }
      else if (info->step % info->printfreq != 0 && info->sumchange == 0) //info->step % info->printfreq != 0 && info->sumchange == 0
      {
        printf("percolate: number of changes on step %d is %d\n", info->step, info->sumchange);
        printf("the average of this map is %.2lf\n",(double)info->sumdigit/(L*L));
      }
      
    }

    //refresh the whole map
    for (i=1; i<=info->M; i++)
	  {
      for (j=1; j<=info->N; j++)
        {
          maps->old[i][j] = maps->new[i][j];
        }
	  
    }
 
    //if there is no change in this step, break the loop
    if(info->sumchange==0){break;}
    info->step++;

  }
  //recorde the end time of this loop, and calculate the average time for each step
   MPI_Barrier(*ring_comm);
   end_time = MPI_Wtime();

   if (topo->rank == 0) printf("the average time for each step is %.6lf \n",(double)(end_time-start_time)/(double)info->step);
    
   if (info->sumchange != 0) printf("percolate: WARNING max steps = %d reached before nchange = 0\n", info->maxstep);
    
    //deliver the final changed old map to smallmap
    for (i=1; i<=info->M; i++)
    {
        for (j=1; j<=info->N; j++){maps->smallmap[i-1][j-1] = maps->old[i][j];}
    }


}

void swap_halos(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm *ring_comm, int* subarray_right_send, int* subarray_right_recv, int* subarray_left_send, int* subarray_left_recv, int* subarray_up_send, int* subarray_up_recv, int* subarray_down_send, int* subarray_down_recv, int* subarray_down_send_exten, int* subarray_down_recv_exten)
{
  //swap the old map' halos with its neighbors by using non-blocking communication methods.
    int i,j;
    
    //initialize send and receive arrays
    for(i = 0; i < info->N+2; i++)
       {
          subarray_up_send[i] = maps->old[1][i];
          subarray_down_send[i] = maps->old[info->M][i];
          subarray_up_recv[i] = 0;
          subarray_down_recv[i] = 0;
       }
      for(j = 0; j < info->M+2; j++)
       {
          subarray_left_send[j] = maps->old[j][1];
          subarray_right_send[j] = maps->old[j][info->N];
          subarray_left_recv[j] = 0;
          subarray_right_recv[j] = 0;    
       }

      //up neighbor
      MPI_Issend(&subarray_up_send[0],info->N+2,MPI_INT,nbr->up_nbr,1,*ring_comm,&comm->request1);
      MPI_Irecv(&subarray_up_recv[0],info->N+2,MPI_INT,nbr->up_nbr,0,*ring_comm,&comm->request2);
      //down neighbor
      MPI_Issend(&subarray_down_send[0],info->N+2,MPI_INT,nbr->down_nbr,0,*ring_comm,&comm->request2);       
      MPI_Irecv(&subarray_down_recv[0],info->N+2,MPI_INT,nbr->down_nbr,1,*ring_comm,&comm->request1);
      //right neighbor
      MPI_Issend(&subarray_right_send[0],info->M+2,MPI_INT,nbr->right_nbr,2,*ring_comm,&comm->request3); 
      MPI_Irecv(&subarray_right_recv[0],info->M+2,MPI_INT,nbr->right_nbr,3,*ring_comm,&comm->request4);
      //left neighbor
      MPI_Issend(&subarray_left_send[0],info->M+2,MPI_INT,nbr->left_nbr,3,*ring_comm,&comm->request4); 
      MPI_Irecv(&subarray_left_recv[0],info->M+2,MPI_INT,nbr->left_nbr,2,*ring_comm,&comm->request3);

      MPI_Wait(&comm->request1,&comm->status);  
      MPI_Wait(&comm->request2,&comm->status);  
      MPI_Wait(&comm->request3,&comm->status);  
      MPI_Wait(&comm->request4,&comm->status); 

      extension_swap(topo, info, maps, nbr, comm, ring_comm, subarray_down_send_exten, subarray_down_recv_exten, subarray_up_send, subarray_up_recv);
      compare_squares(topo, info, maps);

    
    //put receive arrays into the responding positions of old map
      for (i = 0;i< info->N+2; i++)
      { 
        maps->old[0][i] = subarray_up_recv[i];
        maps->old[info->M+1][i] = subarray_down_recv[i];
      }
     
      for (j = 0;j<info->M+2;j++)
      {
        maps->old[j][0] = subarray_left_recv[j];
        maps->old[j][info->N+1] = subarray_right_recv[j];
      }

      
}
void extension_swap(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, struct neighbor* nbr, struct communication *comm, MPI_Comm *ring_comm, int* subarray_down_send_exten, int* subarray_down_recv_exten, int* subarray_up_send, int* subarray_up_recv)
{
    int i,j;

    if (info->extension!=0)
    {
         /*
        if this grid is extended, there are two conditions: 
        1. dims[0] == 1, such as 7,11,31: dims[0] == 1, dims[1] == p.
            by test, in this condition, M is always divided by p. There are also 2 conditions needed to be considered:
            if M > extension, the last real grid row is in the processor rank p-1.
            if M <= extension, the last real grid row is in the processor rank p-1-[extension/M]
        2. dims[0] != 1, and p has non-1 common divisors, such as 35(5*7).
        if M > extension, the last real grid row is in the processors whose rank >= dims[0]*(dims[1]-1)
        if M <= extension, the last real grid row is in the processors whose rank between (dims[1]-1-[extension/M])*dims[0] and (dims[1]-[extension/M])*dims[0]-1
        */
    
        if(info->N==info->NL)
        { 
            if (topo->rank == topo->p-1-info->extension/info->M)
            {
                for(i = 0;i<info->N+2;i++){subarray_down_send_exten[i] = maps->old[info->M-info->extension%info->M][i];subarray_down_recv_exten[i]=0;}
                MPI_Issend(&subarray_down_send_exten[0],info->N+2,MPI_INT,0,4,*ring_comm,&comm->request5);
                MPI_Irecv(&subarray_down_recv_exten[0],info->N+2,MPI_INT,0,4,*ring_comm,&comm->request5);
                MPI_Wait(&comm->request5,&comm->status);
                for(i = 0;i<info->N+2;i++){maps->old[info->M-info->extension%info->M+1][i]=subarray_down_recv_exten[i];}
            } 
            if(topo->rank == 0)
            {
                MPI_Issend(&subarray_up_send[0],info->N+2,MPI_INT,(topo->p-1-info->extension/info->M),4,*ring_comm,&comm->request5);
                MPI_Irecv(&subarray_up_recv[0],info->N+2,MPI_INT,(topo->p-1-info->extension/info->M),4,*ring_comm,&comm->request5);
                MPI_Wait(&comm->request5,&comm->status); 
                for(i = 0;i<info->N+2;i++){maps->old[0][i]=subarray_up_recv[i];}
            }      
        }
        else  
        { 
            if (topo->rank >=(topo->dims[0]-info->extension/info->M)*topo->dims[1]-topo->dims[1] && topo->rank <= (topo->dims[0]-info->extension/info->M)*topo->dims[1]-1)
            {
                for(i = 0;i<info->N+2;i++){subarray_down_send_exten[i] = maps->old[info->M-info->extension%info->M][i];subarray_down_recv_exten[i]=0;}
                MPI_Issend(&subarray_down_send_exten[0],info->N+2,MPI_INT,topo->rank%topo->dims[1],4,*ring_comm,&comm->request5);
                MPI_Irecv(&subarray_down_recv_exten[0],info->N+2,MPI_INT,topo->rank%topo->dims[1],4,*ring_comm,&comm->request5);
                MPI_Wait(&comm->request5,&comm->status);
                for(i = 0;i<info->N+2;i++){maps->old[info->M-info->extension%info->M+1][i]=subarray_down_recv_exten[i];}
            }
            if (topo->rank < topo->dims[1])
            {
                MPI_Issend(&subarray_up_send[0],info->N+2,MPI_INT,((topo->dims[0]-1-info->extension/info->M)*topo->dims[1] + topo->rank),4,*ring_comm,&comm->request5);
                MPI_Irecv(&subarray_up_recv[0],info->N+2,MPI_INT,((topo->dims[0]-1-info->extension/info->M)*topo->dims[1] + topo->rank),4,*ring_comm,&comm->request5);
                MPI_Wait(&comm->request5,&comm->status);
                for(i = 0;i<info->N+2;i++){maps->old[0][i]=subarray_up_recv[i];}
            }    
        }
    }           
}

void compare_squares(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps)
{
    //for each square, compare its neighbors' digits and refresh itself with the biggest digit
    //if the replaced digit is bigger than itself, record a change.
    int i,j;
    for (i=1; i<=info->M; i++)
    {
        for (j=1; j<=info->N; j++)
        {
        info->oldval = maps->old[i][j];
        info->newval = info->oldval;
          if (info->oldval != 0)
          {
          if (maps->old[i-1][j] > info->newval) info->newval = maps->old[i-1][j];
          if (maps->old[i+1][j] > info->newval) info->newval = maps->old[i+1][j];
          if (maps->old[i][j-1] > info->newval) info->newval = maps->old[i][j-1];
          if (maps->old[i][j+1] > info->newval) info->newval = maps->old[i][j+1];

            if (info->newval != info->oldval)
            {
                ++info->nchange;
            }
          }
         maps->new[i][j] = info->newval;
        }
    }
}

double get_sum_of_smallmap_digit(struct topology* topo,struct mapinfo* info, struct dynamicmaps* maps)
{
  //calculate the sum of digits on smallmap, we should consider the smallmap which includes extension part. Its M and N may include halos.
  
  int i,j;

   if (info->extension != 0)
   {
     if(info->N==info->NL)
     {
       if(topo->rank != topo->p - 1 - info->extension/info->M)
       {
         for(i = 1;i<=info->M;i++)
         {
          for(j = 1; j<=(info->N - info->extension);j++)
          {
            info->smalldigit += maps->new[i][j];
          } 
         }
       }
       else
       {
          for(i = 1;i<=(info->M - info->extension%info->M);i++)
          {
          for(j = 1; j<=(info->N - info->extension);j++)
          {
            info->smalldigit += maps->new[i][j];
          } 
          }  
       }    
     }
     else // dims[1]!=1
     {
       if(topo->rank >=(topo->dims[0]-info->extension/info->M)*topo->dims[1]-topo->dims[1] && topo->rank <= (topo->dims[0]-info->extension/info->M)*topo->dims[1]-1)
       {
         for(i = 1;i<=(info->M - info->extension%info->M);i++)
         {
          for(j = 1; j<=info->N;j++)
          {
            info->smalldigit += maps->new[i][j];
          } 
        }
       }
       else
       {
         for(i = 1;i<=info->M;i++)
         {
          for(j = 1; j<=(info->N);j++)
          {
            info->smalldigit += maps->new[i][j];
          } 
        }
         
       }
         
     }
     
   }
   else //no extension
   {
    for(i = 1;i<=info->M;i++)
    {
      for(j = 1; j<=info->N;j++)
      {
        info->smalldigit += maps->new[i][j];
      } 
    }
   }

 
  return info->smalldigit;
}