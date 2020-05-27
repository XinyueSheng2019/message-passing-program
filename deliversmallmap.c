#include <stdio.h>
#include <stdlib.h>
#include "percolate.h"

void deliver_smallmap_to_old(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps)
{
  //this function can deliver the smallmap in each process to a old map with width N+2 and length M+2, and then add halos to each old map.
  
  int i,j;

  //deliver the data of smallmap to old map
   for (i=1; i <= info->M; i++)
  {
    for (j=1; j <= info->N; j++)
    {
        maps->old[i][j] = maps->smallmap[i-1][j-1];
    }
  }
    
    //add 0 to each halo
   for (i=0; i <= info->M+1; i++)  
    {
      // zero the bottom and top halos
      maps->old[i][0]   = 0;
      maps->old[i][info->N+1] = 0;
    }

   for (j=0; j <= info->N+1; j++)  
    {
      // zero the left and right halos
      maps->old[0][j]   = 0;
      maps->old[info->M+1][j] = 0;
    }
}





