#include <stdio.h>
#include <stdlib.h>
#include "percolate.h"

int init_map(struct topology* topo, struct mapinfo* info, struct dynamicmaps* maps, int* argc, char** argv)
{
  //the whole map is created, and when the map is extended, 0 should be set in every extended square.
  int i, j;
  
  if (topo->rank == 0)
  {
    if (*argc != 2)
      {
        printf("Usage: percolate <seed>\n");
        return 1;
      }

    info->rho = 0.411;
    info->seed = atoi(argv[1]);
    rinit(info->seed);
 
    printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, info->rho, info->seed);
    
    info->nhole = 0;

        for (i=0; i < L; i++)
        {
          for (j=0; j < L; j++)
          {
            info->random_number=uni();
          
            if(info->random_number < info->rho)
            {
              maps->map[i][j] = 0;
            }
            else
            {
              info->nhole++;
              maps->map[i][j] = info->nhole;
            }
          }
        }
        //if the extension doesn't equal 0, we need to add 0 to extended squares.
        if (info->extension!= 0)
        {
          for(i = L;i < info->NL; i++)
          {
            for(j = L;j < info->NL; j++)
            {
              maps->map[i][j] = info->nhole;
            }
          }
        }
    printf("percolate: rho = %f, actual density = %f\n", info->rho, 1.0 - ((double) info->nhole)/((double) L*L) );
    return 0;

  }
}