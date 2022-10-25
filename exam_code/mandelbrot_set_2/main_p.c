#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <complex.h> 
#include "zpoint.h"
#include "generate_matrix.h"
#include "checkMB.h"
#include "save_on_file.h"
#define DOUBLE_PRECISION

// it an alis for double
// check if it is defined or not! you can do it in main or from compiler using -D flag


int main (int argc, char* argv[])
{
   #ifdef _OPENMP
    printf("I activate the OPENMP paradigm! Have fun!\n");
   #endif
    
   int Nthreads   = atoi(argv[8]);
   int size_x     = atoi(argv[1]); 
   int size_y     = atoi(argv[2]); 
   float_prec xL  = atof(argv[3]); 
   float_prec yL  = atof(argv[4]); 
   float_prec xR  = atof(argv[5]); 
   float_prec yR  = atof(argv[6]); 
   int I_max      = atoi(argv[7]);
   
   zpoint* data = (struct zpoint*)malloc(size_x*size_y*sizeof( struct zpoint));
   generate_matrix( data,  size_x,  size_y,  xL, yL,  xR,  yR);
   
    for (int k=1;k<=Nthreads; k++) 
   {  
      //time variables
      float_prec start=0;
      float_prec end =0;
      start = omp_get_wtime(); 
      #pragma omp parallel num_threads(k) 
      {  
         #pragma omp for 
         for(int i=0; i<size_x; i++)
         {
            for(int j=0; j<size_y; j++)
            {  
               checkMB(&data[i*size_x+j], I_max);
            }// end for 
         }
 
      } // end parallel region
      end = omp_get_wtime(); 
      printf("\n for %d threads the needing time is: %f \n", k, end-start);
      
   } // end for nthread
   save_on_file( data, size_x, size_y);
   free(data);
   return 0;
}