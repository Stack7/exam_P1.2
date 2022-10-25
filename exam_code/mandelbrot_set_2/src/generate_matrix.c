#include <stdio.h>
#include <stdlib.h>
#include "zpoint.h"
#include "generate_matrix.h"

void generate_matrix( zpoint* data, int size_x, int size_y, double xL,double yL, double xR, double yR)
{  double delta_x = (xR-xL)/size_x;
   double delta_y = (yR-yL)/size_y;
   #pragma omp parallel 
   {
        for(int i=0; i<size_x; i++)
        {
            for(int j=0; j<size_y; j++)
            {
                data[i*size_y+ j].real = xL +i*delta_x ;
                data[i*size_y+ j].im =  yL +j*delta_y ;
                data[i*size_y+ j].counter=0;
            }
        } //end for loop
    }

    
}