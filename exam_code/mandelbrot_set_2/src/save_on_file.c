#include <stdio.h>
#include "zpoint.h"

void save_on_file(struct zpoint* data, int size_x, int size_y)
{
    FILE* mb_set;
    mb_set = fopen("mandelbrot.pgm", "w");
      // Writing Magic Number to the File
    fprintf(mb_set, "P2\n");

      // Writing sizes over x and y
    fprintf(mb_set, "%d %d\n", size_x, size_y);
    fprintf(mb_set, "255\n");
    for(int i =0;i<size_x*size_y; i++)
    {
        fprintf(mb_set, "%d %f %f \n", data[i].counter, data[i].real, data[i].im);
    }
   
    fclose(mb_set);
}