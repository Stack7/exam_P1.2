#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include "zpoint.h"
#include "checkMB.h"

void  checkMB(struct zpoint* data, int I_max)
{
    int interaction =0;
    double xn=0, yn =0;
    double xn2, yn2;
 
   while( (interaction < I_max)&&(xn*xn+yn*yn <4) )
   {     // recursive functionm for z=x+iy;
      xn2=xn*xn - yn*yn;
      yn2=xn*yn;

      xn = data->real+ xn2;
      yn = data->im+2*yn2;
      interaction=interaction+1 ;
   }
   data->counter= (int) interaction/I_max; 
}
