#if !defined(DOUBLE_PRECISION)  
#define float_prec float
#else
#define float_prec double 
#endif
// it an alis for double
// check if it is defined or not! you can do it in main or from co piler using -d


#ifndef ZPOINT
#define ZPOINT
 typedef struct zpoint
    {   
        double real;
        double im;
        int counter;
    } zpoint;
#endif