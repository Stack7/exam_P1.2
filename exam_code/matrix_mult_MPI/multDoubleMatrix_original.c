#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <cblas.h>
void print_mat(double* A, int N_LOC,int N)
{
    for(int k=0;k<N_LOC ;k++)
        {   
            for(int l=0;l<N;l++)
      {  
              printf("%f ",A[k*N+l]);
          }
      printf("\n"); 
    } //end for
  } //end for


void print_all_mat(double* A, int N_LOC,int N, int rest, int Npes, int rank)
{  
    if(!rank)
    {
        print_mat(A, N_LOC,N);//rank 0 print its lines
    
        for(int x = 1; x<Npes;x++)
        {   
            int N_LOC= (x<rest ?  N/Npes + 1 : N/Npes);
            double* print_buff= (double*)malloc(N_LOC*N*sizeof(double));

            MPI_Recv(print_buff, N*N_LOC, MPI_DOUBLE, x, x, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            print_mat(print_buff, N_LOC,N);// the p0 prints the others stuff received
            free(print_buff);
        }
    }else
    {
        MPI_Send(A, N*N_LOC, MPI_DOUBLE, 0,rank, MPI_COMM_WORLD);}//if you aren't 0 sends!
};
//end print_all_mat

//-------------------MAIN-----------------------------

int main (int argc, char* argv[])
{
    int N; 
    if (argc>1)
    {  
        N=atoi(argv[1]);
    }
    int rank, Npes; 
  
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Npes);

    int N_LOC=N/Npes;
    int rest = N % Npes;
    int* size_col =(int*)calloc(Npes , sizeof(int));
    int* size_row =(int*)calloc(Npes , sizeof(int));

    for(int i =0; i<Npes; i++)
    { 
        size_col[i] = (i<rest ?  N_LOC+1 : N_LOC);
        size_row[i] = (i<rest ?  N_LOC+1 : N_LOC);
    }

    if(rank<rest)N_LOC++;
    

    int* recvcounts = (int*)calloc(Npes, sizeof(int));
    int* displacement = (int*)calloc(Npes, sizeof(int));  

    double* A = (double*)calloc(N_LOC*N,sizeof(double));
    double* B = (double*)calloc(N_LOC * N ,sizeof(double));
    double* C = (double*)calloc(N_LOC* N, sizeof(double));
    
    for(int m=0;m<N_LOC;m++) // my process p_n knows inside all the row into it! so just traverse the row
    {
        for(int n=0;n<N;n++)
        {
            A[m*N+n]=n;//rand(rank)%10;
            B[m*N+n]=1.0/(n+1);//(rand(rank)%10);
            C[m*N+n]=0;
        }  
    }   
    if(N<6){ // just to avoid to print too big matrixes
        if(!rank)
        {
            printf("matrix A \n");
        }
        print_all_mat( A, N_LOC, N, rest,  Npes, rank);

        if(!rank) 
        {
            printf("matrixB \n");
        }
            print_all_mat( B,  N_LOC, N, rest,  Npes, rank);   
          
    }
    
    int  j_glob= 0;   
    // my process p_n knows its rows 
    // I need that every processor takes from the ithers their blocks of the column
    for(int i=0; i<Npes;i++)    
    {
        double* recvbuf = (double*)malloc( size_col[i]*N* sizeof(double));
        double* sendbuf=(double*)malloc(N_LOC *  size_col[i]* sizeof(double));

        for(int ii=0;ii<size_row[i]; ii++)  //bc there are the lines for every processor,N_LOC
        { 
            for(int jj=0;jj<N_LOC;jj++) //these are the columns
            {
                sendbuf[jj*size_row[i]+ii]=B[ii+jj*N+j_glob];
            }
        }

        int total_size=0; 
        displacement[0]=0;    
        for(int j=0;j<Npes; j++)
        {  
            //I have to initialise all the vector for every processor for each stap
            recvcounts[j]= size_col[i] *size_col[j];
            if(j>=1){
                displacement[j]=displacement[j-1]+recvcounts[j-1];   
            }
            total_size += recvcounts[j];
        }
                   
        MPI_Allgatherv(sendbuf,N_LOC*size_col[i] ,MPI_DOUBLE, 
                        recvbuf, recvcounts, displacement, MPI_DOUBLE, MPI_COMM_WORLD);  
        


        for(int ii=0;ii<N_LOC;ii++)
        {
            for(int col=0;col<size_col[i];col++)
            { 
                for(int k=0; k<N; k++)
                {
                    C[ii*N+(col+j_glob)]+=A[ii*N+k]*recvbuf[k*size_col[i]+col];
                }
            }
        }
    
        j_glob +=size_col[i];

        free(sendbuf);
        free(recvbuf);
    }  //end for 

    if(N<6){
    if(!rank)
    {
        printf("matrixC \n");
    }
    print_all_mat( C,  N_LOC, N, rest,  Npes, rank); 
    }
//-----------to check if the matrix C calculated bu myself is ok I use cblst
//-----------just to check once, then if it worked it is good!

    double* A1 = (double*)calloc(N*N,sizeof(double));
    double* B1 = (double*)calloc(N*N ,sizeof(double));
    double* C1 = (double*)calloc(N*N, sizeof(double));
    double* D1 = (double*)calloc(N*N, sizeof(double));

    if(!rank)
    {

        for(int i=0; i<N_LOC; i++){
            for(int j=0; j<N;j++){
                A1[i*N + j]=A[i*N+j];
                B1[i*N + j]=B[i*N+j];
                C1[i*N + j]=C[i*N+j];
            }
        }

        int cpROWS=N_LOC; 
        // I need  cpROWS I am starting to receive the matrix from the other processes starting from the blocks of the process 1.
        //the 0 block is already present
        for(int x = 1; x<Npes;x++)
        {   
            int N_LOC= (x<rest ?  N/Npes + 1 : N/Npes);
            

            MPI_Recv(&A1[cpROWS*N], N*N_LOC, MPI_DOUBLE, x, x, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&B1[cpROWS*N], N*N_LOC, MPI_DOUBLE, x, x+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&C1[cpROWS*N], N*N_LOC, MPI_DOUBLE, x, x+2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //x matrix A, x+1 matrixB and ect
            cpROWS+=N_LOC;
        }
    }else
    {
        MPI_Send(A, N*N_LOC, MPI_DOUBLE, 0,rank, MPI_COMM_WORLD);
        MPI_Send(B, N*N_LOC, MPI_DOUBLE, 0,rank+1, MPI_COMM_WORLD);
        MPI_Send(C, N*N_LOC, MPI_DOUBLE, 0,rank+2, MPI_COMM_WORLD);
    }
    
    

    if(!rank)
    {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    N, N, N, 1, A1, N, B1, N, 1, D1, N);
        int i=0;
        for(i=0; i<N*N;i++)
        {
            if(fabs(C1[i]-D1[i])>1e-8) 
            //the tolerance has been introduced to take into account a possible difference due to very far digits from 0
            {
                printf("Matrices are different.\n");
                printf("The values are supposed to be %f and %f\n",D1[i], C1[i]);
                break;
            }   
        }
        if(i==N*N)
        {
            printf("Matrices are equal.. HURRAY..\n");
        }
    }


    free(A1); free(B1); free(C1); free(D1);

    /////////////////////
    
    free(A);
    free(B);
    free(C);
    free(recvcounts);
    free(displacement);

    MPI_Finalize();
    return 0;
}