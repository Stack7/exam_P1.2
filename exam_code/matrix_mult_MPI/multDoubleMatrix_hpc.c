#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <sys/time.h>

double seconds(){

  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
  return sec;
}

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
    int N=0; 
    if (argc>1)
    {  
      N=atoi(argv[1]);
    }
    int rank, Npes; 
  
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Npes);

    // define time variables
    double time_comm = 0; 
    double time_comput =0;
    double time_comm_start, time_comm_end;
    double time_comput_start, time_comput_end;

    double mean_time_comm = 0;
    double mean_time_comput = 0;


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
    double* recvbuf = (double*)malloc( size_col[0]*N* sizeof(double));
    double* sendbuf=(double*)malloc(N_LOC *  size_col[0]* sizeof(double));
    
    for(int m=0;m<N_LOC;m++) // my process p_n knows inside all the row into it! so just traverse the row
    {
        for(int n=0;n<N;n++)
        {
            A[m*N+n]=n;//rand(rank)%10;
            B[m*N+n]=1.0/(n+1);//(rand(rank)%10);
            C[m*N+n]=0;
        }  
    }   

    int  j_glob= 0;   
    // my process p_n knows its rows 
    // I need that every processor takes from the ithers their blocks of the column
    for(int i=0; i<Npes;i++)    
    {
        time_comm_start = seconds();     
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

        //calculating the time dedicated to the communication between nodes  
        time_comm_end =seconds(); 
        time_comm += time_comm_end -time_comm_start;

        time_comput_start=seconds();
        //calculating the time dedicated to COMPUTE the matrix multiplication        
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
        time_comput_end = seconds();
        time_comput +=  time_comput_end- time_comput_start;

    }  //end for 


    // Checking the time of all the processors
    double total_time = time_comm +time_comput; 
    //total time for all the procs
    // double  mean_total_time=0; 
    //mean time per processor
    //each processor sends to the 0 proc its time of communication and computation; then with reduce it will be summed all together

    time_comm/=Npes;
    time_comput/=Npes;
    
    MPI_Reduce(&time_comm, &mean_time_comm, 1, MPI_DOUBLE, MPI_SUM, 0,  MPI_COMM_WORLD);
    MPI_Reduce(&time_comput, &mean_time_comput, 1, MPI_DOUBLE, MPI_SUM, 0,  MPI_COMM_WORLD);

    if(!rank)
    {
        double mean_comp_time = mean_time_comput;
        double mean_comm_time = mean_time_comm;
        double mean_total_time = mean_comm_time + mean_comp_time;
        FILE *fp = fopen("all_times.txt","a");
        fprintf(fp, "%d %d %lf %lf %lf\n",Npes, N*N, mean_total_time, mean_comm_time, mean_comp_time);
        fclose(fp);
    }

    free(sendbuf);
    free(recvbuf);
   
    free(A);
    free(B);
    free(C);
    free(recvcounts);
    free(displacement);

    MPI_Finalize();
    return 0;
}