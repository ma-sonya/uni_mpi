#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
int main(int argc, char** argv){
    int process_Rank, size_Of_Cluster;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size_Of_Cluster);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_Rank);
    struct timeval begin, end;
    int n;
    int** A;
    double* x;
    double* xnew;
    if(process_Rank==0){
        freopen(argv[1],"r",stdin);
        scanf("%d",&n);
        A=(int**)malloc((n)*sizeof(int*));
        x=(double*)malloc(n*sizeof(double));
        xnew=(double*)malloc(n*sizeof(double));
        for(int i=0;i<n;i++){
            A[i]=(int*)malloc((n+1)*sizeof(int));
        }
        gettimeofday(&begin, 0);
        for(int i=0;i<n;i++){
            for(int j=0;j<=n;j++)scanf("%d",A[i]+j);
        }
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //printf("%d\n",n);
    if(process_Rank!=0){

        x=(double*)malloc(n*sizeof(double));
        xnew=(double*)malloc(n*sizeof(double));
    A=(int**)malloc((n)*sizeof(int*));
    for(int i=0;i<n;i++){
            A[i]=(int*)malloc((n+1)*sizeof(int));
        }
    }
    for(int i=0;i<n;i++)
    MPI_Bcast(A[i],(n+1),MPI_INT,0,MPI_COMM_WORLD);
    for(int i=0;i<n;i++)x[i]=0;



        int step=n/size_Of_Cluster+1;
        int l=process_Rank*step;
        int r=(process_Rank+1)*(step);


    for(int t=1;t<=n;t++){
        MPI_Bcast(x,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(r>n)r=n;
        for(int k=l;k<r;k++){
            xnew[k]=A[k][n];
            for(int i=0;i<n;i++)if(i!=k)xnew[k]-=A[k][i]*x[i];
            xnew[k]/=A[k][k];
        }
        MPI_Gather(xnew+l, r-l, MPI_DOUBLE, x+l, r-l, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //for(int i=0;i<n;i++)x[i]=xnew[i];

    }
    if(process_Rank==0){
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;

    double elapsed = seconds + microseconds*1e-6;
   printf("Real time measured: %.3f seconds.\n", elapsed);
    for(int i=0;i<n;i++)printf("%f ",x[i]);
    printf("\n");
    freopen("check.txt","w",stdout);
    for(int i=0;i<n;i++){
        double ans=0;
        for(int j=0;j<n;j++)ans+=x[j]*A[i][j];
        ans-=A[i][n];
        printf("%f ",ans);
    }
    }

    MPI_Finalize();

}
