/*
Part b of Q2:

Parallel Program using Jacobi iterative method
Using delta = 0.01
Convergence to 1e-4 using norm 

Domain x = -1 to +1
       y = -1 to +1

At x=+/-1 alone the boundary condition differs.
At all other boundaries -> phi = 0

Use row-wise or column-wise block decomposition to obtain numerical solution
*/

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<math.h>
#include<string.h>

#define PI 3.14159265359 
#define tol 0.0001   // Error tolerance 1E-4
#define del 0.01     // Define delta = 0.1
enum {n= (int) (2/del + 1)};    // Define n (21 when del=0.1)


// Defines the function q(x,y)
double q_(double x, double y)
{
    // if(x*x+y*y > 2)
    // {
    //     printf("Note: q returns value greater than 2 for x=%f, y=%f", x,y);
    // }
    return (x*x+y*y);
}


int main(int argc, char* argv[])
{
    int i, j;
    double domain_left = -1, domain_right = 1;
    double x,y;
    double proc_norm_err = 0;
    int iter = 0;
    double t1,t2;

    // Initialize MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    t1 = MPI_Wtime();   // Starting time

    double norm_err = 0;    // Error per processor


    // Define and Initialize variables per proc
    int *displs = (int *)malloc(size*sizeof(int)); 
    int *scounts = (int *)malloc(size*sizeof(int));

    int rem = n%size, qu = n/size;
    int offset = 0;
    for(i=0;i<size;i++)
    {
        scounts[i] = qu*n;
        if(rem>0)
        {
            scounts[i] += n;
            rem--;
        }
        displs[i] = offset;
        offset += scounts[i];
    }

    int local_rows = scounts[rank]/n;

    // double loc_phi[local_rows][n];
    double **loc_phi = (double **)calloc(local_rows, sizeof(double*));
    for(int i = 0; i < local_rows; i++) loc_phi[i] = (double *)calloc(n, sizeof(double));

    double loc_right[n]={0}, loc_left[n]={0};
    // double phi_new[local_rows][n] = {0};
    double **phi_new = (double **)calloc(local_rows, sizeof(double*));
    for(int i = 0; i < local_rows; i++) phi_new[i] = (double *)calloc(n, sizeof(double));

    double **q = (double **)calloc(local_rows, sizeof(double*));
    for(int i = 0; i < local_rows; i++) q[i] = (double *)calloc(n, sizeof(double));


    // Defining q(x,y) for all elements 
    for(i=0; i<local_rows;i++)
    {
        x = domain_left + (i+(int)displs[rank]/n)*del;
        for(j=0;j<n;j++)
        {
            y = domain_left + j*del;
            q[i][j] = q_(x,y);
        }
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    // MPI_Scatterv(phi, scounts, displs, MPI_DOUBLE, loc_phi[0], local_rows*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // if(rank==0)
    // {
    //     for(i=0;i<size;i++)
    //         // printf("%d ",displs[i]);
    //         printf("%d ", scounts[i]);
    // }


    // Defining the x=-1 boundary condition
    if(rank==0)
    {
        i=0;
        for(j=1;j<n-1;j++)
        {
            y = domain_left + j*del;
            phi_new[i][j] = sin(2*PI*y);
            loc_phi[i][j] = phi_new[i][j];
        }
    }


    // Performing the Jacobi Iterative Scheme
    while(1)
    {
        iter++;
        norm_err = 0.0;

        // Condition for interior points
        for(i=1;i<local_rows-1;i++)
        {
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (loc_phi[i][j+1] + loc_phi[i][j-1] + loc_phi[i-1][j] + loc_phi[i+1][j] + q[i][j]*del*del)/4;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }
        }

        // Condition for procs that aren't at the ends to compute extremeties
        if(rank!=0 && rank!=size-1)
        {
            i=0;
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (loc_phi[i][j+1] + loc_phi[i][j-1] + loc_left[j] + loc_phi[i+1][j] + q[i][j]*del*del)/4;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }

            i=local_rows-1;
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (loc_phi[i][j+1] + loc_phi[i][j-1] + loc_phi[i-1][j] + loc_right[j] + q[i][j]*del*del)/4;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }
        }

        else if(rank==0)
        {
            i=local_rows-1;
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (loc_phi[i][j+1] + loc_phi[i][j-1] + loc_phi[i-1][j] + loc_right[j] + q[i][j]*del*del)/4;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }
        }

        else if(rank==size-1)
        {
            i=0;
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (loc_phi[i][j+1] + loc_phi[i][j-1] + loc_left[j] + loc_phi[i+1][j] + q[i][j]*del*del)/4;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }

            i = local_rows-1;
            for(j=1;j<n-1;j++)
            {
                phi_new[i][j] = (4*phi_new[i-1][j] - phi_new[i-2][j])/3.0;
                norm_err += (phi_new[i][j]-loc_phi[i][j])*(phi_new[i][j]-loc_phi[i][j]);
            }
        }
        
        // Defining a Global Norm variable and collecting the norm from all procs
        double glob_norm = 0.0;

        MPI_Allreduce(&norm_err, &glob_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if(sqrt(glob_norm) <= tol)
            break;

        MPI_Barrier(MPI_COMM_WORLD);


        
        for(i=0;i<local_rows;i++)
            for(j=0;j<n;j++)
                loc_phi[i][j] = phi_new[i][j];
        // memcpy(loc_phi, phi_new, local_rows*n*sizeof(double));




        // Performing Send-Receives to ensure the end-points are transfered      
        if(rank==0)
        {
            MPI_Send(loc_phi[local_rows-1], n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            // printf("Sent right from rank=0");
            MPI_Recv(loc_right, n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
            // printf("Received right from rank=0");
        }
        else if(rank==size-1)
        {
            MPI_Recv(loc_left, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            // printf("Received left from rank=%d", rank);
            MPI_Send(loc_phi[0], n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            // printf("Sent left from rank=%d", rank);
        }
        else
        {
            MPI_Recv(loc_left, n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            // printf("Received left from rank=%d", rank);
            MPI_Send(loc_phi[local_rows-1], n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            // printf("Sent right from rank=%d", rank);
            MPI_Recv(loc_right, n, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
            // printf("Received right from rank=%d", rank);
            MPI_Send(loc_phi[0], n, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            // printf("Sent left from rank=%d", rank);
        }
    }

    double* loc_phi_flat = (double *)malloc(local_rows*n*sizeof(double));
    for (i = 0; i < local_rows; i++) 
    {
        for (j = 0; j < n; j++) 
        {
            loc_phi_flat[i * n + j] = loc_phi[i][j];
        }
    }

    // for(int k=0; k<size; k++)
    //     if(rank==k)
    //         for(i=0;i<local_rows;i++)
    //         {
    //             for(j=0;j<n;j++)
    //                 printf("%f\t",loc_phi[i][j]);
    //             printf("\n");
    //         }
    // printf("\n");

    double phi[n][n];

    MPI_Gatherv(loc_phi_flat, local_rows*n, MPI_DOUBLE, phi[0], scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    
    t2 = MPI_Wtime();   // End_time

    if(rank==0)
    {
        printf("Iterations: %d\n", iter);
        printf("Time taken: %f\n", t2-t1);
        // printf("Rank=%d\n", rank);

        // for(j=0;j<n;j++)
        // {
        //     for(int k=0; k<n; k++)
        //         printf("%f\t", phi[j][k]);
        //     printf("\n");
        // }
        // printf("\n");

        printf("x,\tphi(x,0)\n");
        for(i=0;i<n;i++)
            printf("%+.4f\t%+.4f\n",domain_left + i*del, phi[i][100]);

        printf("\n");

        printf("y,\tphi(0,y)\n");
        for(i=0;i<n;i++)
            printf("%+.4f\t%+.4f\n",domain_left + i*del, phi[100][i]);
    }

    free(displs);
    free(scounts);
    free(loc_phi);
    free(phi_new);
    free(q);
    free(loc_phi_flat);

    MPI_Finalize();

    return 0;
}