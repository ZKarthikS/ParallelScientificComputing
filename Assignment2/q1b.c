/*
Part b of Q1:

Parallel MPI Program using first order upwind and quick schemes
to solve 1D traveling wave equation
Using delta_x = 0.002
      delta_t = 0.0001

u(0,t) = 0; and u(L,t)=0

Domain x = 0 to 2
       t = 0 to 1
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<mpi.h>


// Defining constants
#define PI 3.14159265359 
#define dx 0.002
#define dt 0.0001
#define dt_save 0.5     // dt_save - time step for saving output
#define t_max 1.0       // t_max - time upto which computation is done
#define c 1.0
enum {nx = (int) (2/dx + 1), nt= (int) (t_max/dt + 1), nt_save= (int) (t_max/dt_save + 1)}; 



int main(int argc, char* argv[])
{
    int i, j;
    double domain_left = 0, domain_right = 2;   // x-Domain limits
    double x, t=0, t_save=0;

    //MPI Initialization:
    MPI_Status status;
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);



    //Domain Decomposition:
    int *displs = (int *)calloc(size,sizeof(int)); 
    int *scounts = (int *)calloc(size,sizeof(int));

    int rem = nx%size, qu = nx/size;
    int offset = 0;
    for(i=0;i<size;i++)
    {
        scounts[i] = qu;
        if(rem>0)
        {
            scounts[i] ++;
            rem--;
        }
        displs[i] = offset;
        offset += scounts[i];
    }

    int local_nx = scounts[rank];   // local_nx - number of elements for each proc
    // printf("Rank=%d, n=%d", rank, local_nx);


    // Defining and initializing variables for first-order-upwind and QUICK schemes:
    double *u_fou = (double *)calloc(local_nx, sizeof(double));
    double *u_fou_new = (double *)calloc(local_nx, sizeof(double));
    double *u_quick = (double *)calloc(local_nx, sizeof(double));
    double *u_quick_new = (double *)calloc(local_nx, sizeof(double));

    double **u_fou_save = (double **)calloc(nt_save, sizeof(double*));
    for(int i = 0; i < nt_save; i++) u_fou_save[i] = (double *)calloc(nx, sizeof(double));

    double **u_quick_save = (double **)calloc(nt_save, sizeof(double*));
    for(int i = 0; i < nt_save; i++) u_quick_save[i] = (double *)calloc(nx, sizeof(double));
    

    for(i=0,x=domain_left+displs[rank]*dx; i<local_nx; x+=dx,i++)
    {
        u_fou[i] = ((x - c*t) <= 0.5 && (x - c*t) >=0 ) ? sin(4*PI*(x-c*t)) : 0,
        u_quick[i] = u_fou[i];
    }

    double l_ext_fou=0; // Additional elements in the extremes for FOU scheme
    double l_ext_quick[2]={0}, r_ext_quick=0; // Additional elements in the extremes for quick scheme



    // Computing the equation
    for(t=0, i=0; t<=t_max; t+=dt, i++)
    {
        // t_save to save the data at a time-step of 0.5s
        if(fabs(t-t_save) < dt/2)
        {
            int k = (int) (t_save/dt_save);

            MPI_Gatherv(u_fou, local_nx, MPI_DOUBLE, u_fou_save[k], scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gatherv(u_quick, local_nx, MPI_DOUBLE, u_quick_save[k], scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            t_save += dt_save;            
        }

        
        // Gathering elements for adjacent procs for extremities
        if(rank!=0)
        {
            MPI_Recv(&l_ext_fou, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(&l_ext_quick, 2, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&u_quick[0], 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank!=size-1)
        {
            MPI_Send(&u_fou[local_nx-1], 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            MPI_Send(&u_quick[local_nx-2], 2, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            MPI_Recv(&r_ext_quick, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
        }

        
        // Left extreme elements computation
        j=0;
        if(rank==0)
        {
            // u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;
            u_fou_new[j] = 0;
            u_quick_new[j] = 0;
        }
        else
        {
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - l_ext_fou)/dx;
            u_quick_new[j] = u_quick[j] - c*dt*(3*u_quick[j] + l_ext_quick[0] + 3*u_quick[j+1] - 7*l_ext_quick[1])/(8*dx);
        }

        j=1;
        if(rank==0)
        {
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;
            u_quick_new[j] = u_fou_new[j];   
        }
        else
        {
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;
            u_quick_new[j] = u_quick[j] - c*dt*(3*u_quick[j] + l_ext_quick[1] + 3*u_quick[j+1] - 7*u_quick[j-1])/(8*dx);
        }


        
        // Non-extreme elements computation
        for(j=2, x=domain_left+displs[rank]*dx+dx; j<local_nx-1; x+=dx,j++)
        {
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;
            u_quick_new[j] = u_quick[j] - c*dt*(3*u_quick[j] + u_quick[j-2] + 3*u_quick[j+1] - 7*u_quick[j-1])/(8*dx);
        }


        // Right extreme element computation
        j=local_nx-1;
        if(rank==size-1)
        {
            u_fou_new[j] = 0;
            u_quick_new[j] = 0;
        }
        else
        {
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;
            u_quick_new[j] = u_quick[j] - c*dt*(3*u_quick[j] + u_quick[j-2] + 3*r_ext_quick - 7*u_quick[j-1])/(8*dx);   
        }

        memcpy(u_quick, u_quick_new, local_nx*sizeof(double));
        memcpy(u_fou, u_fou_new, local_nx*sizeof(double));

    }


    
    // Saving and printing the numerical solution for the saved values
    if(rank==0)
    {
        for(i=0,t_save=0; t_save<=t_max; t_save+=dt_save, i++)
        {
            printf("x,   u_actual, u_fou, u_quick, t=%f\n",t_save);
            for(j=0,x=domain_left; x<domain_right + dx/2; x+=(dx),j+=1)
            {
                printf("%+.4f, %+.4f, %+.4f, %+.4f\n", x,
                        ((x - c*t_save) <= 0.5 && (x - c*t_save) >=0 ) ? sin(4*PI*(x-c*t_save)) : 0,
                        u_fou_save[i][j], u_quick_save[i][j]);
            }
            printf("\n");
        }
    }

    // MPI_Gatherv(u_fou, local_nx, MPI_DOUBLE, u_fou_save[0], scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    free(displs);
    free(scounts);
    free(u_fou);
    free(u_fou_new);
    free(u_fou_save);
    free(u_quick);
    free(u_quick_new);
    free(u_quick_save);

    MPI_Finalize();
    return 0;
}