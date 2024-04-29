/*
del^2 (phi) = -q where q = 2(2-x^2-y^2)
phi(+/-1, y) = 0; phi (x,+/-1) = 0

Part b:
Develop OpenMP program for Gauss-Seidel method using 
(a) Diagonal Approach

Given phi = (x^2-1)(y^2-1)
Using delta = 0.1 to validate the output
Considering allowed percentage error = 1%
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#define del 0.1   // Define delta = 0.1
enum {n= (int) (2/del + 1)};    // Define n (21 when del=0.1). Note:2 = domain size

// Defines the function q(x,y)
double q_(double x, double y)
{
    return (2*(2-x*x-y*y));
}


// Defines the expected value of phi(x,y)
double phi_expected(double x, double y)
{
    return ((x*x-1)*(y*y-1));
}



int main(int argc, char* argv[])
{
    double domain_left = -1, domain_right = 1;
    double q[n][n] = {0};
    double phi[n][n] = {0}, phi_exp[n][n] = {0};
    long int i,j;
    float x,y;
    int thread_count;


    // Getting thread count
    if(argc==2)
        thread_count = strtol(argv[1], NULL, 10);
    else
    {
        printf("A command line argument other than the name of the executable is required. Exiting the program.\n");
        return 1;
    }

    double start = omp_get_wtime();


    // Defining q(x,y) and expected phi(x,y) for all points in the domain
    #pragma omp parallel for /**/ default(none) shared(phi_exp, q, domain_left) private(i,j,x,y) num_threads(thread_count)
    for(i=0; i<n; i++)
    {
        x = domain_left + i*del;
        for(j=0; j<n; j++)
        {
            y = domain_left + j*del;
            q[i][j] = q_(x,y);
            phi_exp[i][j] = phi_expected(x,y);
        }
    }

    // for(i=0;i<n;i++)
    // {
    //     for(j=0;j<n;j++)
    //         printf("%f  ",phi[i][j]);
    //     printf("\n");
    // }

    // Flag is used to exit the loop based on percentage error
    // Iter is used to count the number of iterations
    long int flag=1, iter=0, nd=2*n-1, l=0;
    double allowed_err = 0.01, err;
    long int istart, iend;

    while(flag)
    {
        iter++;
        flag=0;

        for(l=0;l<nd;l++)
        {
            if(l<n)
            {
                istart = 0;
                iend = l;
            }
            else
            {
                istart = l-n;
                iend = n-1;
            }

            #pragma omp parallel for /**/ default(none) shared(phi, phi_exp, q, flag, l, istart, iend, allowed_err) private(i,j,err) num_threads(thread_count)
            for(i=istart; i<=iend; i++)
            {
                j=l-i;

                if(i>=1 && i<n-1 && j>=1 && j<n-1)
                    phi[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] + q[i][j]*del*del)/4;
                else
                    phi[i][j] = 0;

                if(fabs(phi_exp[i][j])>0.001)
                    err = (phi[i][j]-phi_exp[i][j])/phi_exp[i][j];
                else
                    err = phi[i][j];

                if(fabs(err) > allowed_err) flag++;

            }
        }
    }

    double end = omp_get_wtime();

    printf("For p=%d and del=%f, time taken = %f\n",thread_count, del, end-start);

    // printf("%d\n",iter);

    // for(i=0;i<n;i++)
    //     printf("phi[i][15] = %f, phi_exp[i][15] = %f\n", phi[i][15], phi_exp[i][15]);

    // for(i=0;i<n;i++)
    //     printf("%f  ", phi[i][15]);


}