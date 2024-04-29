/*
Part a of Q2:

Serial Program using Jacobi iterative method
Using delta = 0.1
Convergence to 1e-4 using norm 

Domain x = -1 to +1
       y = -1 to +1

At x=1 alone the boundary condition differs.
At all other boundaries -> phi = 0
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>


// Defining Constants
#define PI 3.14159265359 
#define tol 0.0001  // error tolerance
#define del 0.01     // Define delta = 0.1
enum {n= (int) (2/del + 1)};    // Define n (21 when del=0.1)


// Defines the function q(x,y)
double q_(double x, double y)
{
    return (x*x+y*y);
}


int main(int argc, char* argv[])
{
    int i, j;
    double domain_left = -1, domain_right = 1;
    double phi[n][n] = {0};
    double q[n][n] = {0};
    double x,y;
    double norm_err = 0;
    int iter = 0;
    clock_t start, end;
    double cpu_time_used;

    
    start = clock();


    // Defining the x=-1 boundary condition
    i = 0;
    x = domain_left + i*del;
    for(j=1;j<n-1;j++)
    {
        y = domain_left + j*del;
        phi[i][j] = sin(2*PI*y);
        // printf("%f\t", y);
    }
    printf("\n");

    
    // Defining q(x,y) for all elements
    for(int i=0; i<n; i++)
    {
        x = domain_left + i*del;
        for(int j=0; j<n; j++)
        {
            y = domain_left + j*del;
            q[i][j] = q_(x,y);
        }
    }

    // for(i=0;i<n;i++)
    // {
    //     printf("\n");
    //     for(j=0;j<n;j++)
    //         printf("%f\t", q[i][j]);
    // }

    double var=0.0;


    // Performing the GS Red-Black iterative scheme
    while(iter<=2000000)
    {
        iter++;
        if(iter%2)
            norm_err = 0.0;

        for(i=1;i<n-1;i++)
        {
            for(j=1;j<n-1;j++)
            {
                if((iter+i+j)%2)
                {
                    var = (phi[i][j+1] + phi[i][j-1] + phi[i-1][j] + phi[i+1][j] + q[i][j]*del*del)/4;
                    norm_err += (var - phi[i][j])*(var - phi[i][j]);
                    phi[i][j] = var;
                    // norm_err += fabs(var[i][j] - phi[i][j]);//phi[i][j];
                }
            }
        }

        i = n-1;
        for(j=1;j<n-1;j++)
        {
            if((iter+i+j)%2)
            {
                var = (4*phi[i-1][j] - phi[i-2][j])/3;
                norm_err += (var - phi[i][j])*(var - phi[i][j]);
                phi[i][j] = var;
                // norm_err += fabs(var - phi[i][j]);//phi[i][j];
            }
        }

        // for(i=0;i<n;i++)
        // {
        //     printf("\n");
        //     for(j=0;j<n;j++)
        //         printf("%f\t", phi[i][j]);
        // }
        // printf("\n");

        if(iter%2==0)
        {
            if(sqrt(norm_err) <= tol)
                break;
        }

        if(iter==20000000)
            printf("The scheme did not converge");
    }

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\n Time taken when del = %f is t=%f\n", del, cpu_time_used);
    printf("iterations = %d\n",iter/2);

    // for(i=0;i<n;i++)
    // {
    //     printf("\n");
    //     for(j=0;j<n;j++)
    //         printf("%f\t", phi[i][j]);
    // }


    // Printing the solution for x=0.0 and y=0.0:
    printf("x,\tphi(x,0)\n");
    for(i=0;i<n;i++)
        printf("%+.4f\t%+.4f\n",domain_left + i*del, phi[i][100]);

    printf("\n");

    printf("y,\tphi(0,y)\n");
    for(i=0;i<n;i++)
        printf("%+.4f\t%+.4f\n",domain_left + i*del, phi[100][i]);

    printf("\n");


    return 0;
}