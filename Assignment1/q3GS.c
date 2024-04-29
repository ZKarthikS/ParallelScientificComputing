/*
del^2 (phi) = -q where q = 2(2-x^2-y^2)
phi(+/-1, y) = 0; phi (x,+/-1) = 0

Part a:
Assume delta = 0.1
Develop serial Gauss-Seidel program to solve the equations
Report number of iterations taken to be within 1% of the solution

Given phi = (x^2-1)(y^2-1)
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#define del 0.1     // Define delta = 0.1
enum {n= (int) (2/del + 1)};    // Define n (21 when del=0.1)

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
    int i,j;
    float x,y;
    clock_t start, end;
    double cpu_time_used;


    
    start = clock();

    // Defining q(x,y) and expected phi(x,y) for all points in the domain
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

    // Flag is used to exit the loop based on percentage error
    // Iter is used to count the number of iterations
    int flag=1, iter=0;
    double allowed_err = 0.01;

    while(flag)
    {
        iter++;
        flag=0;
        
        for(i=1;i<n-1;i++)
        {
            for(j=1;j<n-1;j++)
            {
                phi[i][j] = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1] + phi[i][j-1] + q[i][j]*del*del)/4;

                double err = (phi[i][j]-phi_exp[i][j])/phi_exp[i][j];

                if(fabs(err) > allowed_err) flag++;
            }
        }
    }

    end = clock();

    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("\n Time taken when del = %f is t=%f\n", del, cpu_time_used);

    printf("%d\n",iter);

    // for(i=0;i<n;i++)
    //     printf("phi[i][15] = %f, phi_exp[i][15] = %f\n", phi[i][15], phi_exp[i][15]);

    // for(i=0;i<n;i++)
    //     printf("%f  ", phi[i][15]);

    // for(i=0;i<n;i++)
    //     printf("%f  ", phi_exp[i][15]);


}