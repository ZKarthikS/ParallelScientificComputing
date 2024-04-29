/*
Solve f(x) = sin(5x) in the domain 0<=x<=3
using four-order accurate Padé scheme for the interior and 
third-order accurate one-sided Padé scheme near the boundaries.

Part a:
Serial Program to compute the derivative using 
tridiagonal LU Decomposition.
Use n=25 (no. of grid points, compute h accordingly)
Note: Grid points range from 0 to n (both inclusive)
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#define n 25     // n=25 defined as a macro

// Function Definition
float f(float x)
{
    return ( sin(5*x) );
}


int main(int argc, char* argv[])
{
    float domain_left=0.0, domain_right=3.0;
    float h=(domain_right-domain_left)/n;
    float a[n+1][n+1] = {0};
    float b[n+1] = {0};
    float x[n+1] = {0};
    float der_f[n+1] = {0};
    float m[n+1][n+1] = {0};
    float d[n+1] = {0};
    int i,j,k;


    // Initialize x:
    for(i=0;i<=n;i++)
        x[i] = domain_left + i*h;

    // Condition near the boundary:

    a[0][0] = 1;
    a[n][n] = 1;
    a[0][1] = 2;
    a[n][n-1] = 2;

    b[0] = (2*f(x[1]) + 0.5*f(x[2]) - 2*f(x[0]))/h;
    b[n] = (2.5*f(x[n]) - 2*f(x[n-1]) - 0.5*f(x[n-2]))/h;


    // Interior Scheme:
    for(i=1;i<n;i++)
    {
        a[i][i] = 4;
        a[i][i-1] = 1;
        a[i][i+1] = 1;

        b[i] = 3*(f(x[i+1])-f(x[i-1]))/h;
    }

    // LU Decomposition:
    // m=L; a=U
    for(k=0;k<n+1;k++)
    {
        m[k][k]=1;
        for(i=k+1;i<n+1;i++)
        {
            m[i][k] = a[i][k]/a[k][k];
            for(j=0;j<n+1;j++)
            {
                if(j<k+1)
                    a[i][j]=0;
                else
                    a[i][j] -= m[i][k]*a[k][j];
            }
        }
    }

    // m[][]*d[] = b[]
    for(i=0;i<n+1;i++)
    {
        d[i] = b[i];
        for(j=0;j<i;j++)
        {
            d[i]-=m[i][j]*d[j];
        }
    }

    // a[][]x[] = d[]
    for(i=n;i>=0;i--)
    {
        der_f[i] = d[i];
        for(j=i+1;j<n+1;j++)
        {
            der_f[i] -= a[i][j]*der_f[j];
        }
        der_f[i] /= a[i][i];
    }

    printf("\n");

    for(i=0;i<=n;i++)
        printf("%f  %f\n",x[i], der_f[i]);
    
    // printf("\n");

    // for(i=0;i<=n;i++)
    //     printf("%f ", der_f[i]);

    // printf("\n");

    // for(i=0;i<=n;i++)
    //     printf("%f ",b[i]);

    // printf("\n");

    // for(i=0;i<=n;i++)
    //     printf("%f ",d[i]);


    // for(i=0;i<=n;i++)
    // {
    //     for(j=0;j<=n;j++)
    //         printf("%f ",m[i][j]);
    //     printf("\n");
    // }
}