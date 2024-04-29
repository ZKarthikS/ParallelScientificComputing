/*
Solve f(x) = sin(5x) in the domain 0<=x<=3
using four-order accurate Padé scheme for the interior and 
third-order accurate one-sided Padé scheme near the boundaries.

Part b:
OpenMP Program to compute the derivative using 
recursive-doubling algorithm.
Use n=100 for plot with p=2 (no. of grid points, compute h accordingly)
Note: Grid points range from 0 to n (both inclusive)

Compare time between p=2,4,8 for n=1000
*/

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#define n 100     // n defined as a macro

// Function Definition
float f(float x)
{
    return ( sin(5*x) );
}


int main(int argc, char* argv[])
{
    float domain_left=0.0, domain_right=3.0;
    float h=(domain_right-domain_left)/n;
    float a[n+1]= {0}, b[n+1]={0}, c[n]={0};
    float a_new[n+1]= {0}, b_new[n+1]={0}, c_new[n]={0}, y_new[n+1]={0};
    float y[n+1] = {0};
    float x[n+1] = {0};
    float der_f[n+1] = {0};
    float alpha[n+1] = {0}, beta[n+1] = {0};
    int i,j,k;
    int left, right;
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

    // Initialize x:
    #pragma omp parallel for default(none) shared(x, h, domain_left) private(i) num_threads(thread_count)
    for(i=0;i<=n;i++)
        x[i] = domain_left + i*h;


    // Condition near the boundary:
    b[0] = 1;
    b[n] = 1;
    c[0] = 2;
    a[n] = 2;

    y[0] = (2*f(x[1]) + 0.5*f(x[2]) - 2*f(x[0]))/h;
    y[n] = (2.5*f(x[n]) - 2*f(x[n-1]) - 0.5*f(x[n-2]))/h;


    // Interior Scheme:
    #pragma omp parallel for default(none) shared(a,b,c,y,x, h, domain_left) private(i) num_threads(thread_count)
    for(i=1;i<n;i++)
    {
        b[i] = 4;
        a[i] = 1;
        c[i] = 1;

        y[i] = 3*(f(x[i+1])-f(x[i-1]))/h;
    }
    // printf("\n Interior Scheme good\n");

    // Recursive Doubling:
    int N = ceil(log2(n+1));
    for(k=1;k<=N;k*=2)
    {
        #pragma omp parallel for /*collapse*/ default(none) shared(a, b, c, alpha, beta, a_new, b_new, c_new, y_new, y, k) private(i,left,right) num_threads(thread_count)
        for(i=0;i<=n;i++)
        {
            left = i-pow(2,k-1);
            right = i+pow(2,k-1);
            
            if(left<0 && right<n)
            {
                alpha[i] = 0;
                a_new[i] = 0;

                beta[i] = -c[i]/b[right];
                c_new[i] = beta[i]*c[right];

                b_new[i] = b[i] + beta[i]*a[right];
                y_new[i] = y[i] + beta[i]*y[right];
            }

            else if(left<0 && right>=n)
            {
                alpha[i] = 0;
                a_new[i] = 0;

                beta[i] = 0;
                c_new[i] = 0;

                b_new[i] = b[i];
                y_new[i] = y[i];
            }

            else if(left>0 && right>=n)
            {
                alpha[i] = -a[i]/b[left];
                a_new[i] = alpha[i]*a[left];

                beta[i] = 0;
                c_new[i] = 0;

                b_new[i] = b[i] + alpha[i]*c[left];
                y_new[i] = y[i] + alpha[i]*y[left];
            }

            else
            {
                alpha[i] = -a[i]/b[left];
                a_new[i] = alpha[i]*a[left];

                beta[i] = -c[i]/b[right];
                c_new[i] = beta[i]*c[right];

                b_new[i] = b[i] + alpha[i]*c[left] + beta[i]*a[right];
                y_new[i] = y[i] + alpha[i]*y[left] + beta[i]*y[right];
            }
        }

        #pragma omp parallel for /*collapse*/ default(none) shared(a, b, c, a_new, b_new, c_new, y_new, y) private(i) num_threads(thread_count)
        for(i=0;i<=n;i++)
        {
            a[i] = a_new[i];
            b[i] = b_new[i];
            if(i<n)
                c[i] = c_new[i];
            y[i] = y_new[i];
        }
    }

    #pragma omp parallel for /*collapse*/ default(none) shared(b, y, der_f) private(i) num_threads(thread_count)
    for(i=0;i<=n;i++)
    {
        der_f[i] = y[i]/b[i];
    }

    double end = omp_get_wtime();

    printf("%f is the time taken for p=%d threads\n", end-start, thread_count);

    // printf("\n");

    // for(i=0;i<=n;i++)
    //     printf("%f  %f\n",x[i], der_f[i]);

    // for(i=0;i<=n;i++)
    //     printf("%f ", der_f[i]);
    
    // printf("\n");
}