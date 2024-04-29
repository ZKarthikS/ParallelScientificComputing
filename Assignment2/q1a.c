/*
Part a of Q1:

Serial Program using first order upwind and quick schemes
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


// Defining Constants
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
    double domain_left = 0, domain_right = 2;
    double x, t=0, t_save=0;


    // Defining and initializing variables for first-order-upwind and QUICK schemes:
    double u_fou[nx] = {0}, u_fou_new[nx] = {0};
    double u_quick[nx] = {0}, u_quick_new[nx] = {0};

    double u_fou_save[nt_save][nx] = {0};
    double u_quick_save[nt_save][nx] = {0};

    for(i=0,x=domain_left; x<=domain_right; x+=dx,i++)
    {
        if(x>=0.5)
            break;
        u_fou[i] = sin(4*PI*x);
        u_quick[i] = u_fou[i];
    }



    // Computing the equation
    for(t=0,i=0; t<=t_max; t+=dt,i++)
    {
        // t_save to save the data at a time-step of 0.5s
        if(fabs(t-t_save) < dt/2)
        {
            int k = (int) (t_save/dt_save);
            
            memcpy(u_fou_save[k], u_fou, nx*sizeof(double));
            memcpy(u_quick_save[k], u_quick, nx*sizeof(double));

            t_save += dt_save;
        }


        for(j=1,x=domain_left+dx; x<domain_right; x+=dx,j++)
        {
            if(x==domain_left+dx)
            {
                u_quick_new[j] = u_quick[j] - c*dt*(u_quick[j] - u_quick[j-1])/dx;
                continue;
            }
            
            u_fou_new[j] = u_fou[j] - c*dt*(u_fou[j] - u_fou[j-1])/dx;

            u_quick_new[j] = u_quick[j] - c*dt*(3*u_quick[j] + u_quick[j-2] + 3*u_quick[j+1] - 7*u_quick[j-1])/(8*dx);
        }

        memcpy(u_quick, u_quick_new, nx*sizeof(double));
        memcpy(u_fou, u_fou_new, nx*sizeof(double));

    }


    // Saving and printing the numerical solution for the saved values
    for(i=0,t_save=0; t_save<=t_max; t_save+=dt_save, i++)
    {
        printf("x,   u_actual, u_fou, u_quick, t=%f\n",t_save);
        // for(j=0,x=domain_left; fabs(x-domain_right) < dx/2; x+=(dx),j+=1)
        for(j=0,x=domain_left; x<domain_right+dx/2; x+=(dx),j+=1)
        {
            printf("%+.4f, %+.4f, %+.4f, %+.4f\n", x,
                    ((x - c*t_save) <= 0.5 && (x - c*t_save) >=0 ) ? sin(4*PI*(x-c*t_save)) : 0,
                    u_fou_save[i][j], u_quick_save[i][j]);
        }
    }


    return 0;
}