/*
 ============================================================================
 Name        : Shooting.c
 Author      : Ivan Syzonenko
 Version     :
 Copyright   : Free to copy for everyone
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define _X_YES_

void allocate(const double a, const double b, const double *h, double **x, double **u, double **u_prime, unsigned int *totiter);

void rungekutta(const double a, const double b, const double h, const double alpha0, const double alpha1, const double alpha_r_old, 
			/*OUT ->*/	double *x, double *u, double *u_prime, const unsigned int totiter);

// #ifdef _X_YES_
	double p2 (const double x);
	double p1 (const double x);
	double p0 (const double x);
// #else
// 	double p2 ();
// 	double p1 ();
// 	double p0 ();
// #endif
	double RHS_fx_func(const double x);
	double f1(const double x, const double u, const double v);
	double f2(const double x, const double u, const double v);






int main(void)
{
	clock_t begin, end;
	double 	a = 0, // lower integration limit
			b = 1, // upper integration limit
			h = 0.0000001, // runge-kutta step
	// below are our boundary coefficients
			alpha0 = 1, // for lower limit near first derivative
			alpha1 = -1, // for lower limit near zero derivative
			beta0 = 0, // for upper limit near first derivative
			beta1 = 1; // for upper limit near zero derivative

	unsigned int totiter = 0;
	double *x, *u, *u_prime;
	allocate(a, b, &h, &x, &u, &u_prime, &totiter);


	double	temp_error = 0,
			temp_error_old = 0;

	double 	alpha_r,
			alpha_r_old = 1;

	unsigned int stopper = 1;
		double tmp = 0;//really dummy variable

	begin = clock();
	rungekutta(a, b, h, alpha0, alpha1, alpha_r_old, x, u, u_prime, totiter);
	temp_error = u[totiter-1]*beta1 + u_prime[totiter-1]*beta0;
	alpha_r  = 42;
	

	while ( abs(temp_error) > 1E-8 && stopper < 100 )//usually should finish in several passes(<5). 
	{
		temp_error_old = temp_error;
		rungekutta(a, b, h, alpha0, alpha1, alpha_r, x, u, u_prime, totiter);
		temp_error = u[totiter-1]*beta1 + u_prime[totiter-1]*beta0;
		tmp = alpha_r;
		alpha_r = alpha_r - temp_error * (alpha_r - alpha_r_old)/(temp_error - temp_error_old);
		alpha_r_old = tmp;
		++stopper;
	}
	end = clock();
	printf("Total time %f\n", (double)(end - begin) / CLOCKS_PER_SEC);

	printf("\nTotal RK iterations = %d\n", stopper);

	// DONE WITH COMPUTATION. GOOD TIME TO CHECK THE RESULT

	printf("Now checking result against function: %s\n\
			recommended input is :\n\
			p2 = %s\n\
			p1 = %s\n\
			p0 = %s\n\
			RHS_fx_func = %s\n",
			"0.25 * (2 * x[i] * x[i] - x[i] - 1)", 
			"1", "4", "8", "6.75 + (3.75 + 0.5 * x) * x");
	double *check = (double*)malloc(totiter * sizeof(double));
	for (int i = 0; i < totiter; ++i)
	{
		check[i] = 0.25 * (2 * x[i] * x[i] - x[i] - 1);
	}

	double total_error = 0;

	for (int i = 0; i < totiter; ++i)
	{
		total_error += fabs(check[i] - u[i]);
//		printf("%5.10f - %5.10f = %5.10f\n", check[i], u[i], check[i] - u[i]);
	}

	printf("Total error is : %5.5f\n", total_error);

	printf("Done with all calculations. Exiting...\n");

	return EXIT_SUCCESS;
}


void rungekutta(const double a, const double b, const double h, const double alpha0, const double alpha1, const double alpha_r,
			/*OUT ->*/	 double *x, double *u,  double *u_prime, const unsigned int totiter)
{
	
	double *v = u_prime;
	x[0] = a;

	if (alpha1 == 0)
	{
		u[0] = alpha_r;
		v[0] = 0;
	}
	else
	{
		u[0] = -(alpha0/alpha1)*alpha_r;
		v[0] = alpha_r;
	}

	double 	k[4], l[4];

	for(int i=0; i < totiter - 1; ++i)
    {
		k[0] = h * f1(x[i], u[i] , v[i]);
		l[0] = h * f2(x[i], u[i] , v[i]);

		k[1] = h * f1(x[i] + h/2.0, u[i] + k[0]/2.0 , v[i] + l[0]/2.0);
		l[1] = h * f2(x[i] + h/2.0, u[i] + k[0]/2.0 , v[i] + l[0]/2.0);

		k[2] = h * f1(x[i] + h/2.0, u[i] + k[1]/2.0 , v[i] + l[1]/2.0);
		l[2] = h * f2(x[i] + h/2.0, u[i] + k[1]/2.0 , v[i] + l[1]/2.0);

		k[3] = h * f1(x[i] + h, u[i] + k[2], v[i] + l[2]);
		l[3] = h * f2(x[i] + h, u[i] + k[2], v[i] + l[2]);

		u[i+1] = u[i] + (k[0] + 2*k[1] + 2*k[2] + k[3])/6.0;
		v[i+1] = v[i] + (l[0]+ 2*l[1] + 2*l[2] + l[3])/6.0;
    }
}


void allocate(const double a, const double b, const double *h, double **x, double **u, double **u_prime, unsigned int *totiter)
{
	if (*h != 0)
	{
		*totiter = (unsigned int)round(abs(b - a)/ *h);
	}
	else
	{
		printf("Attempted to divide by zero\nExiting...\n\n");
		exit(EXIT_FAILURE);
	}

	*x		 = (double*)malloc( *totiter * sizeof(double));
	*u		 = (double*)malloc( *totiter * sizeof(double));
	*u_prime = (double*)malloc( *totiter * sizeof(double));

	for (int i = 0; i < *totiter-1; ++i)
		(*x)[i+1] = (*x)[i] + *h;
}

double f1(const double x, const double u, const double v)
{
	return v;
}


double f2(const double x, const double u, const double v)
{
    return ( RHS_fx_func(x) - p1(x) * v - p2(x)*u )/p0(x);
}

// #ifdef _X_YES_
double p2 (const double x)
// #else
// double p2 ()
// #endif
{
	return 1;
}

// #ifdef _X_YES_
double p1 (const double x)
// #else
// double p1 ()
// #endif
{
	return 4;
}

// #ifdef _X_YES_
double p0 (const double x)
// #else
// double p0 ()
// #endif
{
	return 8;
} 

// #ifdef _X_YES_
double RHS_fx_func(const double x)
// #else
// double RHS_fx_func()
// #endif
{
	return 6.75 + (3.75 + 0.5 * x) * x;
}
