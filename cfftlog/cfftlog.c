#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <string.h>

#include <time.h>

#include <fftw3.h>

#include "utils.h"
#include "utils_complex.h"
#include "cfftlog.h"

void cfftlog(double *x, double *fx, long N, config *config, int ell, double *y, double *Fy) {

	if(N % 2) {printf("Please use even number of x !\n"); exit(0);}
	long halfN = N/2;

	double x0, y0;
	x0 = x[0];

	double dlnx;
	dlnx = log(x[1]/x0);

	// Only calculate the m>=0 part
	double eta_m[halfN+1];
	long i;
	for(i=0; i<=halfN; i++) {eta_m[i] = 2*M_PI / dlnx / N * i;}

	double complex gl[halfN+1];
	
	switch(config->derivative) {
		case 0: g_l((double)ell, config->nu, eta_m, gl, halfN+1); break;
		case 1: g_l_1((double)ell, config->nu, eta_m, gl, halfN+1); break;
		case 2: g_l_2((double)ell, config->nu, eta_m, gl, halfN+1); break;
		default: printf("Integral Not Supported! Please choose config->derivative from [0,1,2].\n");
	}
	// printf("g2[0]: %.15e+I*(%.15e)\n", creal(g2[0]),cimag(g2[0]));

	// calculate y arrays
	for(i=0; i<N; i++) {y[i] = (ell+1.) / x[N-1-i];}
	y0 = y[0];

	// biased input func
	double *fb;
	fb = malloc(N* sizeof(double));
	for(i=0; i<N; i++) {
		fb[i] = fx[i] / pow(x[i], config->nu) ;
	}

	fftw_complex *out;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	plan_forward = fftw_plan_dft_r2c_1d(N, fb, out, FFTW_ESTIMATE);

	fftw_execute(plan_forward);

	c_window(out, config->c_window_width, halfN);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	for(i=0; i<=halfN; i++) {
		out[i] *= cpow(x0*y0, -I*eta_m[i]) * gl[i] ;
		out[i] = conj(out[i]);
	}

	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));
	double *out_ifft;
	out_ifft = malloc(sizeof(double) * N );
	plan_backward = fftw_plan_dft_c2r_1d(N, out, out_ifft, FFTW_ESTIMATE);

	fftw_execute(plan_backward);

	for(i=0; i<N; i++) {
		Fy[i] = out_ifft[i] * sqrt(M_PI) / (4.*N * pow(y[i], config->nu));
	}

	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);
	free(out_ifft);
}

void cfftlog_ells(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy) {

	if(N % 2) {printf("Please use even number of x !\n"); exit(0);}
	long halfN = N/2;

	double x0, y0;
	x0 = x[0];

	double dlnx;
	dlnx = log(x[1]/x0);

	// Only calculate the m>=0 part
	double eta_m[halfN+1];
	long i, j;
	for(i=0; i<=halfN; i++) {eta_m[i] = 2*M_PI / dlnx / N * i;}

	double complex gl[halfN+1];
	
	// biased input func
	double *fb;
	fb = malloc(N* sizeof(double));
	for(i=0; i<N; i++) {
		fb[i] = fx[i] / pow(x[i], config->nu) ;
	}

	fftw_complex *out, *out_vary;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	out_vary = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (halfN+1) );
	plan_forward = fftw_plan_dft_r2c_1d(N, fb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);

	c_window(out, config->c_window_width, halfN);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	double *out_ifft;
	out_ifft = malloc(sizeof(double) * N );
	plan_backward = fftw_plan_dft_c2r_1d(N, out_vary, out_ifft, FFTW_ESTIMATE);

	for(j=0; j<Nell; j++){
		switch(config->derivative) {
			case 0: g_l((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 1: g_l_1((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			case 2: g_l_2((double)ell[j], config->nu, eta_m, gl, halfN+1); break;
			default: printf("Integral Not Supported! Please choose config->derivative from [0,1,2].\n");
		}

		// calculate y arrays
		for(i=0; i<N; i++) {y[j][i] = (ell[j]+1.) / x[N-1-i];}
		y0 = y[j][0];

		for(i=0; i<=halfN; i++) {
			out_vary[i] = conj(out[i] * cpow(x0*y0, -I*eta_m[i]) * gl[i]) ;
		}

		fftw_execute(plan_backward);

		for(i=0; i<N; i++) {
			Fy[j][i] = out_ifft[i] * sqrt(M_PI) / (4.*N * pow(y[j][i], config->nu));
		}
	}
	fftw_destroy_plan(plan_forward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);
	fftw_free(out_vary);
	free(out_ifft);
}
