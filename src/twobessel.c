#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "utils.h"
#include "twobessel.h"

#ifdef __cplusplus
#include "compatibility.h"
#endif

#include <fftw3.h>

void two_sph_bessel(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double *r1, double *r2, double **result) {
	// printf("fk1k2[0][0]:%.15e\n", fk1k2[0][0]);
	if(N1 % 2) {printf("Please use even number of k1 !\n"); exit(0);}
	if(N2 % 2) {printf("Please use even number of k2 !\n"); exit(0);}
	long halfN1 = N1/2, halfN2 = N2/2;

	double k10, k20, r10,r20;
	k10 = k1[0];
	k20 = k2[0];

	double dlnk1, dlnk2;
	dlnk1 = log(k1[1]/k10);
	dlnk2 = log(k2[1]/k20);

	// Only calculate the m>=0 part
	double eta_m[halfN1+1], eta_n[halfN2+1];
	long i,j;
	for(i=0; i<=halfN1; i++) {eta_m[i] = 2*M_PI / dlnk1 / N1 * i;}
	for(j=0; j<=halfN2; j++) {eta_n[j] = 2*M_PI / dlnk2 / N2 * j;}

	complex_t g1[halfN1+1], g2[halfN2+1];
	g_l(config->l1, config->nu1, eta_m, g1, halfN1+1);
	g_l(config->l2, config->nu2, eta_n, g2, halfN2+1);

	// printf("g2[0]: %.15e+I*(%.15e)\n", creal(g2[0]),cimag(g2[0]));

	// calculate r1,r2 arrays
	for(i=0; i<N1; i++) {r1[i] = 1. / k1[N1-1-i];}
	for(i=0; i<N2; i++) {r2[i] = 1. / k2[N2-1-i];}
	r10 = r1[0];
	r20 = r2[0];

	// biased input func
	double *Pb;
	Pb = (double *)malloc(N1 * N2* sizeof(double));
	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			Pb[i*N2+j] = fk1k2[i][j] / pow(k1[i], config->nu1) / pow(k2[j], config->nu2) ;
		}
	}

	fftw_complex *out;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*(halfN2+1) );
	plan_forward = fftw_plan_dft_r2c_2d(N1, N2, Pb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);
	fftw_destroy_plan(plan_forward);
	free(Pb);

	c_window_2d(out, config->c_window_width, halfN1, halfN2);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	long ij;
	for(i=0; i<=halfN1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i]);
			out[ij][1] *= cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i]);
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i] ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	for(i=halfN1+1; i<N1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i]) * conj(g1[N1-i]));
			out[ij][1] *= cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i]) * conj(g1[N1-i]));
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i]) * conj(g1[N1-i]) ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));
	double *out_ifft;
	out_ifft = (double *)malloc(sizeof(double) * N1*N2 );
	plan_backward = fftw_plan_dft_c2r_2d(N1, N2, out, out_ifft, FFTW_ESTIMATE);

	fftw_execute(plan_backward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);

	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			result[i][j] = out_ifft[i*N2 + j] * M_PI / (16.*N1*N2 * pow(r2[j], config->nu2) * pow(r1[i], config->nu1));
		}
	}

	free(out_ifft);
}

void two_sph_bessel_binave(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double smooth_dlnr, int dimension, double *r1, double *r2, double **result) {
	// printf("fk1k2[0][0]:%.15e\n", fk1k2[0][0]);
	if(N1 % 2) {printf("Please use even number of k1 !\n"); exit(0);}
	if(N2 % 2) {printf("Please use even number of k2 !\n"); exit(0);}
	long halfN1 = N1/2, halfN2 = N2/2;

	double k10, k20, r10,r20;
	k10 = k1[0];
	k20 = k2[0];

	double dlnk1, dlnk2;
	dlnk1 = log(k1[1]/k10);
	dlnk2 = log(k2[1]/k20);

	// Only calculate the m>=0 part
	double eta_m[halfN1+1], eta_n[halfN2+1];
	long i,j;
	for(i=0; i<=halfN1; i++) {eta_m[i] = 2*M_PI / dlnk1 / N1 * i;}
	for(j=0; j<=halfN2; j++) {eta_n[j] = 2*M_PI / dlnk2 / N2 * j;}

	complex_t g1[halfN1+1], g2[halfN2+1];
	g_l_smooth(config->l1, config->nu1, eta_m, g1, halfN1+1, smooth_dlnr, dimension);
	g_l_smooth(config->l2, config->nu2, eta_n, g2, halfN2+1, smooth_dlnr, dimension);

	double s_d_lambda = (exp(dimension*smooth_dlnr) - 1.) / dimension;
	// printf("g2[0]: %.15e+I*(%.15e)\n", creal(g2[0]),cimag(g2[0]));

	// calculate r1,r2 arrays
	for(i=0; i<N1; i++) {r1[i] = 1. / k1[N1-1-i];}
	for(i=0; i<N2; i++) {r2[i] = 1. / k2[N2-1-i];}
	r10 = r1[0];
	r20 = r2[0];

	// biased input func
	double *Pb;
	Pb = (double *)malloc(N1 * N2* sizeof(double));
	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			Pb[i*N2+j] = fk1k2[i][j] / pow(k1[i], config->nu1) / pow(k2[j], config->nu2) ;
		}
	}

	fftw_complex *out;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*(halfN2+1) );
	plan_forward = fftw_plan_dft_r2c_2d(N1, N2, Pb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);
	fftw_destroy_plan(plan_forward);
	free(Pb);

	c_window_2d(out, config->c_window_width, halfN1, halfN2);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	long ij;
	for(i=0; i<=halfN1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i] );
			out[ij][1] *= cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i] );
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i] ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	for(i=halfN1+1; i<N1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * creal(g1[N1-i]) + cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * cimag(g1[N1-i]);
			out[ij][1] *= -creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * cimag(g1[N1-i]) + cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * creal(g1[N1-i]);
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i]) * conj(g1[N1-i]) ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));
	double *out_ifft;
	out_ifft = (double *)malloc(sizeof(double) * N1*N2 );
	plan_backward = fftw_plan_dft_c2r_2d(N1, N2, out, out_ifft, FFTW_ESTIMATE);

	fftw_execute(plan_backward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);

	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			result[i][j] = out_ifft[i*N2 + j] * M_PI / (16.*N1*N2 * pow(r2[j], config->nu2) * pow(r1[i], config->nu1)) / s_d_lambda/s_d_lambda;
		}
	}

	free(out_ifft);
}

void two_Bessel_binave(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double smooth_dlnr, int dimension, double *r1, double *r2, double **result) {
	// printf("fk1k2[0][0]:%.15e\n", fk1k2[0][0]);
	if(N1 % 2) {printf("Please use even number of k1 !\n"); exit(0);}
	if(N2 % 2) {printf("Please use even number of k2 !\n"); exit(0);}
	long halfN1 = N1/2, halfN2 = N2/2;

	double k10, k20, r10,r20;
	k10 = k1[0];
	k20 = k2[0];

	double dlnk1, dlnk2;
	dlnk1 = log(k1[1]/k10);
	dlnk2 = log(k2[1]/k20);

	// Only calculate the m>=0 part
	double eta_m[halfN1+1], eta_n[halfN2+1];
	long i,j;
	for(i=0; i<=halfN1; i++) {eta_m[i] = 2*M_PI / dlnk1 / N1 * i;}
	for(j=0; j<=halfN2; j++) {eta_n[j] = 2*M_PI / dlnk2 / N2 * j;}

	complex_t g1[halfN1+1], g2[halfN2+1];
	g_l_smooth(config->l1 -0.5, config->nu1, eta_m, g1, halfN1+1, smooth_dlnr, dimension+0.5);
	g_l_smooth(config->l2 -0.5, config->nu2, eta_n, g2, halfN2+1, smooth_dlnr, dimension+0.5);

	double s_d_lambda = (exp(dimension*smooth_dlnr) - 1.) / dimension;
	// printf("g2[0]: %.15e+I*(%.15e)\n", creal(g1[0]),cimag(g1[0]));

	// calculate r1,r2 arrays
	for(i=0; i<N1; i++) {r1[i] = 1. / k1[N1-1-i];}
	for(i=0; i<N2; i++) {r2[i] = 1. / k2[N2-1-i];}
	r10 = r1[0];
	r20 = r2[0];

	// biased input func
	double *Pb;
	Pb = (double *)malloc(N1 * N2* sizeof(double));
	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			Pb[i*N2+j] = fk1k2[i][j] / pow(k1[i], config->nu1 -0.5) / pow(k2[j], config->nu2 -0.5) ;
		}
	}

	fftw_complex *out;
	fftw_plan plan_forward, plan_backward;
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*(halfN2+1) );
	plan_forward = fftw_plan_dft_r2c_2d(N1, N2, Pb, out, FFTW_ESTIMATE);
	fftw_execute(plan_forward);
	fftw_destroy_plan(plan_forward);
	free(Pb);

	c_window_2d(out, config->c_window_width, halfN1, halfN2);
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));

	long ij;
	for(i=0; i<=halfN1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i]);
			out[ij][1] *= cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i]);
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, -I*eta_m[i]) * g1[i] ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	for(i=halfN1+1; i<N1; i++) {
		for(j=0; j<=halfN2; j++) {
			ij = i*(halfN2+1) + j;
#ifdef __cplusplus
			out[ij][0] *= creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * creal(g1[N1-i]) + cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * cimag(g1[N1-i]);
			out[ij][1] *= -creal(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * cimag(g1[N1-i]) + cimag(cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i])) * creal(g1[N1-i]);
			out[ij][0] = out[ij][0];
			out[ij][1] = -out[ij][1];
#else
			out[ij] *= cpow(k20*r20, -I*eta_n[j]) * g2[j] * cpow(k10*r10, I*eta_m[N1-i]) * conj(g1[N1-i]) ;
			out[ij] = conj(out[ij]);
#endif
		}
	}
	// printf("out[1]:%.15e+i*(%.15e)\n", creal(out[1]), cimag(out[1]));
	double *out_ifft;
	out_ifft = (double *)malloc(sizeof(double) * N1*N2 );
	plan_backward = fftw_plan_dft_c2r_2d(N1, N2, out, out_ifft, FFTW_ESTIMATE);

	fftw_execute(plan_backward);
	fftw_destroy_plan(plan_backward);
	fftw_free(out);

	for(i=0; i<N1; i++) {
		for(j=0; j<N2; j++) {
			result[i][j] = out_ifft[i*N2 + j] / (8.*N1*N2 * pow(r2[j], config->nu2 -0.5) * pow(r1[i], config->nu1 -0.5)) / s_d_lambda/s_d_lambda;
		}
	}

	free(out_ifft);
}
