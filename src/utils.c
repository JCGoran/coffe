#include <stdlib.h>
#include <math.h>

#include "utils.h"

#ifdef __cplusplus
#include "compatibility.h"
#endif

void mk_diag_g_to_ng(double *in, long N, double dlnk, double **out) {
	long i, j;
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++) {
			out[i][j] = 0.;
		}
		out[i][i] = in[i] / dlnk;
	}
}

void extrap_log_linear(double *fk, int N_origin, int N_extra, double *large_fk) {
	double dln_left, dln_right;
	int i;

	dln_left = log(fk[1]/fk[0]);
	// printf("fk[0],fk[1]: %.15e,%.15e,%.15e,%.15e,%.15e\n", fk[0],fk[1],fk[2],fk[3],fk[4]);
	if(fk[0]<=0.) {
		for(i=0; i<N_extra; i++) {
			large_fk[i] = 0.;
		}
	}
	else{
		for(i=0; i<N_extra; i++) {
			large_fk[i] = exp(log(fk[0]) + (i - N_extra) * dln_left);
		}
	}

	for(i=N_extra; i< N_extra+N_origin; i++) {
		large_fk[i] = fk[i - N_extra];
	}

	dln_right = log(fk[N_origin-1]/fk[N_origin-2]);
	if(fk[N_origin-1]<=0.) {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = 0.;
		}
	}
	else {
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i] = exp(log(fk[N_origin-1]) + (i - N_extra - N_origin +1) * dln_right);
		}
	}
}

void extrap_log_bilinear(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	double dln_left, dln_right;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		dln_left = log(fk[i-N_extra][1]/fk[i-N_extra][0]);
		dln_right = log(fk[i-N_extra][N_origin-1]/fk[i-N_extra][N_origin-2]);
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = exp(log(fk[i-N_extra][0]) + (j - N_extra) * dln_left);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = exp(log(fk[i-N_extra][N_origin-1]) + (j - N_extra - N_origin +1) * dln_right);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
	double dln_up, dln_down;
	for(j=0; j<N_origin+2*N_extra; j++) {
		dln_up = log(large_fk[N_extra+1][j]/large_fk[N_extra][j]);
		dln_down = log(large_fk[N_extra+N_origin-1][j]/large_fk[N_extra+N_origin-2][j]);
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = exp(log(large_fk[N_extra][j]) + (i - N_extra) * dln_up);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = exp(log(large_fk[N_extra+N_origin-1][j]) + (i - N_extra - N_origin +1) * dln_down);
			if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}


void extrap_bilinear(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	double dleft, dright;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		dleft = fk[i-N_extra][1]-fk[i-N_extra][0];
		dright = fk[i-N_extra][N_origin-1]-fk[i-N_extra][N_origin-2];
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = fk[i-N_extra][0] + (j - N_extra) * dleft;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][N_origin-1] + (j - N_extra - N_origin +1) * dright;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
	double dup, ddown;
	for(j=0; j<N_origin+2*N_extra; j++) {
		dup = large_fk[N_extra+1][j]-large_fk[N_extra][j];
		ddown = large_fk[N_extra+N_origin-1][j]-large_fk[N_extra+N_origin-2][j];
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = large_fk[N_extra][j] + (i - N_extra) * dup;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = large_fk[N_extra+N_origin-1][j] + (i - N_extra - N_origin +1) * ddown;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}

void extrap_2dzeros(double **fk, int N_origin, int N_extra, double **large_fk) {
	int i,j;
	for(i=N_extra; i<N_origin+N_extra; i++) {
		for(j=0; j<N_extra; j++) {
			large_fk[i][j] = 0.;
		}
		for(j=N_extra; j< N_extra+N_origin; j++) {
			large_fk[i][j] = fk[i-N_extra][j - N_extra];
			// if(isnan(large_fk[i][j])) {printf("large_fk nan at: %d,%d\n", i,j);large_fk[i][j]=0.;}
		}
		for(j=N_extra+N_origin; j< 2*N_extra+N_origin; j++) {
			large_fk[i][j] = 0.;
		}
	}
	for(j=0; j<N_origin+2*N_extra; j++) {
		for(i=0; i<N_extra; i++) {
			large_fk[i][j] = 0.;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
		for(i=N_extra+N_origin; i< 2*N_extra+N_origin; i++) {
			large_fk[i][j] = 0.;
			// if(isnan(large_fk[i][j])) {large_fk[i][j]=0.;}
		}
	}
}

// void resample_fourier_gauss(double *k, double *fk, config *config, double *k_sample, double *fk_sample) {
// 	long i;
// 	double dlnk = log(k[sizeof(k)-1]/k[0]) / (config->Nk_sample-1.);
// 	for(i=0; i<config->Nk_sample; i++) {
// 		k_sample[i] = k[0] * exp(i*dlnk);
// 		fk_sample[i] = 
// 	}
// }

complex_t gamma_lanczos(complex_t z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(creal(z) < 0.5) {return M_PI / (csin(M_PI*z)*gamma_lanczos(1. - z));}
	z -= 1;
	complex_t x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	complex_t t = z + 7.5;
	return sqrt(2*M_PI) * cpow(t, z+0.5) * cexp(-t) * x;
}


complex_t lngamma_lanczos(complex_t z) {
/* Lanczos coefficients for g = 7 */
	static double p[] = {
		0.99999999999980993227684700473478,
		676.520368121885098567009190444019,
		-1259.13921672240287047156078755283,
		771.3234287776530788486528258894,
		-176.61502916214059906584551354,
		12.507343278686904814458936853,
		-0.13857109526572011689554707,
		9.984369578019570859563e-6,
		1.50563273514931155834e-7};

	if(creal(z) < 0.5) {return log(M_PI) -clog(csin(M_PI*z)) - lngamma_lanczos(1. - z);}
	z -= 1;
	complex_t x = p[0];
	for(int n = 1; n < 9; n++){ x += p[n] / (z + (double)(n));}

	complex_t t = z + 7.5;
	return log(2*M_PI) /2.  + (z+0.5)*clog(t) -t + clog(x);
}

void g_l(double l, double nu, double *eta, complex_t *gl, long N) {
/* z = nu + I*eta
Calculate g_l = zln2 + lngamma( (l+nu)/2 + I*eta/2 ) - lngamma( (3+l-nu)/2 - I*eta/2 ) */
	long i;
	complex_t z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(z*log(2.) + lngamma_lanczos((l+z)/2.) - lngamma_lanczos((3.+l-z)/2.));
	}
}

void g_l_smooth(double l, double nu, double *eta, complex_t *gl, long N, double smooth_dlnr, double alpha_pow) {
/* z = nu + I*eta
Calculate g_l_smooth = g_l * exp((2.-z)*smooth_dlnr) / (2.-z) */
	long i;
	complex_t z;
	for(i=0; i<N; i++) {
		z = nu+I*eta[i];
		gl[i] = cexp(z*log(2.) + lngamma_lanczos((l+z)/2.) - lngamma_lanczos((3.+l-z)/2.));
		gl[i] *= (cexp((alpha_pow-z)*smooth_dlnr)-1.)/(alpha_pow-z);
	}
}

void c_window_2d(fftw_complex *out, double c_window_width, long halfN1, long halfN2) {
	// 'out' is N1*(halfN2+1) complex array
	long Ncut1, Ncut2;
	long N1;
	N1 = 2*halfN1;
	Ncut1 = (long)(halfN1 * c_window_width);
	Ncut2 = (long)(halfN2 * c_window_width);
	long i,j;
	double W;
	for(j=0; j<=Ncut2; j++) { // window for right-side
		W = (double)(j)/Ncut2 - 1./(2.*M_PI) * sin(2.*j*M_PI/Ncut2);
		for(i=0; i<N1; i++){
#ifdef __cplusplus
			out[i*(halfN2+1)+ halfN2-j][0] *= W;
			out[i*(halfN2+1)+ halfN2-j][1] *= W;
#else
			out[i*(halfN2+1)+ halfN2-j] *= W;
#endif
		}
	}
	for(i=0; i<=Ncut1; i++){ // window for center part (equivalent for up-down edges when rearanged)
		W = (double)(i)/Ncut1 - 1./(2.*M_PI) * sin(2.*i*M_PI/Ncut1);
		for(j=0; j<=halfN2; j++) {
#ifdef __cplusplus
			out[(halfN1-i)*(halfN2+1) + j][0] *= W;
			out[(halfN1-i)*(halfN2+1) + j][1] *= W;
			out[(halfN1+1+i)*(halfN2+1) + j][0] *= W;
			out[(halfN1+1+i)*(halfN2+1) + j][1] *= W;
#else
			out[(halfN1-i)*(halfN2+1) + j] *= W;
			out[(halfN1+1+i)*(halfN2+1) + j] *= W;
#endif
		}
	}
}

