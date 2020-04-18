#include <complex.h>
#include <fftw3.h>

typedef struct preconfig {
	double pre_k1min, pre_k1max;
	double pre_k2min, pre_k2max;
} preconfig;

void mk_diag_g_to_ng(double *in, long N, double dlnk, double **out);

void extrap_log_linear(double *fk, int N_origin, int N_extra, double *large_fk);

void extrap_log_bilinear(double **fk, int N_origin, int N_extra, double **large_fk);

void extrap_bilinear(double **fk, int N_origin, int N_extra, double **large_fk);

void extrap_2dzeros(double **fk, int N_origin, int N_extra, double **large_fk);

void g_l(double l, double nu, double *eta, double complex *gl, long N);

void g_l_smooth(double l, double nu, double *eta, double complex *gl, long N, double smooth_dlnr, double alpha_pow);

void c_window_2d(double complex *out, double c_window_width, long halfN1, long halfN2);

// void resample_fourier_gauss(double *k, double *fk, config *config);

double complex gamma_lanczos(double complex z);
double complex lngamma_lanczos(double complex z);