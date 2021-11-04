#ifdef __cplusplus
#include <complex>
typedef std::complex<double> complex_t;
using namespace std::complex_literals;
#define I 1i
#else
#include <complex.h>
typedef double complex complex_t;
#endif

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

void g_l(double l, double nu, double *eta, complex_t *gl, long N);

void g_l_smooth(double l, double nu, double *eta, complex_t *gl, long N, double smooth_dlnr, double alpha_pow);

void c_window_2d(fftw_complex *out, double c_window_width, long halfN1, long halfN2);

// void resample_fourier_gauss(double *k, double *fk, config *config);

complex_t gamma_lanczos(complex_t z);
complex_t lngamma_lanczos(complex_t z);
