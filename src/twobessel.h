typedef struct config {
	double l1, l2;
	double nu1, nu2;
	double c_window_width;
	long Nk_sample;
} config;
void two_sph_bessel(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double *r1, double *r2, double **result);
void two_sph_bessel_binave(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double smooth_dlnr, int dimension, double *r1, double *r2, double **result);
void two_Bessel_binave(double *k1, double *k2, double **fk1k2, long N1, long N2, config *config, double smooth_dlnr, int dimension, double *r1, double *r2, double **result);
