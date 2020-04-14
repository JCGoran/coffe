static int coffe_test_corrfunc(
    const struct coffe_corrfunc_t *cf
)
{
}

static int coffe_test_multipoles(
    const struct coffe_multipoles_t *mp
)
{
}

int main(void)
{
    struct coffe_parameters_t par;
    coffe_parse_default_parameters(&par);

    par.divergent = 1;
    par.nonzero_terms[8].n = 4, par.nonzero_terms[8].l = 0;

    #ifdef _OPENMP
    par.nthreads = omp_get_num_procs();
    #endif

    struct coffe_background_t bg;
    coffe_background_init(&par, &bg);
    coffe_test_background(&bg);

    struct coffe_integrals_t integrals[10];
    coffe_integrals_init(&par, &bg, integrals);
    coffe_test_integrals(integrals);

    struct coffe_multipoles_t mp;
    coffe_multipoles_init(&par, &bg, integrals, &mp);

    return 0;


    struct coffe_corrfunc_ang_t cf_ang;
    struct coffe_corrfunc_t cf;
    struct coffe_average_multipoles_t ramp;
    struct coffe_covariance_t cov_mp;
    struct coffe_covariance_t cov_ramp;
    struct coffe_corrfunc2d_t cf2d;

    coffe_output_init(
        &par, &bg,
        &cf_ang, &cf,
        &mp, &ramp,
        &cov_mp, &cov_ramp,
        &cf2d
    );

    return 0;
}
