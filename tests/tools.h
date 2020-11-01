#ifndef COFFE_TESTS_PRINT_SUCCESS
#define COFFE_TESTS_PRINT_SUCCESS \
    fprintf(stderr, "All tests in `%s` passed successfully\n", __func__);
#endif

/* our desired accuracy */
#ifndef COFFE_REL_EPSILON
#define COFFE_REL_EPSILON (5e-4)
#endif
#ifndef COFFE_ABS_EPSILON
#define COFFE_ABS_EPSILON (1e-10)
#endif

/* 1 if equal, else 0 */
static int approx_equal_const_epsilon(
    const double a,
    const double b
)
{
    const double rel_epsilon = COFFE_REL_EPSILON, abs_epsilon = COFFE_ABS_EPSILON;
    if (a != 0 && b != 0){
        if (fabs(a - b) / fabs(a) < rel_epsilon)
            return 1;
        else{
            fprintf(stderr, "ERROR: %e != %e to %e\n", a, b, rel_epsilon);
            return 0;
        }
    }
    else{
        if (fabs(a - b) < abs_epsilon)
            return 1;
        else{
            fprintf(stderr, "ERROR: %e != %e to %e\n", a, b, abs_epsilon);
            return 0;
        }
    }
}

/* 1 if equal, else 0 */
static int approx_equal(
    const double a,
    const double b,
    const double rel_epsilon,
    const double abs_epsilon
)
{
    if (a != 0 && b != 0){
        if (fabs(a - b) / fabs(a) < rel_epsilon)
            return 1;
        else{
            fprintf(stderr, "ERROR: %e != %e to %e\n", a, b, rel_epsilon);
            return 0;
        }
    }
    else{
        if (fabs(a - b) < abs_epsilon)
            return 1;
        else{
            fprintf(stderr, "ERROR: %e != %e to %e\n", a, b, abs_epsilon);
            return 0;
        }
    }
}
/* a version of assert that DOESN'T exit on error */
void weak_assert(
    const int value,
    int *error_flag
)
{
    if (!value){
        if (!(*error_flag))
            *error_flag = 1;
    }
}
