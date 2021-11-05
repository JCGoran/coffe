/* boilerplate for CPP aliasing */
#ifndef COFFE_COMPATIBILITY
#define COFFE_COMPATIBILITY

static double creal(complex_t value)
{
    return std::real(value);
}

static double cimag(complex_t value)
{
    return std::imag(value);
}

static complex_t csin(complex_t value)
{
    return std::sin(value);
}

static complex_t ccos(complex_t value)
{
    return std::cos(value);
}

static complex_t cpow(complex_t value, complex_t exponent)
{
    return std::pow(value, exponent);
}

static complex_t cexp(complex_t value)
{
    return std::exp(value);
}

static complex_t clog(complex_t value)
{
    return std::log(value);
}

static complex_t csqrt(complex_t value)
{
    return std::sqrt(value);
}

#endif
