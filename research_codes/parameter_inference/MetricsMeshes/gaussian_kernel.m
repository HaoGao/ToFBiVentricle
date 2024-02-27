function K = gaussian_kernel(r2, s)

    K =  exp(-r2./ (s * s));