# cython: boundscheck=False, wraparound=False, initializedcheck=False
"""
Port of gamma functions for calculating Fisher scores.
"""
'''

# CONSTANTS TAKEN FROM CEPHES under UNK setting.
cdef int sgngam = 0
DEF sizedouble = 8
DEF LOGPI = 1.14472988584940017414
DEF MAXLGM = 2.035093e36
DEF LS2PI = 0.91893853320467274178
DEF MAXLOG = 8.8029691931113054295988E1
DEF MACHEP = 1.11022302462515654042E-16
DEF big = 4.503599627370496e15
DEF biginv = 2.22044604925031308085e-16
# Multiply a phred score by LOG10E_X5_INV to get CHI2


cdef double_t[5] A = [8.11614167470508450300E-4,
                      -5.95061904284301438324E-4,
                      7.93650340457716943945E-4,
                      -2.77777777730099687205E-3,
                      8.33333333333331927722E-2]
cdef double_t[6] B = [-1.37825152569120859100E3,
                      -3.88016315134637840924E4,
                      -3.31612992738871184744E5,
                      -1.16237097492762307383E6,
                      -1.72173700820839662146E6,
                      -8.53555664245765465627E5]
cdef double_t[6] C = [-3.51815701436523470549E2, -1.70642106651881159223E4,
                      -2.20528590553854454839E5, -1.13933444367982507207E6,
                      -2.53252307177582951285E6, -2.01889141433532773231E6]

cdef inline double_t polevl(double x, double[] coef, int N) nogil:
    cdef double_t ans
    cdef int i
    cdef double_t * exps = <double_t *> malloc(sizedouble * (N + 1))
    ans = coef[N]
    exps[N] = 1.
    for i from N - 1 >= i > 0:
        exps[i] = x * exps[i + 1]
        ans += coef[i] * exps[i]
    free(exps)
    return ans


cdef inline double_t p1evl(double_t x, double_t[] coef, int N) nogil:
    cdef double_t ans
    cdef double_t * exps = <double_t *> malloc(sizedouble * (N + 1))
    cdef int i
    ans = x + coef[N]
    exps[N - 1] = x

    for i from N - 2 >= i > 0:
        exps[i] = x * exps[i + 1]
        ans += coef[i] * exps[i]
    free(exps)
    return ans


cdef inline double_t lgam(double_t x) nogil:
    cdef double_t p, q, u, w, z
    cdef int i

    sgngam = 1

    if x < -34.0:
        q = -x
        w = lgam(q)  # /* note this modifies sgngam! */
        p = floor(q)
        if(p == q):
            return INFINITY
        z = q - p
        if z > 0.5:
            p += 1.0
            z = p - q
        z = q * sin(PI * z)
        return LOGPI - log(z) - w if(z != 0.0) else INFINITY

    elif x < 13.0:
        z = 1.0
        p = 0.0
        u = x
        while u >= 3.0:
            p -= 1.0
            u = x + p
            z *= u
        while u < 2.0:
            if u == 0.0:
                return INFINITY
            z /= u
            p += 1.0
            u = x + p
        if z < 0.0:
            sgngam = -1
            z = -z
        else:
            sgngam = 1
        if u == 2.0:
            return log(z)
        p -= 2.0
        x = x + p
        p = x * polevl(x, B, 5) / p1evl(x, C, 6)
        return log(z) + p
    elif x > MAXLGM:
        return INFINITY

    q = (x - 0.5) * log(x) - x + LS2PI
    if x > 1.0e8:
        return q

    p = 1.0/ (x * x)
    if x >= 1000.0:
        q += ((7.9365079365079365079365e-4 * p -
               2.7777777777777777777778e-3) * p +
              0.0833333333333333333333) / x
    else:
        q += polevl(p, A, 4) / x
    return q

cdef inline double_t igam(double_t a, double_t x) nogil:
    cdef double_t ans, ax, c, r

    if(x <= 0 or a <= 0):
        return 0.0

    if(x > 1.0 and x > a):
        return 1.0 - igamc(a,x)

    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:
        return 0.0
    ax = exp(ax)

    # /* power series */
    r = a
    c = 1.0
    ans = 1.0

    while c / ans > MACHEP:
        r += 1.0
        c *= x / r
        ans += c

    return ans * ax / a


cdef inline double_t igamc(double_t a, double_t x) nogil:
    cdef double_t ans, ax, c, yc, r, t, y, z
    cdef double_t pk, pkm1, pkm2, qk, qkm1, qkm2
    if x <= 0 or a <= 0:
        return 1.0
    elif x < 1.0 or x < a:
        return 1.0 - igam(a,x)

    ax = a * log(x) - x - lgam(a)
    if ax < -MAXLOG:
        return 0.0
    ax = exp(ax)

    \''' continued fraction \'''
    y = 1.0 - a
    z = x + y + 1.0
    c = 0.0
    pkm2 = 1.0
    qkm2 = x
    pkm1 = x + 1.0
    qkm1 = z * x
    ans = pkm1/qkm1

    while t > MACHEP:
        c += 1.0
        y += 1.0
        z += 2.0
        yc = y * c
        pk = pkm1 * z - pkm2 * yc
        qk = qkm1 * z - qkm2 * yc
        if qk != 0:
            r = pk/qk
            t = fabs((ans - r) / r)
            ans = r
        else:
            t = 1.0
        pkm2 = pkm1
        pkm1 = pk
        qkm2 = qkm1
        qkm1 = qk
        if fabs(pk) > big:
            pkm2 *= biginv
            pkm1 *= biginv
            qkm2 *= biginv
            qkm1 *= biginv
    return ans * ax

'''

DEF LOG10E_X5_INV = 0.46051701859880917
# Multiply a phred score by this to convert

cdef inline double_t CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return phredInt * LOG10E_X5_INV


cdef inline double_t INV_CHI2_FROM_PHRED(int32_t phredInt) nogil:
    return -2 * c_log(1 - c_pow(10, (phredInt * -1.)))
