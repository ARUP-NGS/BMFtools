/*                            igami()
 *
 *      Inverse of complemented imcomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igami();
 *
 * x = igami(a, y);
 *
 *
 *
 * DESCRIPTION:
 *
 * Given y, the function finds x such that
 *
 *  igamc(a, x) = y.
 *
 * Starting with the approximate value
 *
 *         3
 *  x = a t
 *
 *  where
 *
 *  t = 1 - d - ndtri(y) sqrt(d)
 *
 * and
 *
 *  d = 1/9a,
 *
 * the routine performs up to 10 Newton iterations to find the
 * root of igamc(a,x) - y = 0.
 *
 *
 * ACCURACY:
 *
 * Tested for a ranging from 0.5 to 30 and x from 0 to 0.5.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,0.5         3400       8.8e-16     1.3e-16
 *    IEEE      0,0.5        10000       1.1e-14     1.0e-15
 *
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

#include <tgmath.h>
#include "igamc_cephes.h"


#if 1
double MACHEP =  1.11022302462515654042E-16;   /* 2**-53 */
#else
double MACHEP =  1.38777878078144567553E-17;   /* 2**-56 */
#endif
double UFLOWTHRESH =  2.22507385850720138309E-308; /* 2**-1022 */
#ifdef DENORMAL
double MAXLOG =  7.09782712893383996732E2;     /* log(MAXNUM) */
/* double MINLOG = -7.44440071921381262314E2; */     /* log(2**-1074) */
double MINLOG = -7.451332191019412076235E2;     /* log(2**-1075) */
#else
double MAXLOG =  7.08396418532264106224E2;     /* log 2**1022 */
double MINLOG = -7.08396418532264106224E2;     /* log 2**-1022 */
#endif
double MAXNUM =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
double PI     =  3.14159265358979323846;       /* pi */
double PIO2   =  1.57079632679489661923;       /* pi/2 */
double PIO4   =  7.85398163397448309616E-1;    /* pi/4 */
double SQRT2  =  1.41421356237309504880;       /* sqrt(2) */
double SQRTH  =  7.07106781186547524401E-1;    /* sqrt(2)/2 */
double LOG2E  =  1.4426950408889634073599;     /* 1/log(2) */
double SQ2OPI =  7.9788456080286535587989E-1;  /* sqrt( 2/pi ) */
double LOGE2  =  6.93147180559945309417E-1;    /* log(2) */
double LOGSQ2 =  3.46573590279972654709E-1;    /* log(2)/2 */
double THPIO4 =  2.35619449019234492885;       /* 3*pi/4 */
double TWOOPI =  6.36619772367581343075535E-1; /* 2/pi */
#ifdef MINUSZERO
double NEGZERO = -0.0;
#else
double NEGZERO = 0.0;
#endif

#define MAXGAM 171.624376956302725
#define MAXLGM 2.556348e305
static double LS2PI  =  0.91893853320467274178;

static inline double ndtri (double y0);
double igami(double a, double y0);
double lgam(double x);
double igaml(double a, double x);
static double stirf(double);
double polevl(double x, void *p, int n);
double p1evl(double x, void *p, int n);
static double stirf (double);
//extern int isnan(double x);
//extern int isfinite(double x);

double imgamil(double a, double y0)
{
double x0, x1, x, yl, yh, y, d, lgm, dithresh;
int i, dir;

/* bound the solution */
x0 = MAXNUM;
yl = 0.0;
x1 = 0.0;
yh = 1.0;
dithresh = 4.0 * MACHEP;

/* approximation to inverse function */
d = 1.0/(9.0*a);
y = (1.0 - d - ndtri(y0) * sqrt(d));
x = a * y * y * y;

lgm = lgam(a);

for(i=0; i<10; i++)
    {
    if(x > x0 || x < x1)
        goto ihalve;
    y = igamc(a,x);
    if(y < yl || y > yh)
        goto ihalve;
    if(y < y0)
        {
        x0 = x;
        yl = y;
        }
    else
        {
        x1 = x;
        yh = y;
        }
/* compute the derivative of the function at this point */
    d = (a - 1.0) * log(x0) - x0 - lgm;
    if(d < -MAXLOG)
        goto ihalve;
    d = -exp(d);
/* compute the step to the next approximation of x */
    d = (y - y0)/d;
    x = x - d;
    if(i < 3)
        continue;
    if(fabs(d/x) < dithresh)
        goto done;
    }

/* Resort to interval halving if Newton iteration did not converge. */
ihalve:

d = 0.0625;
if(x0 == MAXNUM)
    {
    if(x <= 0.0)
        x = 1.0;
    while(x0 == MAXNUM)
        {
        x = (1.0 + d) * x;
        y = igamc(a, x);
        if(y < y0)
            {
            x0 = x;
            yl = y;
            break;
            }
        d = d + d;
        }
    }
d = 0.5;
dir = 0;

for(i=0; i<400; i++)
    {
    x = x1  +  d * (x0 - x1);
    y = igamc(a, x);
    lgm = (x0 - x1)/(x1 + x0);
    if(fabs(lgm) < dithresh)
        break;
    lgm = (y - y0)/y0;
    if(fabs(lgm) < dithresh)
        break;
    if(x <= 0.0)
        break;
    if(y > y0)
        {
        x1 = x;
        yh = y;
        if(dir < 0)
            {
            dir = 0;
            d = 0.5;
            }
        else if(dir > 1)
            d = 0.5 * d + 0.5;
        else
            d = (y0 - yl)/(yh - yl);
        dir += 1;
        }
    else
        {
        x0 = x;
        yl = y;
        if(dir > 0)
            {
            dir = 0;
            d = 0.5;
            }
        else if(dir < -1)
            d = 0.5 * d;
        else
            d = (y0 - yl)/(yh - yl);
        dir -= 1;
        }
    }
if(x == 0.0)
    mtherr("igami", UNDERFLOW);

done:
return(x);
}

/*                            igaml.c
 *
 *    Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igaml();
 *
 * y = igaml(a, x);
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         4000       4.4e-15     6.3e-16
 *    IEEE      0,30        10000       3.6e-14     5.1e-15
 *
 */
/*                            igamc()
 *
 *    Complemented incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igamc();
 *
 * y = igamc(a, x);
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0,30         2000       2.7e-15     4.0e-16
 *    IEEE      0,30        60000       1.4e-12     6.3e-15
 *
 */

/*
Cephes Math Library Release 2.3:  March, 1995
Copyright 1985, 1995 by Stephen L. Moshier
*/

#define BIG 9.223372036854775808e18

double igamc(double a, double x)
{
double ans, c, yc, ax, y, z, r, t;
double pk, pkm1, pkm2, qk, qkm1, qkm2;

if((x <= 0.0) || (a <= 0.0))
    return(1.0);

if((x < 1.0) || (x < a))
    return(1.0 - igaml(a,x));

ax = a * log(x) - x - lgam(a);
if(ax < MINLOG)
    {
#if !NDEBUG
    mtherr("igamc", UNDERFLOW);
#endif
    //return(0.0);
    return(MINLOG);
    }
ax = exp(ax);

/* continued fraction */
y = 1.0 - a;
z = x + y + 1.0;
c = 0.0;
pkm2 = 1.0;
qkm2 = x;
pkm1 = x + 1.0;
qkm1 = z * x;
ans = pkm1/qkm1;

do
    {
    c += 1.0;
    y += 1.0;
    z += 2.0;
    yc = y * c;
    pk = pkm1 * z  -  pkm2 * yc;
    qk = qkm1 * z  -  qkm2 * yc;
    if(qk != 0.0)
        {
        r = pk/qk;
        t = fabs((ans - r)/r);
        ans = r;
        }
    else
        t = 1.0;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
    if(fabs(pk) > BIG)
        {
        pkm2 /= BIG;
        pkm1 /= BIG;
        qkm2 /= BIG;
        qkm1 /= BIG;
        }
    }
while(t > MACHEP);
return(ans * ax);
}



/* left tail of incomplete gamma function:
 *
 *          inf.      k
 *   a  -x   -       x
 *  x  e     >   ----------
 *           -     -
 *          k=0   | (a+k+1)
 *
 */

double igaml(double a, double x)
{
double ans, ax, c, r;

if((x <= 0.0) || (a <= 0.0))
    return(0.0);

if((x > 1.0) && (x > a))
    return(1.0 - igamc(a,x));

ax = a * log(x) - x - lgam(a);
if(ax < MINLOG)
    {
    mtherr("igaml", UNDERFLOW);
    //return(0.0);
    return(MINLOG);
    }
ax = exp(ax);

/* power series */
r = a;
c = 1.0;
ans = 1.0;

do
    {
    r += 1.0;
    c *= x/r;
    ans += c;
    }
while(c/ans > MACHEP);

return(ans * ax/a);
}

/*                            ndtri.c
 *
 *    Inverse of Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, ndtri();
 *
 * x = ndtri(y);
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the argument, x, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to x) is equal to y.
 *
 *
 * For small arguments 0 < y < exp(-2), the program computes
 * z = sqrt(-2 log(y));  then the approximation is
 * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z) .
 * For larger arguments,  x/sqrt(2 pi) = w + w^3 R(w^2)/S(w^2)) ,
 * where w = y - 0.5 .
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain        # trials      peak         rms
 *  Arguments uniformly distributed:
 *    IEEE       0, 1           5000       7.8e-19     9.9e-20
 *  Arguments exponentially distributed:
 *    IEEE     exp(-11355),-1  30000       1.7e-19     4.3e-20
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition    value returned
 * ndtri domain      x <= 0        -MAXNUM
 * ndtri domain      x >= 1         MAXNUM
 *
 */


/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

/*
 * See the header file igamc_cephes.h for the implementation of ndtri.
 */

/*                            gamma.c
 *
 *    Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, gamma();
 * extern int sgngam;
 *
 * y = gamma(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngam.
 * This variable is also filled in by the logarithmic gamma
 * function lgam().
 *
 * Arguments |x| <= 13 are reduced by recurrence and the function
 * approximated by a rational function of degree 7/8 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -40,+40      10000       3.6e-19     7.9e-20
 *    IEEE    -1755,+1755   10000       4.8e-18     6.5e-19
 *
 * Accuracy for large arguments is dominated by error in powl().
 *
 */
/*                            lgam()
 *
 *    Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 * extern int sgngam;
 *
 * y = lgam(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngam.
 *
 * For arguments greater than 33, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGML (10^4928) return MAXNUM.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    IEEE         -40, 40        100000     2.2e-19     4.6e-20
 *    IEEE    10^-2000,10^+2000    20000     1.6e-19     3.3e-20
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 */

/*                            gamma.c    */
/*    gamma function    */

/*
Copyright 1994 by Stephen L. Moshier
*/


/*
gamma(x+2)  = gamma(x+2) P(x)/Q(x)
0 <= x <= 1
Relative error
n=7, d=8
Peak error =  1.83e-20
Relative error spread =  8.4e-23
*/
/*static double LOGPI = 1.14472988584940017414L;*/

/* Stirling's formula for the gamma function
gamma(x) = sqrt(2 pi) x^(x-.5) exp(-x) (1 + 1/x P(1/x))
z(x) = x
13 <= x <= 1024
Relative error
n=8, d=0
Peak error =  9.44e-21
Relative error spread =  8.8e-4
*/

int sgngam = 0;
extern int sgngam;

#ifdef __cplusplus
double gamma(double x) throw()
#else
double gamma(double x)
#endif
{
double p, q, z;
int i;

sgngam = 1;
#ifdef NANS
if( isnan(x) )
    return(x);
#endif
#ifdef INFINITIES
#ifdef NANS
if( x == INFINITY )
    return(x);
if( x == -INFINITY )
    return(NAN);
#else
if( !isfinite(x) )
    return(x);
#endif
#endif
q = fabs(x);

if( q > 33.0 )
    {
    if( x < 0.0 )
        {
        p = floor(q);
        if( p == q )
            {
#ifdef NANS
gamnan:
            mtherr( "gamma", DOMAIN );
            return (NAN);
#else
            goto goverf;
#endif
            }
        i = p;
        if( (i & 1) == 0 )
            sgngam = -1;
        z = q - p;
        if( z > 0.5 )
            {
            p += 1.0;
            z = q - p;
            }
        z = q * sin( PI * z );
        if( z == 0.0 )
            {
#ifdef INFINITIES
            return( sgngam * INFINITY);
#else
goverf:
            mtherr( "gamma", OVERFLOW );
            return( sgngam * MAXNUM);
#endif
            }
        z = fabs(z);
        z = PI/(z * stirf(q) );
        }
    else
        {
        z = stirf(x);
        }
    return( sgngam * z );
    }

z = 1.0;
while( x >= 3.0 )
    {
    x -= 1.0;
    z *= x;
    }

while( x < 0.0 )
    {
    if( x > -1.E-9 )
        goto small;
    z /= x;
    x += 1.0;
    }

while( x < 2.0 )
    {
    if( x < 1.e-9 )
        goto small;
    z /= x;
    x += 1.0;
    }

if( x == 2.0 )
    return(z);

x -= 2.0;
p = polevl( x, P, 6 );
q = polevl( x, Q, 7 );
return( z * p / q );

small:
if( x == 0.0 )
    {
#ifdef INFINITIES
#ifdef NANS
      goto gamnan;
#else
      return( INFINITY );
#endif
#else
    mtherr( "gamma", SING );
    return( MAXNUM );
#endif
    }
else
    return( z/((1.0 + 0.5772156649015329 * x) * x) );
}



/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */

static double A[] = {
 8.11614167470508450300E-4,
-5.95061904284301438324E-4,
 7.93650340457716943945E-4,
-2.77777777730099687205E-3,
 8.33333333333331927722E-2
};
static double B[] = {
-1.37825152569120859100E3,
-3.88016315134637840924E4,
-3.31612992738871184744E5,
-1.16237097492762307383E6,
-1.72173700820839662146E6,
-8.53555664245765465627E5
};
static double C[] = {
/* 1.00000000000000000000E0, */
-3.51815701436523470549E2,
-1.70642106651881159223E4,
-2.20528590553854454839E5,
-1.13933444367982507207E6,
-2.53252307177582951285E6,
-2.01889141433532773231E6
};


/* log gamma(x) = (x - 0.5) * log(x) - x + LS2PI + 1/x A(1/x^2)
 * x >= 8
 * Peak relative error 1.51e-21
 * Relative spread of error peaks 5.67e-21
 */

/* log(sqrt(2*pi)) */


/* Logarithm of gamma function */


double lgam(double x)
{
    double p, q, u, w, z;
    int i;

    sgngam = 1;
    #ifdef NANS
    if( isnan(x) )
        return(x);
    #endif

    #ifdef INFINITIES
    if( !isfinite(x) )
        return(INFINITY);
    #endif

    if( x < -34.0 )
        {
        q = -x;
        w = lgam(q); /* note this modifies sgngam! */
        p = floor(q);
        if( p == q )
            {
    lgsing:
    #ifdef INFINITIES
            mtherr( "lgam", SING );
            return (INFINITY);
    #else
            goto loverf;
    #endif
            }
        i = p;
        if( (i & 1) == 0 )
            sgngam = -1;
        else
            sgngam = 1;
        z = q - p;
        if( z > 0.5 )
            {
            p += 1.0;
            z = p - q;
            }
        z = q * sin( PI * z );
        if( z == 0.0 )
            goto lgsing;
    /*    z = log(PI) - log( z ) - w;*/
        z = LOGPI - log( z ) - w;
        return( z );
        }

    if( x < 13.0 )
        {
        z = 1.0;
        p = 0.0;
        u = x;
        while( u >= 3.0 )
            {
                p -= 1.0;
                u = x + p;
                z *= u;
            }
        while( u < 2.0 )
            {
                if( u == 0.0 )
                    goto lgsing;
                z /= u;
                p += 1.0;
                u = x + p;
            }
        if( z < 0.0 )
            {
                sgngam = -1;
                z = -z;
            }
        else
            sgngam = 1;
        if( u == 2.0 )
            return( log(z) );
        p -= 2.0;
        x = x + p;
        p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
        return( log(z) + p );
        }

    if( x > MAXLGM )
        {
    #ifdef INFINITIES
        return( sgngam * INFINITY );
    #else
    loverf:
        mtherr( "lgam", OVERFLOW );
        return( sgngam * MAXNUM );
    #endif
        }

    q = ( x - 0.5 ) * log(x) - x + LS2PI;
    if( x > 1.0e8 )
        return( q );

    p = 1.0/(x*x);
    if( x >= 1000.0 )
        q += ((   7.9365079365079365079365e-4 * p
            - 2.7777777777777777777778e-3) *p
            + 0.0833333333333333333333) / x;
    else
        q += polevl( p, A, 4 ) / x;
    return( q );
}
