/*							igam.c
 *
 *	Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igam();
 *
 * y = igam( a, x );
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
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30       200000       3.6e-14     2.9e-15
 *    IEEE      0,100      300000       9.9e-14     1.5e-14
 */
/*							igamc()
 *
 *	Complemented incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * double a, x, y, igamc();
 *
 * y = igamc( a, x );
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
 * ACCURACY:
 *
 * Tested at random a, x.
 *                a         x                      Relative error:
 * arithmetic   domain   domain     # trials      peak         rms
 *    IEEE     0.5,100   0,100      200000       1.9e-14     1.7e-15
 *    IEEE     0.01,0.5  0,100      200000       1.4e-13     1.6e-15
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include "mconf.h"
#ifndef ANSIPROT
double lgam(), exp(), log(), fabs(), igam(), igamc();
#endif
#ifdef UNK
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
#ifdef INFINITIES
double INFINITY = 1.0/0.0;  /* 99e999; */
#else
double INFINITY =  1.79769313486231570815E308;    /* 2**1024*(1-MACHEP) */
#endif
#ifdef NANS
double NAN = 1.0/0.0 - 1.0/0.0;
#else
double NAN = 0.0;
#endif
#ifdef MINUSZERO
double NEGZERO = -0.0;
#else
double NEGZERO = 0.0;
#endif
#endif
extern double MACHEP;
static double MAXLOG =  7.09782712893383996732E2;
static double big = 4.503599627370496e15;
static double biginv =  2.22044604925031308085e-16;

double igamc( a, x )
double a, x;
{
double ans, ax, c, yc, r, t, y, z;
double pk, pkm1, pkm2, qk, qkm1, qkm2;

if( (x <= 0) || ( a <= 0) )
	return( 1.0 );

if( (x < 1.0) || (x < a) )
	return( 1.0 - igam(a,x) );

ax = a * log(x) - x - lgam(a);
if( ax < -MAXLOG )
	{
	mtherr( "igamc", UNDERFLOW );
	return( 0.0 );
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
	if( qk != 0 )
		{
		r = pk/qk;
		t = fabs( (ans - r)/r );
		ans = r;
		}
	else
		t = 1.0;
	pkm2 = pkm1;
	pkm1 = pk;
	qkm2 = qkm1;
	qkm1 = qk;
	if( fabs(pk) > big )
		{
		pkm2 *= biginv;
		pkm1 *= biginv;
		qkm2 *= biginv;
		qkm1 *= biginv;
		}
	}
while( t > MACHEP );

return( ans * ax );
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

double igam( a, x )
double a, x;
{
double ans, ax, c, r;

if( (x <= 0) || ( a <= 0) )
	return( 0.0 );

if( (x > 1.0) && (x > a ) )
	return( 1.0 - igamc(a,x) );

/* Compute  x**a * exp(-x) / gamma(a)  */
ax = a * log(x) - x - lgam(a);
if( ax < -MAXLOG )
	{
	mtherr( "igam", UNDERFLOW );
	return( 0.0 );
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
while( c/ans > MACHEP );

return( ans * ax/a );
}
