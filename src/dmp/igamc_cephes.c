/*							igami()
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

static double P[] = {
  1.60119522476751861407E-4,
  1.19135147006586384913E-3,
  1.04213797561761569935E-2,
  4.76367800457137231464E-2,
  2.07448227648435975150E-1,
  4.94214826801497100753E-1,
  9.99999999999999996796E-1
};
static double Q[] = {
-2.31581873324120129819E-5,
 5.39605580493303397842E-4,
-4.45641913851797240494E-3,
 1.18139785222060435552E-2,
 3.58236398605498653373E-2,
-2.34591795718243348568E-1,
 7.14304917030273074085E-2,
 1.00000000000000000320E0
};
#define MAXGAM 171.624376956302725
#define MAXLGM 2.556348e305
static double LS2PI  =  0.91893853320467274178;

double igamc(double a, double x);
double ndtri (double y0);
double igami(double a, double y0);
double lgam(double x);
double igaml(double a, double x);
static double stirf(double);
double polevl(double x, void *p, int n);
double p1evl(double x, void *p, int n);
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

/*							igaml.c
 *
 *	Incomplete gamma integral
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

double igamc(a, x)
double a, x;
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

/*							ndtri.c
 *
 *	Inverse of Normal distribution function
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

/*							gamma.c
 *
 *	Gamma function
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
/*							lgam()
 *
 *	Natural logarithm of gamma function
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

/*							gamma.c	*/
/*	gamma function	*/

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
#if UNK
static double STIR[9] = {
 7.147391378143610789273E-4L,
-2.363848809501759061727E-5L,
-5.950237554056330156018E-4L,
 6.989332260623193171870E-5L,
 7.840334842744753003862E-4L,
-2.294719747873185405699E-4L,
-2.681327161876304418288E-3L,
 3.472222222230075327854E-3L,
 8.333333333333331800504E-2L,
};
#endif
#if IBMPC
static short STIR[] = {
0x6ede,0x69f7,0x54e3,0xbb5d,0x3ff4, XPD
0xc395,0x0295,0x4443,0xc64b,0xbfef, XPD
0xba6f,0x7c59,0x5e47,0x9bfb,0xbff4, XPD
0x5704,0x1a39,0xb11d,0x9293,0x3ff1, XPD
0x30b7,0x1a21,0x98b2,0xcd87,0x3ff4, XPD
0xbef3,0x7023,0x6a08,0xf09e,0xbff2, XPD
0x3a1c,0x5ac8,0x3478,0xafb9,0xbff6, XPD
0xc3c9,0x906e,0x38e3,0xe38e,0x3ff6, XPD
0xa1d5,0xaaaa,0xaaaa,0xaaaa,0x3ffb, XPD
};
#endif
#if MIEEE
static long STIR[27] = {
0x3ff40000,0xbb5d54e3,0x69f76ede,
0xbfef0000,0xc64b4443,0x0295c395,
0xbff40000,0x9bfb5e47,0x7c59ba6f,
0x3ff10000,0x9293b11d,0x1a395704,
0x3ff40000,0xcd8798b2,0x1a2130b7,
0xbff20000,0xf09e6a08,0x7023bef3,
0xbff60000,0xafb93478,0x5ac83a1c,
0x3ff60000,0xe38e38e3,0x906ec3c9,
0x3ffb0000,0xaaaaaaaa,0xaaaaa1d5,
};
#endif
#define MAXSTIR 1024.0
static double SQTPI = 2.50662827463100050242E0L;

/* 1/gamma(x) = z P(z)
 * z(x) = 1/x
 * 0 < x < 0.03125
 * Peak relative error 4.2e-23
 */
#if UNK
static double S[9] = {
-1.193945051381510095614E-3L,
 7.220599478036909672331E-3L,
-9.622023360406271645744E-3L,
-4.219773360705915470089E-2L,
 1.665386113720805206758E-1L,
-4.200263503403344054473E-2L,
-6.558780715202540684668E-1L,
 5.772156649015328608253E-1L,
 1.000000000000000000000E0L,
};
#endif
#if IBMPC
static short S[] = {
0xbaeb,0xd6d3,0x25e5,0x9c7e,0xbff5, XPD
0xfe9a,0xceb4,0xc74e,0xec9a,0x3ff7, XPD
0x9225,0xdfef,0xb0e9,0x9da5,0xbff8, XPD
0x10b0,0xec17,0x87dc,0xacd7,0xbffa, XPD
0x6b8d,0x7515,0x1905,0xaa89,0x3ffc, XPD
0xf183,0x126b,0xf47d,0xac0a,0xbffa, XPD
0x7bf6,0x57d1,0xa013,0xa7e7,0xbffe, XPD
0xc7a9,0x7db0,0x67e3,0x93c4,0x3ffe, XPD
0x0000,0x0000,0x0000,0x8000,0x3fff, XPD
};
#endif
#if MIEEE
static long S[27] = {
0xbff50000,0x9c7e25e5,0xd6d3baeb,
0x3ff70000,0xec9ac74e,0xceb4fe9a,
0xbff80000,0x9da5b0e9,0xdfef9225,
0xbffa0000,0xacd787dc,0xec1710b0,
0x3ffc0000,0xaa891905,0x75156b8d,
0xbffa0000,0xac0af47d,0x126bf183,
0xbffe0000,0xa7e7a013,0x57d17bf6,
0x3ffe0000,0x93c467e3,0x7db0c7a9,
0x3fff0000,0x80000000,0x00000000,
};
#endif
/* 1/gamma(-x) = z P(z)
 * z(x) = 1/x
 * 0 < x < 0.03125
 * Peak relative error 5.16e-23
 * Relative error spread =  2.5e-24
 */
#if UNK
static double SN[9] = {
 1.133374167243894382010E-3L,
 7.220837261893170325704E-3L,
 9.621911155035976733706E-3L,
-4.219773343731191721664E-2L,
-1.665386113944413519335E-1L,
-4.200263503402112910504E-2L,
 6.558780715202536547116E-1L,
 5.772156649015328608727E-1L,
-1.000000000000000000000E0L,
};
#endif
#if IBMPC
static short SN[] = {
0x5dd1,0x02de,0xb9f7,0x948d,0x3ff5, XPD
0x989b,0xdd68,0xc5f1,0xec9c,0x3ff7, XPD
0x2ca1,0x18f0,0x386f,0x9da5,0x3ff8, XPD
0x783f,0x41dd,0x87d1,0xacd7,0xbffa, XPD
0x7a5b,0xd76d,0x1905,0xaa89,0xbffc, XPD
0x7f64,0x1234,0xf47d,0xac0a,0xbffa, XPD
0x5e26,0x57d1,0xa013,0xa7e7,0x3ffe, XPD
0xc7aa,0x7db0,0x67e3,0x93c4,0x3ffe, XPD
0x0000,0x0000,0x0000,0x8000,0xbfff, XPD
};
#endif
#if MIEEE
static long SN[27] = {
0x3ff50000,0x948db9f7,0x02de5dd1,
0x3ff70000,0xec9cc5f1,0xdd68989b,
0x3ff80000,0x9da5386f,0x18f02ca1,
0xbffa0000,0xacd787d1,0x41dd783f,
0xbffc0000,0xaa891905,0xd76d7a5b,
0xbffa0000,0xac0af47d,0x12347f64,
0x3ffe0000,0xa7e7a013,0x57d15e26,
0x3ffe0000,0x93c467e3,0x7db0c7aa,
0xbfff0000,0x80000000,0x00000000,
};
#endif

int sgngaml = 0;
extern int sgngaml;


static double stirf (double);

/* Gamma function computed by Stirling's formula.
 */
static double stirf(x)
double x;
{
double y, w, v;

w = 1.0/x;
/* For large x, use rational coefficients from the analytical expansion.  */
if(x > 1024.0)
	w = (((((6.97281375836585777429E-5 * w
		+ 7.84039221720066627474E-4) * w
		- 2.29472093621399176955E-4) * w
		- 2.68132716049382716049E-3) * w
		+ 3.47222222222222222222E-3) * w
		+ 8.33333333333333333333E-2) * w
		+ 1.0;
else
	w = 1.0 + w * polevl(w, STIR, 8);
y = expl(x);
if(x > MAXSTIR)
	{ /* Avoid overflow in pow() */
	v = pow(x, 0.5 * x - 0.25);
	y = v * (v / y);
	}
else
	{
	y = pow(x, x - 0.5) / y;
	}
y = SQTPI * y * w;
return(y);
}



double gamma(double x)
{
double p, q, z;
int i;

sgngaml = 1;
#ifdef NANS
if(isnan(x))
	return(NAN);
#endif
#ifdef INFINITIES
if(x == INFINITY)
	return(INFINITY);
#ifdef NANS
if(x == -INFINITY)
	goto gamnan;
#endif
#endif
q = fabsl(x);

if(q > 13.0L)
	{
	if(q > MAXGAM)
		goto goverf;
	if(x < 0.0L)
		{
		p = floor(q);
		if(p == q)
			{
gamnan:
#ifdef NANS
			mtherr("gamma", DOMAIN);
			return (NAN);
#else
			goto goverf;
#endif
			}
		i = p;
		if((i & 1) == 0)
			sgngaml = -1;
		z = q - p;
		if(z > 0.5L)
			{
			p += 1.0L;
			z = q - p;
			}
		z = q * sin(PI * z);
		z = fabsl(z) * stirf(q);
		if(z <= PI/MAXNUM)
			{
goverf:
#ifdef INFINITIES
			return(sgngaml * INFINITY);
#else
			mtherr("gamma", OVERFLOW);
			return(sgngaml * MAXNUM);
#endif
			}
		z = PI/z;
		}
	else
		{
		z = stirf(x);
		}
	return(sgngaml * z);
	}

z = 1.0L;
while(x >= 3.0L)
	{
	x -= 1.0L;
	z *= x;
	}

while(x < -0.03125L)
	{
	z /= x;
	x += 1.0L;
	}

if(x <= 0.03125L)
	goto small;

while(x < 2.0L)
	{
	z /= x;
	x += 1.0L;
	}

if(x == 2.0L)
	return(z);

x -= 2.0L;
p = polevl(x, P, 7);
q = polevl(x, Q, 8);
return(z * p / q);

small:
if(x == 0.0L)
	{
	  goto gamnan;
	}
else
	{
	if(x < 0.0L)
		{
		x = -x;
		q = z / (x * polevl(x, SN, 8));
		}
	else
		q = z / (x * polevl(x, S, 8));
	}
return q;
}



/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */


/* log gamma(x) = (x - 0.5) * log(x) - x + LS2PI + 1/x A(1/x^2)
 * x >= 8
 * Peak relative error 1.51e-21
 * Relative spread of error peaks 5.67e-21
 */
#if UNK
static double A[7] = {
 4.885026142432270781165E-3L,
-1.880801938119376907179E-3L,
 8.412723297322498080632E-4L,
-5.952345851765688514613E-4L,
 7.936507795855070755671E-4L,
-2.777777777750349603440E-3L,
 8.333333333333331447505E-2L,
};
#endif
#if IBMPC
static short A[] = {
0xd984,0xcc08,0x91c2,0xa012,0x3ff7, XPD
0x3d91,0x0304,0x3da1,0xf685,0xbff5, XPD
0x3bdc,0xaad1,0xd492,0xdc88,0x3ff4, XPD
0x8b20,0x9fce,0x844e,0x9c09,0xbff4, XPD
0xf8f2,0x30e5,0x0092,0xd00d,0x3ff4, XPD
0x4d88,0x03a8,0x60b6,0xb60b,0xbff6, XPD
0x9fcc,0xaaaa,0xaaaa,0xaaaa,0x3ffb, XPD
};
#endif
#if MIEEE
static long A[21] = {
0x3ff70000,0xa01291c2,0xcc08d984,
0xbff50000,0xf6853da1,0x03043d91,
0x3ff40000,0xdc88d492,0xaad13bdc,
0xbff40000,0x9c09844e,0x9fce8b20,
0x3ff40000,0xd00d0092,0x30e5f8f2,
0xbff60000,0xb60b60b6,0x03a84d88,
0x3ffb0000,0xaaaaaaaa,0xaaaa9fcc,
};
#endif

/* log gamma(x+2) = x B(x)/C(x)
 * 0 <= x <= 1
 * Peak relative error 7.16e-22
 * Relative spread of error peaks 4.78e-20
 */
#if UNK
static double B[7] = {
-2.163690827643812857640E3L,
-8.723871522843511459790E4L,
-1.104326814691464261197E6L,
-6.111225012005214299996E6L,
-1.625568062543700591014E7L,
-2.003937418103815175475E7L,
-8.875666783650703802159E6L,
};
static double C[7] = {
/* 1.000000000000000000000E0L,*/
-5.139481484435370143617E2L,
-3.403570840534304670537E4L,
-6.227441164066219501697E5L,
-4.814940379411882186630E6L,
-1.785433287045078156959E7L,
-3.138646407656182662088E7L,
-2.099336717757895876142E7L,
};
#endif
#if IBMPC
static short B[] = {
0x9557,0x4995,0x0da1,0x873b,0xc00a, XPD
0xfe44,0x9af8,0x5b8c,0xaa63,0xc00f, XPD
0x5aa8,0x7cf5,0x3684,0x86ce,0xc013, XPD
0x259a,0x258c,0xf206,0xba7f,0xc015, XPD
0xbe18,0x1ca3,0xc0a0,0xf80a,0xc016, XPD
0x168f,0x2c42,0x6717,0x98e3,0xc017, XPD
0x2051,0x9d55,0x92c8,0x876e,0xc016, XPD
};
static short C[] = {
/*0x0000,0x0000,0x0000,0x8000,0x3fff, XPD*/
0xaa77,0xcf2f,0xae76,0x807c,0xc008, XPD
0xb280,0x0d74,0xb55a,0x84f3,0xc00e, XPD
0xa505,0xcd30,0x81dc,0x9809,0xc012, XPD
0x3369,0x4246,0xb8c2,0x92f0,0xc015, XPD
0x63cf,0x6aee,0xbe6f,0x8837,0xc017, XPD
0x26bb,0xccc7,0xb009,0xef75,0xc017, XPD
0x462b,0xbae8,0xab96,0xa02a,0xc017, XPD
};
#endif
#if MIEEE
static long B[21] = {
0xc00a0000,0x873b0da1,0x49959557,
0xc00f0000,0xaa635b8c,0x9af8fe44,
0xc0130000,0x86ce3684,0x7cf55aa8,
0xc0150000,0xba7ff206,0x258c259a,
0xc0160000,0xf80ac0a0,0x1ca3be18,
0xc0170000,0x98e36717,0x2c42168f,
0xc0160000,0x876e92c8,0x9d552051,
};
static long C[21] = {
/*0x3fff0000,0x80000000,0x00000000,*/
0xc0080000,0x807cae76,0xcf2faa77,
0xc00e0000,0x84f3b55a,0x0d74b280,
0xc0120000,0x980981dc,0xcd30a505,
0xc0150000,0x92f0b8c2,0x42463369,
0xc0170000,0x8837be6f,0x6aee63cf,
0xc0170000,0xef75b009,0xccc726bb,
0xc0170000,0xa02aab96,0xbae8462b,
};
#endif

/* log(sqrt(2*pi)) */



/* Logarithm of gamma function */


double lgam(double x)
{
double p, q, w, z, f, nx;
int i;

sgngaml = 1;
#ifdef NANS
if(isnan(x))
	return(NAN);
#endif
#ifdef INFINITIES
if(!isfinite(x))
	return(INFINITY);
#endif
if(x < -34.0L)
	{
	q = -x;
	w = lgam(q); /* note this modifies sgngam! */
	p = floor(q);
	if(p == q)
		{
#ifdef INFINITIES
		mtherr("lgam", SING);
		return (INFINITY);
#else
		goto loverf;
#endif
		}
	i = p;
	if((i & 1) == 0)
		sgngaml = -1;
	else
		sgngaml = 1;
	z = q - p;
	if(z > 0.5L)
		{
		p += 1.0L;
		z = p - q;
		}
	z = q * sin(PI * z);
	if(z == 0.0L)
		goto loverf;
/*	z = LOGPI - log(z) - w; */
	z = log(PI/z) - w;
	return(z);
	}

if(x < 13.0L)
	{
	z = 1.0L;
	nx = floor(x +  0.5L);
	f = x - nx;
	while(x >= 3.0L)
		{
		nx -= 1.0L;
		x = nx + f;
		z *= x;
		}
	while(x < 2.0L)
		{
		if(fabsl(x) <= 0.03125)
			goto lsmall;
		z /= nx +  f;
		nx += 1.0L;
		x = nx + f;
		}
	if(z < 0.0L)
		{
		sgngaml = -1;
		z = -z;
		}
	else
		sgngaml = 1;
	if(x == 2.0L)
		return(log(z));
	x = (nx - 2.0L) + f;
	p = x * polevl(x, B, 6) / p1evl(x, C, 7);
	return(log(z) + p);
	}

if(x > MAXLGM)
	{
loverf:
#ifdef INFINITIES
	return(sgngaml * INFINITY);
#else
	mtherr("lgam", OVERFLOW);
	return(sgngaml * MAXNUM);
#endif
	}

q = (x - 0.5L) * log(x) - x + LS2PI;
if(x > 1.0e10L)
	return(q);
p = 1.0L/(x*x);
q += polevl(p, A, 6) / x;
return(q);


lsmall:
if(x == 0.0L)
	goto loverf;
if(x < 0.0L)
	{
	x = -x;
	q = z / (x * polevl(x, SN, 8));
	}
else
	q = z / (x * polevl(x, S, 8));
if(q < 0.0L)
	{
	sgngaml = -1;
	q = -q;
	}
else
	sgngaml = 1;
q = log(q);
return(q);
}

