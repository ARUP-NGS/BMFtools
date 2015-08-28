/*							igamil()
 *
 *      Inverse of complemented imcomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float128_t a, x, y, igamil();
 *
 * x = igamil(a, y);
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
#include <quadmath.h>
#include "igamc_cephes.h"

static float128_t LS2PI  =  0.91893853320467274178L;

//Always set it this way.
//#if UNK
/* almost 2^16384 */
float128_t MAXNUML = 1.189731495357231765021263853E4932L;
/* 2^-64 */
float128_t MACHEPL = 5.42101086242752217003726400434970855712890625E-20L;
/* log(MAXNUML) */
float128_t MAXLOGL =  1.1356523406294143949492E4L;
#ifdef DENORMAL
/* log(smallest denormal number = 2^-16446) */
float128_t MINLOGL = -1.13994985314888605586758E4L;
#else
/* log(underflow threshold = 2^(-16382)) */
float128_t MINLOGL = -1.1355137111933024058873E4L;
#endif
float128_t LOGE2L  = 6.9314718055994530941723E-1L;
float128_t LOG2EL  = 1.4426950408889634073599E0L;
float128_t PIL     = 3.1415926535897932384626L;
float128_t PIO2L   = 1.5707963267948966192313L;
float128_t PIO4L   = 7.8539816339744830961566E-1L;
#ifdef INFINITIES
float128_t NANL = 0.0L / 0.0L;
float128_t INFINITYL = 1.0L / 0.0L;
#else
float128_t INFINITYL = 1.189731495357231765021263853E4932L;
float128_t NANL = 0.0L;
#endif
//#endif
#if IBMPC
short MAXNUML[] = {0xffff,0xffff,0xffff,0xffff,0x7ffe, XPD};
short MAXLOGL[] = {0x79ab,0xd1cf,0x17f7,0xb172,0x400c, XPD};
#ifdef INFINITIES
short INFINITYL[] = {0,0,0,0x8000,0x7fff, XPD};
short NANL[] = {0,0,0,0xc000,0x7fff, XPD};
#else
short INFINITYL[] = {0xffff,0xffff,0xffff,0xffff,0x7ffe, XPD};
float128_t NANL = 0.0L;
#endif
#ifdef DENORMAL
short MINLOGL[] = {0xbaaa,0x09e2,0xfe7f,0xb21d,0xc00c, XPD};
#else
short MINLOGL[] = {0xeb2f,0x1210,0x8c67,0xb16c,0xc00c, XPD};
#endif
short MACHEPL[] = {0x0000,0x0000,0x0000,0x8000,0x3fbf, XPD};
short LOGE2L[]  = {0x79ac,0xd1cf,0x17f7,0xb172,0x3ffe, XPD};
short LOG2EL[]  = {0xf0bc,0x5c17,0x3b29,0xb8aa,0x3fff, XPD};
short PIL[]     = {0xc235,0x2168,0xdaa2,0xc90f,0x4000, XPD};
short PIO2L[]   = {0xc235,0x2168,0xdaa2,0xc90f,0x3fff, XPD};
short PIO4L[]   = {0xc235,0x2168,0xdaa2,0xc90f,0x3ffe, XPD};
#endif
#if MIEEE
long MAXNUML[] = {0x7ffe0000,0xffffffff,0xffffffff};
long MAXLOGL[] = {0x400c0000,0xb17217f7,0xd1cf79ab};
#ifdef INFINITIES
long INFINITY[] = {0x7fff0000,0x80000000,0x00000000};
long NANL[] = {0x7fff0000,0xffffffff,0xffffffff};
#else
long INFINITYL[] = {0x7ffe0000,0xffffffff,0xffffffff};
float128_t NANL = 0.0L;
#endif
#ifdef DENORMAL
long MINLOGL[] = {0xc00c0000,0xb21dfe7f,0x09e2baaa};
#else
long MINLOGL[] = {0xc00c0000,0xb16c8c67,0x1210eb2f};
#endif
long MACHEPL[] = {0x3fbf0000,0x80000000,0x00000000};
long LOGE2L[]  = {0x3ffe0000,0xb17217f7,0xd1cf79ac};
long LOG2EL[]  = {0x3fff0000,0xb8aa3b29,0x5c17f0bc};
long PIL[]     = {0x40000000,0xc90fdaa2,0x2168c235};
long PIO2L[]   = {0x3fff0000,0xc90fdaa2,0x2168c235};
long PIO4L[]   = {0x3ffe0000,0xc90fdaa2,0x2168c235};
#endif

#ifdef MINUSZERO
float128_t NEGZEROL = -0.0L;
#else
float128_t NEGZEROL = 0.0L;
#endif
#define MAXLGM 2.556348e305

float128_t igamcl(float128_t a, float128_t x);
float128_t ndtril (float128_t y0);
float128_t igamil(float128_t a, float128_t y0);
float128_t lgaml(float128_t x);
float128_t igaml(float128_t a, float128_t x);
static float128_t stirf(float128_t);
float128_t polevll(float128_t x, void *p, int n);
float128_t p1evll(float128_t x, void *p, int n);

extern int isnanf128(float128_t x);
extern int isfinitef128(float128_t x);

float128_t imgamil(float128_t a, float128_t y0)
{
float128_t x0, x1, x, yl, yh, y, d, lgm, dithresh;
int i, dir;

/* bound the solution */
x0 = MAXNUML;
yl = 0.0L;
x1 = 0.0L;
yh = 1.0L;
dithresh = 4.0 * MACHEPL;

/* approximation to inverse function */
d = 1.0L/(9.0L*a);
y = (1.0L - d - ndtril(y0) * sqrt(d));
x = a * y * y * y;

lgm = lgaml(a);

for(i=0; i<10; i++)
	{
	if(x > x0 || x < x1)
		goto ihalve;
	y = igamcl(a,x);
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
	d = (a - 1.0L) * log(x0) - x0 - lgm;
	if(d < -MAXLOGL)
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

d = 0.0625L;
if(x0 == MAXNUML)
	{
	if(x <= 0.0L)
		x = 1.0L;
	while(x0 == MAXNUML)
		{
		x = (1.0L + d) * x;
		y = igamcl(a, x);
		if(y < y0)
			{
			x0 = x;
			yl = y;
			break;
			}
		d = d + d;
		}
	}
d = 0.5L;
dir = 0;

for(i=0; i<400; i++)
	{
	x = x1  +  d * (x0 - x1);
	y = igamcl(a, x);
	lgm = (x0 - x1)/(x1 + x0);
	if(fabs(lgm) < dithresh)
		break;
	lgm = (y - y0)/y0;
	if(fabs(lgm) < dithresh)
		break;
	if(x <= 0.0L)
		break;
	if(y > y0)
		{
		x1 = x;
		yh = y;
		if(dir < 0)
			{
			dir = 0;
			d = 0.5L;
			}
		else if(dir > 1)
			d = 0.5L * d + 0.5L;
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
			d = 0.5L;
			}
		else if(dir < -1)
			d = 0.5L * d;
		else
			d = (y0 - yl)/(yh - yl);
		dir -= 1;
		}
	}
if(x == 0.0L)
	mtherr("igamil", UNDERFLOW);

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
 * float128_t a, x, y, igaml();
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
/*							igamcl()
 *
 *	Complemented incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float128_t a, x, y, igamcl();
 *
 * y = igamcl(a, x);
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

#define BIG 9.223372036854775808e18L
#define MAXGAML 1755.455L
extern float128_t MACHEPL, MINLOGL;

float128_t igamcl(a, x)
float128_t a, x;
{
float128_t ans, c, yc, ax, y, z, r, t;
float128_t pk, pkm1, pkm2, qk, qkm1, qkm2;

if((x <= 0.0L) || (a <= 0.0L))
	return(1.0L);

if((x < 1.0L) || (x < a))
	return(1.0L - igaml(a,x));

ax = a * log(x) - x - lgaml(a);
if(ax < MINLOGL)
	{
	mtherr("igamcl", UNDERFLOW);
	return(0.0L);
	}
ax = exp(ax);

/* continued fraction */
y = 1.0L - a;
z = x + y + 1.0L;
c = 0.0L;
pkm2 = 1.0L;
qkm2 = x;
pkm1 = x + 1.0L;
qkm1 = z * x;
ans = pkm1/qkm1;

do
	{
	c += 1.0L;
	y += 1.0L;
	z += 2.0L;
	yc = y * c;
	pk = pkm1 * z  -  pkm2 * yc;
	qk = qkm1 * z  -  qkm2 * yc;
	if(qk != 0.0L)
		{
		r = pk/qk;
		t = fabs((ans - r)/r);
		ans = r;
		}
	else
		t = 1.0L;
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
while(t > MACHEPL);
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

float128_t igaml(float128_t a, float128_t x)
{
float128_t ans, ax, c, r;

if((x <= 0.0L) || (a <= 0.0L))
	return(0.0L);

if((x > 1.0L) && (x > a))
	return(1.0L - igamcl(a,x));

ax = a * log(x) - x - lgaml(a);
if(ax < MINLOGL)
	{
	mtherr("igaml", UNDERFLOW);
	return(0.0L);
	}
ax = exp(ax);

/* power series */
r = a;
c = 1.0L;
ans = 1.0L;

do
	{
	r += 1.0L;
	c *= x/r;
	ans += c;
	}
while(c/ans > MACHEPL);

return(ans * ax/a);
}

/*							ndtril.c
 *
 *	Inverse of Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * float128_t x, y, ndtril();
 *
 * x = ndtril(y);
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
 * ndtril domain      x <= 0        -MAXNUML
 * ndtril domain      x >= 1         MAXNUML
 *
 */


/*
Cephes Math Library Release 2.3:  January, 1995
Copyright 1984, 1995 by Stephen L. Moshier
*/

extern float128_t MAXNUML;

/* ndtri(y+0.5)/sqrt(2 pi) = y + y^3 R(y^2)
   0 <= y <= 3/8
   Peak relative error 6.8e-21.  */
#if UNK
/* sqrt(2pi) */
static float128_t s2pi = 2.506628274631000502416E0L;
static float128_t P0[8] = {
 8.779679420055069160496E-3L,
-7.649544967784380691785E-1L,
 2.971493676711545292135E0L,
-4.144980036933753828858E0L,
 2.765359913000830285937E0L,
-9.570456817794268907847E-1L,
 1.659219375097958322098E-1L,
-1.140013969885358273307E-2L,
};
static float128_t Q0[7] = {
/* 1.000000000000000000000E0L, */
-5.303846964603721860329E0L,
 9.908875375256718220854E0L,
-9.031318655459381388888E0L,
 4.496118508523213950686E0L,
-1.250016921424819972516E0L,
 1.823840725000038842075E-1L,
-1.088633151006419263153E-2L,
};
#endif
#if IBMPC
static unsigned short s2p[] = {
0x2cb3,0xb138,0x98ff,0xa06c,0x4000, XPD
};
#define s2pi *(float128_t *)s2p
static short P0[] = {
0xb006,0x9fc1,0xa4fe,0x8fd8,0x3ff8, XPD
0x6f8a,0x976e,0x0ed2,0xc3d4,0xbffe, XPD
0xf1f1,0x6fcc,0xf3d0,0xbe2c,0x4000, XPD
0xccfb,0xa681,0xad2c,0x84a3,0xc001, XPD
0x9a0d,0x0082,0xa825,0xb0fb,0x4000, XPD
0x13d1,0x054a,0xf220,0xf500,0xbffe, XPD
0xcee9,0x2c92,0x70bd,0xa9e7,0x3ffc, XPD
0x5fee,0x4a42,0xa6cb,0xbac7,0xbff8, XPD
};
static short Q0[] = {
/* 0x0000,0x0000,0x0000,0x8000,0x3fff, XPD */
0x841e,0xfec7,0x1d44,0xa9b9,0xc001, XPD
0x97e6,0xcde0,0xc0e7,0x9e8a,0x4002, XPD
0x66f9,0x8f3e,0x47fd,0x9080,0xc002, XPD
0x212f,0x2185,0x33ec,0x8fe0,0x4001, XPD
0x8e73,0x7bac,0x8df2,0xa000,0xbfff, XPD
0xc143,0xcb94,0xe3ea,0xbac2,0x3ffc, XPD
0x25d9,0xc8f3,0x9573,0xb25c,0xbff8, XPD
};
#endif
#if MIEEE
static unsigned long s2p[] = {
0x40000000,0xa06c98ff,0xb1382cb3,
};
#define s2pi *(float128_t *)s2p
static long P0[24] = {
0x3ff80000,0x8fd8a4fe,0x9fc1b006,
0xbffe0000,0xc3d40ed2,0x976e6f8a,
0x40000000,0xbe2cf3d0,0x6fccf1f1,
0xc0010000,0x84a3ad2c,0xa681ccfb,
0x40000000,0xb0fba825,0x00829a0d,
0xbffe0000,0xf500f220,0x054a13d1,
0x3ffc0000,0xa9e770bd,0x2c92cee9,
0xbff80000,0xbac7a6cb,0x4a425fee,
};
static long Q0[21] = {
/* 0x3fff0000,0x80000000,0x00000000, */
0xc0010000,0xa9b91d44,0xfec7841e,
0x40020000,0x9e8ac0e7,0xcde097e6,
0xc0020000,0x908047fd,0x8f3e66f9,
0x40010000,0x8fe033ec,0x2185212f,
0xbfff0000,0xa0008df2,0x7bac8e73,
0x3ffc0000,0xbac2e3ea,0xcb94c143,
0xbff80000,0xb25c9573,0xc8f325d9,
};
#endif

/* Approximation for interval z = sqrt(-2 log y) between 2 and 8
 */
/*  ndtri(p) = z - ln(z)/z - 1/z P1(1/z)/Q1(1/z)
    z = sqrt(-2 ln(p))
    2 <= z <= 8, i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
    Peak relative error 5.3e-21  */
#if UNK
static float128_t P1[10] = {
 4.302849750435552180717E0L,
 4.360209451837096682600E1L,
 9.454613328844768318162E1L,
 9.336735653151873871756E1L,
 5.305046472191852391737E1L,
 1.775851836288460008093E1L,
 3.640308340137013109859E0L,
 3.691354900171224122390E-1L,
 1.403530274998072987187E-2L,
 1.377145111380960566197E-4L,
};
static float128_t Q1[9] = {
/* 1.000000000000000000000E0L, */
 2.001425109170530136741E1L,
 7.079893963891488254284E1L,
 8.033277265194672063478E1L,
 5.034715121553662712917E1L,
 1.779820137342627204153E1L,
 3.845554944954699547539E0L,
 3.993627390181238962857E-1L,
 1.526870689522191191380E-2L,
 1.498700676286675466900E-4L,
};
#endif
#if IBMPC
static short P1[] = {
0x6105,0xb71e,0xf1f5,0x89b0,0x4001, XPD
0x461d,0x2604,0x8b77,0xae68,0x4004, XPD
0x8b33,0x4a47,0x9ec8,0xbd17,0x4005, XPD
0xa0b2,0xc1b0,0x1627,0xbabc,0x4005, XPD
0x9901,0x28f7,0xad06,0xd433,0x4004, XPD
0xddcb,0x5009,0x7213,0x8e11,0x4003, XPD
0x2432,0x0fa6,0xcfd5,0xe8fa,0x4000, XPD
0x3e24,0xd53c,0x53b2,0xbcff,0x3ffd, XPD
0x4058,0x3d75,0x5393,0xe5f4,0x3ff8, XPD
0x1789,0xf50a,0x7524,0x9067,0x3ff2, XPD
};
static short Q1[] = {
/* 0x0000,0x0000,0x0000,0x8000,0x3fff, XPD */
0xd901,0x2673,0x2fad,0xa01d,0x4003, XPD
0x24f5,0xc93c,0x0e9d,0x8d99,0x4005, XPD
0x8cda,0x523a,0x612d,0xa0aa,0x4005, XPD
0x602c,0xb5fc,0x7b9b,0xc963,0x4004, XPD
0xac72,0xd3e7,0xb766,0x8e62,0x4003, XPD
0x048e,0xe34c,0x927c,0xf61d,0x4000, XPD
0x6d88,0xa5cc,0x45de,0xcc79,0x3ffd, XPD
0xe6d1,0x199a,0x9931,0xfa29,0x3ff8, XPD
0x4c7d,0x3675,0x70a0,0x9d26,0x3ff2, XPD
};
#endif
#if MIEEE
static long P1[30] = {
0x40010000,0x89b0f1f5,0xb71e6105,
0x40040000,0xae688b77,0x2604461d,
0x40050000,0xbd179ec8,0x4a478b33,
0x40050000,0xbabc1627,0xc1b0a0b2,
0x40040000,0xd433ad06,0x28f79901,
0x40030000,0x8e117213,0x5009ddcb,
0x40000000,0xe8facfd5,0x0fa62432,
0x3ffd0000,0xbcff53b2,0xd53c3e24,
0x3ff80000,0xe5f45393,0x3d754058,
0x3ff20000,0x90677524,0xf50a1789,
};
static long Q1[27] = {
/* 0x3fff0000,0x80000000,0x00000000, */
0x40030000,0xa01d2fad,0x2673d901,
0x40050000,0x8d990e9d,0xc93c24f5,
0x40050000,0xa0aa612d,0x523a8cda,
0x40040000,0xc9637b9b,0xb5fc602c,
0x40030000,0x8e62b766,0xd3e7ac72,
0x40000000,0xf61d927c,0xe34c048e,
0x3ffd0000,0xcc7945de,0xa5cc6d88,
0x3ff80000,0xfa299931,0x199ae6d1,
0x3ff20000,0x9d2670a0,0x36754c7d,
};
#endif

/* ndtri(x) = z - ln(z)/z - 1/z P2(1/z)/Q2(1/z)
   z = sqrt(-2 ln(y))
   8 <= z <= 32
   i.e., y between exp(-32) = 1.27e-14 and exp(-512) = 4.38e-223
   Peak relative error 1.0e-21  */
#if UNK
static float128_t P2[8] = {
 3.244525725312906932464E0L,
 6.856256488128415760904E0L,
 3.765479340423144482796E0L,
 1.240893301734538935324E0L,
 1.740282292791367834724E-1L,
 9.082834200993107441750E-3L,
 1.617870121822776093899E-4L,
 7.377405643054504178605E-7L,
};
static float128_t Q2[7] = {
/* 1.000000000000000000000E0L, */
 6.021509481727510630722E0L,
 3.528463857156936773982E0L,
 1.289185315656302878699E0L,
 1.874290142615703609510E-1L,
 9.867655920899636109122E-3L,
 1.760452434084258930442E-4L,
 8.028288500688538331773E-7L,
};
#endif
#if IBMPC
static short P2[] = {
0xafb1,0x4ff9,0x4f3a,0xcfa6,0x4000, XPD
0xbd81,0xaffa,0x7401,0xdb66,0x4001, XPD
0x3a32,0x3863,0x9d0f,0xf0fd,0x4000, XPD
0x300e,0x633d,0x977a,0x9ed5,0x3fff, XPD
0xea3a,0x56b6,0x74c5,0xb234,0x3ffc, XPD
0x38c6,0x49d2,0x2af6,0x94d0,0x3ff8, XPD
0xc85d,0xe17d,0x5ed1,0xa9a5,0x3ff2, XPD
0x536c,0x808b,0x2542,0xc609,0x3fea, XPD
};
static short Q2[] = {
/* 0x0000,0x0000,0x0000,0x8000,0x3fff, XPD */
0xaabd,0x125a,0x34a7,0xc0b0,0x4001, XPD
0x0ded,0xe6da,0x5a11,0xe1d2,0x4000, XPD
0xc742,0x9d16,0x0640,0xa504,0x3fff, XPD
0xea1e,0x4cc2,0x643a,0xbfed,0x3ffc, XPD
0x7a9b,0xfaff,0xf2dd,0xa1ab,0x3ff8, XPD
0xfd90,0x4688,0xc902,0xb898,0x3ff2, XPD
0xf003,0x032a,0xfa7e,0xd781,0x3fea, XPD
};
#endif
#if MIEEE
static long P2[24] = {
0x40000000,0xcfa64f3a,0x4ff9afb1,
0x40010000,0xdb667401,0xaffabd81,
0x40000000,0xf0fd9d0f,0x38633a32,
0x3fff0000,0x9ed5977a,0x633d300e,
0x3ffc0000,0xb23474c5,0x56b6ea3a,
0x3ff80000,0x94d02af6,0x49d238c6,
0x3ff20000,0xa9a55ed1,0xe17dc85d,
0x3fea0000,0xc6092542,0x808b536c,
};
static long Q2[21] = {
/* 0x3fff0000,0x80000000,0x00000000, */
0x40010000,0xc0b034a7,0x125aaabd,
0x40000000,0xe1d25a11,0xe6da0ded,
0x3fff0000,0xa5040640,0x9d16c742,
0x3ffc0000,0xbfed643a,0x4cc2ea1e,
0x3ff80000,0xa1abf2dd,0xfaff7a9b,
0x3ff20000,0xb898c902,0x4688fd90,
0x3fea0000,0xd781fa7e,0x032af003,
};
#endif

/*  ndtri(x) = z - ln(z)/z - 1/z P3(1/z)/Q3(1/z)
    32 < z < 2048/13
    Peak relative error 1.4e-20  */
#if UNK
static float128_t P3[8] = {
 2.020331091302772535752E0L,
 2.133020661587413053144E0L,
 2.114822217898707063183E-1L,
-6.500909615246067985872E-3L,
-7.279315200737344309241E-4L,
-1.275404675610280787619E-5L,
-6.433966387613344714022E-8L,
-7.772828380948163386917E-11L,
};
static float128_t Q3[7] = {
/* 1.000000000000000000000E0L, */
 2.278210997153449199574E0L,
 2.345321838870438196534E-1L,
-6.916708899719964982855E-3L,
-7.908542088737858288849E-4L,
-1.387652389480217178984E-5L,
-7.001476867559193780666E-8L,
-8.458494263787680376729E-11L,
};
#endif
#if IBMPC
static short P3[] = {
0x87b2,0x0f31,0x1ac7,0x814d,0x4000, XPD
0x491c,0xcd74,0x6917,0x8883,0x4000, XPD
0x935e,0x1776,0xcba9,0xd88e,0x3ffc, XPD
0xbafd,0x8abb,0x9518,0xd505,0xbff7, XPD
0xc87e,0x2ed3,0xa84a,0xbed2,0xbff4, XPD
0x0094,0xa402,0x36b5,0xd5fa,0xbfee, XPD
0xbc53,0x0fc3,0x1ab2,0x8a2b,0xbfe7, XPD
0x30b4,0x71c0,0x223d,0xaaed,0xbfdd, XPD
};
static short Q3[] = {
/* 0x0000,0x0000,0x0000,0x8000,0x3fff, XPD */
0xdfc1,0x8a57,0x357f,0x91ce,0x4000, XPD
0xcc4f,0x9e03,0x346e,0xf029,0x3ffc, XPD
0x38b1,0x9788,0x8f42,0xe2a5,0xbff7, XPD
0xb281,0x2117,0x53da,0xcf51,0xbff4, XPD
0xf2ab,0x1d42,0x3760,0xe8cf,0xbfee, XPD
0x741b,0xf14f,0x06b0,0x965b,0xbfe7, XPD
0x37c2,0xa91f,0x16ea,0xba01,0xbfdd, XPD
};
#endif
#if MIEEE
static long P3[24] = {
0x40000000,0x814d1ac7,0x0f3187b2,
0x40000000,0x88836917,0xcd74491c,
0x3ffc0000,0xd88ecba9,0x1776935e,
0xbff70000,0xd5059518,0x8abbbafd,
0xbff40000,0xbed2a84a,0x2ed3c87e,
0xbfee0000,0xd5fa36b5,0xa4020094,
0xbfe70000,0x8a2b1ab2,0x0fc3bc53,
0xbfdd0000,0xaaed223d,0x71c030b4,
};
static long Q3[21] = {
/* 0x3fff0000,0x80000000,0x00000000, */
0x40000000,0x91ce357f,0x8a57dfc1,
0x3ffc0000,0xf029346e,0x9e03cc4f,
0xbff70000,0xe2a58f42,0x978838b1,
0xbff40000,0xcf5153da,0x2117b281,
0xbfee0000,0xe8cf3760,0x1d42f2ab,
0xbfe70000,0x965b06b0,0xf14f741b,
0xbfdd0000,0xba0116ea,0xa91f37c2,
};
#endif
/*
 * See the header file igamc_cephes.h for the implementation of ndtril.
 */

/*							gammaf128.c
 *
 *	Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * float128_t x, y, gammaf128();
 * extern int sgngam;
 *
 * y = gammaf128(x);
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
/*							lgaml()
 *
 *	Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * float128_t x, y, lgaml();
 * extern int sgngam;
 *
 * y = lgaml(x);
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
 * Arguments greater than MAXLGML (10^4928) return MAXNUML.
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
#if UNK
static float128_t P[8] = {
 4.212760487471622013093E-5L,
 4.542931960608009155600E-4L,
 4.092666828394035500949E-3L,
 2.385363243461108252554E-2L,
 1.113062816019361559013E-1L,
 3.629515436640239168939E-1L,
 8.378004301573126728826E-1L,
 1.000000000000000000009E0L,
};
static float128_t Q[9] = {
-1.397148517476170440917E-5L,
 2.346584059160635244282E-4L,
-1.237799246653152231188E-3L,
-7.955933682494738320586E-4L,
 2.773706565840072979165E-2L,
-4.633887671244534213831E-2L,
-2.243510905670329164562E-1L,
 4.150160950588455434583E-1L,
 9.999999999999999999908E-1L,
};
#endif
#if IBMPC
static short P[] = {
0x434a,0x3f22,0x2bda,0xb0b2,0x3ff0, XPD
0xf5aa,0xe82f,0x335b,0xee2e,0x3ff3, XPD
0xbe6c,0x3757,0xc717,0x861b,0x3ff7, XPD
0x7f43,0x5196,0xb166,0xc368,0x3ff9, XPD
0x9549,0x8eb5,0x8c3a,0xe3f4,0x3ffb, XPD
0x8d75,0x23af,0xc8e4,0xb9d4,0x3ffd, XPD
0x29cf,0x19b3,0x16c8,0xd67a,0x3ffe, XPD
0x0000,0x0000,0x0000,0x8000,0x3fff, XPD
};
static short Q[] = {
0x5473,0x2de8,0x1268,0xea67,0xbfee, XPD
0x334b,0xc2f0,0xa2dd,0xf60e,0x3ff2, XPD
0xbeed,0x1853,0xa691,0xa23d,0xbff5, XPD
0x296e,0x7cb1,0x5dfd,0xd08f,0xbff4, XPD
0x0417,0x7989,0xd7bc,0xe338,0x3ff9, XPD
0x3295,0x3698,0xd580,0xbdcd,0xbffa, XPD
0x75ef,0x3ab7,0x4ad3,0xe5bc,0xbffc, XPD
0xe458,0x2ec7,0xfd57,0xd47c,0x3ffd, XPD
0x0000,0x0000,0x0000,0x8000,0x3fff, XPD
};
#endif
#if MIEEE
static long P[24] = {
0x3ff00000,0xb0b22bda,0x3f22434a,
0x3ff30000,0xee2e335b,0xe82ff5aa,
0x3ff70000,0x861bc717,0x3757be6c,
0x3ff90000,0xc368b166,0x51967f43,
0x3ffb0000,0xe3f48c3a,0x8eb59549,
0x3ffd0000,0xb9d4c8e4,0x23af8d75,
0x3ffe0000,0xd67a16c8,0x19b329cf,
0x3fff0000,0x80000000,0x00000000,
};
static long Q[27] = {
0xbfee0000,0xea671268,0x2de85473,
0x3ff20000,0xf60ea2dd,0xc2f0334b,
0xbff50000,0xa23da691,0x1853beed,
0xbff40000,0xd08f5dfd,0x7cb1296e,
0x3ff90000,0xe338d7bc,0x79890417,
0xbffa0000,0xbdcdd580,0x36983295,
0xbffc0000,0xe5bc4ad3,0x3ab775ef,
0x3ffd0000,0xd47cfd57,0x2ec7e458,
0x3fff0000,0x80000000,0x00000000,
};
#endif
/*
static float128_t P[] = {
-3.01525602666895735709e0L,
-3.25157411956062339893e1L,
-2.92929976820724030353e2L,
-1.70730828800510297666e3L,
-7.96667499622741999770e3L,
-2.59780216007146401957e4L,
-5.99650230220855581642e4L,
-7.15743521530849602425e4L
};
static float128_t Q[] = {
 1.00000000000000000000e0L,
-1.67955233807178858919e1L,
 8.85946791747759881659e1L,
 5.69440799097468430177e1L,
-1.98526250512761318471e3L,
 3.31667508019495079814e3L,
 1.60577839621734713377e4L,
-2.97045081369399940529e4L,
-7.15743521530849602412e4L
};
*/
#define MAXGAML 1755.455L
/*static float128_t LOGPI = 1.14472988584940017414L;*/

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
static float128_t STIR[9] = {
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
#define MAXSTIR 1024.0L
static float128_t SQTPI = 2.50662827463100050242E0L;

/* 1/gamma(x) = z P(z)
 * z(x) = 1/x
 * 0 < x < 0.03125
 * Peak relative error 4.2e-23
 */
#if UNK
static float128_t S[9] = {
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
static float128_t SN[9] = {
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
extern float128_t MAXLOGL, MAXNUML, PIL;


static float128_t stirf (float128_t);
#ifdef INFINITIES
extern float128_t INFINITYL;
#endif
#ifdef NANS
extern float128_t NANL;
#endif

/* Gamma function computed by Stirling's formula.
 */
static float128_t stirf(x)
float128_t x;
{
float128_t y, w, v;

w = 1.0L/x;
/* For large x, use rational coefficients from the analytical expansion.  */
if(x > 1024.0L)
	w = (((((6.97281375836585777429E-5L * w
		+ 7.84039221720066627474E-4L) * w
		- 2.29472093621399176955E-4L) * w
		- 2.68132716049382716049E-3L) * w
		+ 3.47222222222222222222E-3L) * w
		+ 8.33333333333333333333E-2L) * w
		+ 1.0L;
else
	w = 1.0L + w * polevll(w, STIR, 8);
y = expl(x);
if(x > MAXSTIR)
	{ /* Avoid overflow in pow() */
	v = powl(x, 0.5L * x - 0.25L);
	y = v * (v / y);
	}
else
	{
	y = powl(x, x - 0.5L) / y;
	}
y = SQTPI * y * w;
return(y);
}



float128_t gammaf128(float128_t x)
{
float128_t p, q, z;
int i;

sgngaml = 1;
#ifdef NANS
if(isnanf128(x))
	return(NANL);
#endif
#ifdef INFINITIES
if(x == INFINITYL)
	return(INFINITYL);
#ifdef NANS
if(x == -INFINITYL)
	goto gamnan;
#endif
#endif
q = fabsl(x);

if(q > 13.0L)
	{
	if(q > MAXGAML)
		goto goverf;
	if(x < 0.0L)
		{
		p = floorl(q);
		if(p == q)
			{
gamnan:
#ifdef NANS
			mtherr("gammaf128", DOMAIN);
			return (NANL);
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
		z = q * sinl(PIL * z);
		z = fabsl(z) * stirf(q);
		if(z <= PIL/MAXNUML)
			{
goverf:
#ifdef INFINITIES
			return(sgngaml * INFINITYL);
#else
			mtherr("gammaf128", OVERFLOW);
			return(sgngaml * MAXNUML);
#endif
			}
		z = PIL/z;
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
p = polevll(x, P, 7);
q = polevll(x, Q, 8);
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
		q = z / (x * polevll(x, SN, 8));
		}
	else
		q = z / (x * polevll(x, S, 8));
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
static float128_t A[7] = {
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
static float128_t B[7] = {
-2.163690827643812857640E3L,
-8.723871522843511459790E4L,
-1.104326814691464261197E6L,
-6.111225012005214299996E6L,
-1.625568062543700591014E7L,
-2.003937418103815175475E7L,
-8.875666783650703802159E6L,
};
static float128_t C[7] = {
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


float128_t lgaml(float128_t x)
{
float128_t p, q, w, z, f, nx;
int i;

sgngaml = 1;
#ifdef NANS
if(isnanf128(x))
	return(NANL);
#endif
#ifdef INFINITIES
if(!isfinitef128(x))
	return(INFINITYL);
#endif
if(x < -34.0L)
	{
	q = -x;
	w = lgaml(q); /* note this modifies sgngam! */
	p = floorl(q);
	if(p == q)
		{
#ifdef INFINITIES
		mtherr("lgaml", SING);
		return (INFINITYL);
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
	z = q * sinl(PIL * z);
	if(z == 0.0L)
		goto loverf;
/*	z = LOGPI - logl(z) - w; */
	z = logl(PIL/z) - w;
	return(z);
	}

if(x < 13.0L)
	{
	z = 1.0L;
	nx = floorl(x +  0.5L);
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
		return(logl(z));
	x = (nx - 2.0L) + f;
	p = x * polevll(x, B, 6) / p1evll(x, C, 7);
	return(logl(z) + p);
	}

if(x > MAXLGM)
	{
loverf:
#ifdef INFINITIES
	return(sgngaml * INFINITYL);
#else
	mtherr("lgaml", OVERFLOW);
	return(sgngaml * MAXNUML);
#endif
	}

q = (x - 0.5L) * logl(x) - x + LS2PI;
if(x > 1.0e10L)
	return(q);
p = 1.0L/(x*x);
q += polevll(p, A, 6) / x;
return(q);


lsmall:
if(x == 0.0L)
	goto loverf;
if(x < 0.0L)
	{
	x = -x;
	q = z / (x * polevll(x, SN, 8));
	}
else
	q = z / (x * polevll(x, S, 8));
if(q < 0.0L)
	{
	sgngaml = -1;
	q = -q;
	}
else
	sgngaml = 1;
q = logl(q);
return(q);
}


float128_t ndtril(float128_t y0)
{
float128_t x, y, z, y2, x0, x1;
int code;

if(y0 <= 0.0L)
	{
	mtherr("ndtril", DOMAIN);
	return(-MAXNUML);
	}
if(y0 >= 1.0L)
	{
	mtherr("ndtri", DOMAIN);
	return(MAXNUML);
	}
code = 1;
y = y0;
if(y > (1.0L - 0.13533528323661269189L)) /* 0.135... = exp(-2) */
	{
	y = 1.0L - y;
	code = 0;
	}

if(y > 0.13533528323661269189L)
	{
	y = y - 0.5L;
	y2 = y * y;
	x = y + y * (y2 * polevll(y2, P0, 7)/p1evll(y2, Q0, 7));
	x = x * s2pi;
	return(x);
	}

x = sqrt(-2.0L * log(y));
x0 = x - log(x)/x;
z = 1.0L/x;
if(x < 8.0L)
	x1 = z * polevll(z, P1, 9)/p1evll(z, Q1, 9);
else if(x < 32.0L)
	x1 = z * polevll(z, P2, 7)/p1evll(z, Q2, 7);
else
	x1 = z * polevll(z, P3, 7)/p1evll(z, Q3, 7);
x = x0 - x1;
if(code != 0)
	x = -x;
return(x);
}
