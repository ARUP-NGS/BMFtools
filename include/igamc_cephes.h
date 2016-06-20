#ifndef IGAMC_CEPHES_H
#define IGAMC_CEPHES_H

#include <inttypes.h>
#include <stdint.h>
#include <stddef.h>
#include <tgmath.h>
#include "mtherr.c"
#include "dlib/compiler_util.h"

#ifndef MAX_PV
#define MAX_PV 3117
#endif

#ifdef __cplusplus
extern "C" {
#else
#endif
extern double MAXNUM;

double igamc(double a, double x);

// Declarations needed to inline ndtri.
#ifdef UNK
/* sqrt(2pi) */
static double s2pi = 2.50662827463100050242E0;
#endif

#ifdef DEC
static unsigned short s2p[] = {0040440,0066230,0177661,0034055};
#define s2pi *(double *)s2p
#endif

#ifdef IBMPC
static unsigned short s2p[] = {0x2706,0x1ff6,0x0d93,0x4004};
#define s2pi *(double *)s2p
#endif

#ifdef MIEEE
static unsigned short s2p[] = {
0x4004,0x0d93,0x1ff6,0x2706
};
#define s2pi *(double *)s2p
#endif

/* approximation for 0 <= |y - 0.5| <= 3/8 */
#ifdef UNK
static double P0[5] = {
-5.99633501014107895267E1,
 9.80010754185999661536E1,
-5.66762857469070293439E1,
 1.39312609387279679503E1,
-1.23916583867381258016E0,
};
static double Q0[8] = {
/* 1.00000000000000000000E0,*/
 1.95448858338141759834E0,
 4.67627912898881538453E0,
 8.63602421390890590575E1,
-2.25462687854119370527E2,
 2.00260212380060660359E2,
-8.20372256168333339912E1,
 1.59056225126211695515E1,
-1.18331621121330003142E0,
};
#endif
#ifdef DEC
static unsigned short P0[20] = {
0141557,0155170,0071360,0120550,
0041704,0000214,0172417,0067307,
0141542,0132204,0040066,0156723,
0041136,0163161,0157276,0007747,
0140236,0116374,0073666,0051764,
};
static unsigned short Q0[32] = {
/*0040200,0000000,0000000,0000000,*/
0040372,0026256,0110403,0123707,
0040625,0122024,0020277,0026661,
0041654,0134161,0124134,0007244,
0142141,0073162,0133021,0131371,
0042110,0041235,0043516,0057767,
0141644,0011417,0036155,0137305,
0041176,0076556,0004043,0125430,
0140227,0073347,0152776,0067251,
};
#endif
#ifdef IBMPC
static unsigned short P0[20] = {
0x142d,0x0e5e,0xfb4f,0xc04d,
0xedd9,0x9ea1,0x8011,0x4058,
0xdbba,0x8806,0x5690,0xc04c,
0xc1fd,0x3bd7,0xdcce,0x402b,
0xca7e,0x8ef6,0xd39f,0xbff3,
};
static unsigned short Q0[36] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0x74f9,0xd220,0x4595,0x3fff,
0xe5b6,0x8417,0xb482,0x4012,
0x81d4,0x350b,0x970e,0x4055,
0x365f,0x56c2,0x2ece,0xc06c,
0xcbff,0xa8e9,0x0853,0x4069,
0xb7d9,0xe78d,0x8261,0xc054,
0x7563,0xc104,0xcfad,0x402f,
0xcdd5,0xfabf,0xeedc,0xbff2,
};
#endif
#ifdef MIEEE
static unsigned short P0[20] = {
0xc04d,0xfb4f,0x0e5e,0x142d,
0x4058,0x8011,0x9ea1,0xedd9,
0xc04c,0x5690,0x8806,0xdbba,
0x402b,0xdcce,0x3bd7,0xc1fd,
0xbff3,0xd39f,0x8ef6,0xca7e,
};
static unsigned short Q0[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x3fff,0x4595,0xd220,0x74f9,
0x4012,0xb482,0x8417,0xe5b6,
0x4055,0x970e,0x350b,0x81d4,
0xc06c,0x2ece,0x56c2,0x365f,
0x4069,0x0853,0xa8e9,0xcbff,
0xc054,0x8261,0xe78d,0xb7d9,
0x402f,0xcfad,0xc104,0x7563,
0xbff2,0xeedc,0xfabf,0xcdd5,
};
#endif


/* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
 * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
 */
#ifdef UNK
static double P1[9] = {
 4.05544892305962419923E0,
 3.15251094599893866154E1,
 5.71628192246421288162E1,
 4.40805073893200834700E1,
 1.46849561928858024014E1,
 2.18663306850790267539E0,
-1.40256079171354495875E-1,
-3.50424626827848203418E-2,
-8.57456785154685413611E-4,
};
static double Q1[8] = {
/*  1.00000000000000000000E0,*/
 1.57799883256466749731E1,
 4.53907635128879210584E1,
 4.13172038254672030440E1,
 1.50425385692907503408E1,
 2.50464946208309415979E0,
-1.42182922854787788574E-1,
-3.80806407691578277194E-2,
-9.33259480895457427372E-4,
};
#endif
#ifdef DEC
static unsigned short P1[36] = {
0040601,0143074,0150744,0073326,
0041374,0031554,0113253,0146016,
0041544,0123272,0012463,0176771,
0041460,0051160,0103560,0156511,
0041152,0172624,0117772,0030755,
0040413,0170713,0151545,0176413,
0137417,0117512,0022154,0131671,
0137017,0104257,0071432,0007072,
0135540,0143363,0063137,0036166,
};
static unsigned short Q1[32] = {
/*0040200,0000000,0000000,0000000,*/
0041174,0075325,0004736,0120326,
0041465,0110044,0047561,0045567,
0041445,0042321,0012142,0030340,
0041160,0127074,0166076,0141051,
0040440,0046055,0040745,0150400,
0137421,0114146,0067330,0010621,
0137033,0175162,0025555,0114351,
0135564,0122773,0145750,0030357,
};
#endif
#ifdef IBMPC
static unsigned short P1[36] = {
0x8edb,0x9a3c,0x38c7,0x4010,
0x7982,0x92d5,0x866d,0x403f,
0x7fbf,0x42a6,0x94d7,0x404c,
0x1ba9,0x10ee,0x0a4e,0x4046,
0x463e,0x93ff,0x5eb2,0x402d,
0xbfa1,0x7a6c,0x7e39,0x4001,
0x9677,0x448d,0xf3e9,0xbfc1,
0x41c7,0xee63,0xf115,0xbfa1,
0xe78f,0x6ccb,0x18de,0xbf4c,
};
static unsigned short Q1[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xd41b,0xa13b,0x8f5a,0x402f,
0x296f,0x89ee,0xb204,0x4046,
0x461c,0x228c,0xa89a,0x4044,
0xd845,0x9d87,0x15c7,0x402e,
0xba20,0xa83c,0x0985,0x4004,
0x0232,0xcddb,0x330c,0xbfc2,
0xb31d,0x456d,0x7f4e,0xbfa3,
0x061e,0x797d,0x94bf,0xbf4e,
};
#endif
#ifdef MIEEE
static unsigned short P1[36] = {
0x4010,0x38c7,0x9a3c,0x8edb,
0x403f,0x866d,0x92d5,0x7982,
0x404c,0x94d7,0x42a6,0x7fbf,
0x4046,0x0a4e,0x10ee,0x1ba9,
0x402d,0x5eb2,0x93ff,0x463e,
0x4001,0x7e39,0x7a6c,0xbfa1,
0xbfc1,0xf3e9,0x448d,0x9677,
0xbfa1,0xf115,0xee63,0x41c7,
0xbf4c,0x18de,0x6ccb,0xe78f,
};
static unsigned short Q1[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x402f,0x8f5a,0xa13b,0xd41b,
0x4046,0xb204,0x89ee,0x296f,
0x4044,0xa89a,0x228c,0x461c,
0x402e,0x15c7,0x9d87,0xd845,
0x4004,0x0985,0xa83c,0xba20,
0xbfc2,0x330c,0xcddb,0x0232,
0xbfa3,0x7f4e,0x456d,0xb31d,
0xbf4e,0x94bf,0x797d,0x061e,
};
#endif

/* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
 * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
 */

#ifdef UNK
static double P2[9] = {
  3.23774891776946035970E0,
  6.91522889068984211695E0,
  3.93881025292474443415E0,
  1.33303460815807542389E0,
  2.01485389549179081538E-1,
  1.23716634817820021358E-2,
  3.01581553508235416007E-4,
  2.65806974686737550832E-6,
  6.23974539184983293730E-9,
};
static double Q2[8] = {
/*  1.00000000000000000000E0,*/
  6.02427039364742014255E0,
  3.67983563856160859403E0,
  1.37702099489081330271E0,
  2.16236993594496635890E-1,
  1.34204006088543189037E-2,
  3.28014464682127739104E-4,
  2.89247864745380683936E-6,
  6.79019408009981274425E-9,
};
#endif
#ifdef DEC
static unsigned short P2[36] = {
0040517,0033507,0036236,0125641,
0040735,0044616,0014473,0140133,
0040574,0012567,0114535,0102541,
0040252,0120340,0143474,0150135,
0037516,0051057,0115361,0031211,
0036512,0131204,0101511,0125144,
0035236,0016627,0043160,0140216,
0033462,0060512,0060141,0010641,
0031326,0062541,0101304,0077706,
};
static unsigned short Q2[32] = {
/*0040200,0000000,0000000,0000000,*/
0040700,0143322,0132137,0040501,
0040553,0101155,0053221,0140257,
0040260,0041071,0052573,0010004,
0037535,0066472,0177261,0162330,
0036533,0160475,0066666,0036132,
0035253,0174533,0027771,0044027,
0033502,0016147,0117666,0063671,
0031351,0047455,0141663,0054751,
};
#endif
#ifdef IBMPC
static unsigned short P2[36] = {
0xd574,0xe793,0xe6e8,0x4009,
0x780b,0xc327,0xa931,0x401b,
0xb0ac,0xf32b,0x82ae,0x400f,
0x9a0c,0x18e7,0x541c,0x3ff5,
0x2651,0xf35e,0xca45,0x3fc9,
0x354d,0x9069,0x5650,0x3f89,
0x1812,0xe8ce,0xc3b2,0x3f33,
0x2234,0x4c0c,0x4c29,0x3ec6,
0x8ff9,0x3058,0xccac,0x3e3a,
};
static unsigned short Q2[32] = {
/*0x0000,0x0000,0x0000,0x3ff0,*/
0xe828,0x568b,0x18da,0x4018,
0x3816,0xaad2,0x704d,0x400d,
0x6200,0x2aaf,0x0847,0x3ff6,
0x3c9b,0x5fd6,0xada7,0x3fcb,
0xc78b,0xadb6,0x7c27,0x3f8b,
0x2903,0x65ff,0x7f2b,0x3f35,
0xccf7,0xf3f6,0x438c,0x3ec8,
0x6b3d,0xb876,0x29e5,0x3e3d,
};
#endif
#ifdef MIEEE
static unsigned short P2[36] = {
0x4009,0xe6e8,0xe793,0xd574,
0x401b,0xa931,0xc327,0x780b,
0x400f,0x82ae,0xf32b,0xb0ac,
0x3ff5,0x541c,0x18e7,0x9a0c,
0x3fc9,0xca45,0xf35e,0x2651,
0x3f89,0x5650,0x9069,0x354d,
0x3f33,0xc3b2,0xe8ce,0x1812,
0x3ec6,0x4c29,0x4c0c,0x2234,
0x3e3a,0xccac,0x3058,0x8ff9,
};
static unsigned short Q2[32] = {
/*0x3ff0,0x0000,0x0000,0x0000,*/
0x4018,0x18da,0x568b,0xe828,
0x400d,0x704d,0xaad2,0x3816,
0x3ff6,0x0847,0x2aaf,0x6200,
0x3fcb,0xada7,0x5fd6,0x3c9b,
0x3f8b,0x7c27,0xadb6,0xc78b,
0x3f35,0x7f2b,0x65ff,0x2903,
0x3ec8,0x438c,0xf3f6,0xccf7,
0x3e3d,0x29e5,0xb876,0x6b3d,
};
#endif


/* Polynomial evaluator:
 *  P[0] x^n  +  P[1] x^(n-1)  +  ...  +  P[n]
 */
static inline double polevl(double x, void *p, int n)
{
register double y;
register double *P = (double *)p;

y = *P++;
do
    {
    y = y * x + *P++;
    }
while( --n );
return(y);
}



/* Polynomial evaluator:
 *  x^n  +  P[0] x^(n-1)  +  P[1] x^(n-2)  +  ...  +  P[n]
 */
static inline double p1evl(double x, void *p, int n)
{
register double y;
register double *P = (double *)p;

n -= 1;
y = x + *P++;
do
    {
    y = y * x + *P++;
    }
while( --n );
return( y );
}

static inline double ndtri(double y0)
{
double x, y, z, y2, x0, x1;
int code;

if( y0 <= 0.0 )
    {
    mtherr( "ndtri", DOMAIN );
    return( -MAXNUM );
    }
if( y0 >= 1.0 )
    {
    mtherr( "ndtri", DOMAIN );
    return( MAXNUM );
    }
code = 1;
y = y0;
if( y > (1.0 - 0.13533528323661269189) ) /* 0.135... = exp(-2) */
    {
    y = 1.0 - y;
    code = 0;
    }

if( y > 0.13533528323661269189 )
    {
    y = y - 0.5;
    y2 = y * y;
    x = y + y * (y2 * polevl( y2, P0, 4)/p1evl( y2, Q0, 8 ));
    x = x * s2pi;
    return(x);
    }

x = sqrt( -2.0 * log(y) );
x0 = x - log(x)/x;

z = 1.0/x;
if( x < 8.0 ) /* y > exp(-32) = 1.2664165549e-14 */
    x1 = z * polevl( z, P1, 8 )/p1evl( z, Q1, 8 );
else
    x1 = z * polevl( z, P2, 8 )/p1evl( z, Q2, 8 );
x = x0 - x1;
if( code != 0 )
    x = -x;
return( x );
}



#ifdef UNK
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
static double LOGPI = 1.14472988584940017414;
#endif

#ifdef DEC
static unsigned short P[] = {
0035047,0162701,0146301,0005234,
0035634,0023437,0032065,0176530,
0036452,0137157,0047330,0122574,
0037103,0017310,0143041,0017232,
0037524,0066516,0162563,0164605,
0037775,0004671,0146237,0014222,
0040200,0000000,0000000,0000000
};
static unsigned short Q[] = {
0134302,0041724,0020006,0116565,
0035415,0072121,0044251,0025634,
0136222,0003447,0035205,0121114,
0036501,0107552,0154335,0104271,
0037022,0135717,0014776,0171471,
0137560,0034324,0165024,0037021,
0037222,0045046,0047151,0161213,
0040200,0000000,0000000,0000000
};
#define MAXGAM 34.84425627277176174
static unsigned short LPI[4] = {
0040222,0103202,0043475,0006750,
};
#define LOGPI *(double *)LPI
#endif

#ifdef IBMPC
static unsigned short UNUSED(P[]) = {
0x2153,0x3998,0xfcb8,0x3f24,
0xbfab,0xe686,0x84e3,0x3f53,
0x14b0,0xe9db,0x57cd,0x3f85,
0x23d3,0x18c4,0x63d9,0x3fa8,
0x7d31,0xdcae,0x8da9,0x3fca,
0xe312,0x3993,0xa137,0x3fdf,
0x0000,0x0000,0x0000,0x3ff0
};
static unsigned short UNUSED(Q[]) = {
0xd3af,0x8400,0x487a,0xbef8,
0x2573,0x2915,0xae8a,0x3f41,
0xb44a,0xe750,0x40e4,0xbf72,
0xb117,0x5b1b,0x31ed,0x3f88,
0xde67,0xe33f,0x5779,0x3fa2,
0x87c2,0x9d42,0x071a,0xbfce,
0x3c51,0xc9cd,0x4944,0x3fb2,
0x0000,0x0000,0x0000,0x3ff0
};
#define MAXGAM 171.624376956302725
static unsigned short UNUSED(LPI[4]) = {
0xa1bd,0x48e7,0x50d0,0x3ff2,
};
#define LOGPI *(double *)LPI
#endif

#ifdef MIEEE
static unsigned short P[] = {
0x3f24,0xfcb8,0x3998,0x2153,
0x3f53,0x84e3,0xe686,0xbfab,
0x3f85,0x57cd,0xe9db,0x14b0,
0x3fa8,0x63d9,0x18c4,0x23d3,
0x3fca,0x8da9,0xdcae,0x7d31,
0x3fdf,0xa137,0x3993,0xe312,
0x3ff0,0x0000,0x0000,0x0000
};
static unsigned short Q[] = {
0xbef8,0x487a,0x8400,0xd3af,
0x3f41,0xae8a,0x2915,0x2573,
0xbf72,0x40e4,0xe750,0xb44a,
0x3f88,0x31ed,0x5b1b,0xb117,
0x3fa2,0x5779,0xe33f,0xde67,
0xbfce,0x071a,0x9d42,0x87c2,
0x3fb2,0x4944,0xc9cd,0x3c51,
0x3ff0,0x0000,0x0000,0x0000
};
#define MAXGAM 171.624376956302725
static unsigned short LPI[4] = {
0x3ff2,0x50d0,0x48e7,0xa1bd,
};
#define LOGPI *(double *)LPI
#endif

/* Stirling's formula for the gamma function */
#if UNK
static double STIR[5] = {
 7.87311395793093628397E-4,
-2.29549961613378126380E-4,
-2.68132617805781232825E-3,
 3.47222221605458667310E-3,
 8.33333333333482257126E-2,
};
#define MAXSTIR 143.01608
static double SQTPI = 2.50662827463100050242E0;
#endif
#if DEC
static unsigned short STIR[20] = {
0035516,0061622,0144553,0112224,
0135160,0131531,0037460,0165740,
0136057,0134460,0037242,0077270,
0036143,0107070,0156306,0027751,
0037252,0125252,0125252,0146064,
};
#define MAXSTIR 26.77
static unsigned short SQT[4] = {
0040440,0066230,0177661,0034055,
};
#define SQTPI *(double *)SQT
#endif
#if IBMPC
static unsigned short STIR[20] = {
0x7293,0x592d,0xcc72,0x3f49,
0x1d7c,0x27e6,0x166b,0xbf2e,
0x4fd7,0x07d4,0xf726,0xbf65,
0xc5fd,0x1b98,0x71c7,0x3f6c,
0x5986,0x5555,0x5555,0x3fb5,
};
#define MAXSTIR 143.01608
static unsigned short SQT[4] = {
0x2706,0x1ff6,0x0d93,0x4004,
};
#define SQTPI *(double *)SQT
#endif
#if MIEEE
static unsigned short STIR[20] = {
0x3f49,0xcc72,0x592d,0x7293,
0xbf2e,0x166b,0x27e6,0x1d7c,
0xbf65,0xf726,0x07d4,0x4fd7,
0x3f6c,0x71c7,0x1b98,0xc5fd,
0x3fb5,0x5555,0x5555,0x5986,
};
#define MAXSTIR 143.01608
static unsigned short SQT[4] = {
0x4004,0x0d93,0x1ff6,0x2706,
};
#define SQTPI *(double *)SQT
#endif


CONST static inline double stirf(double x)
{
double y, w, v;

w = 1.0/x;
w = 1.0 + w * polevl( w, STIR, 4 );
y = exp(x);
if( x > MAXSTIR )
    { /* Avoid overflow in pow() */
    v = pow( x, 0.5 * x - 0.25 );
    y = v * (v / y);
    }
else
    {
    y = pow( x, x - 0.5 ) / y;
    }
y = SQTPI * y * w;
return( y );
}

//Multiply a phred score by this to convert a -10log_10(x) to a -2log_e(x)
#define LOG10E_X5_INV 0.460517018598809136803598290936872841520220297725754595206665580193514521935470496L
#define LOG10E_X5_1_2 0.230258509299404568401799145468436420760110148862877297603332790096757260967735248L
//such as in the following macro:
#define LOG10_TO_CHI2(x) ((x) * LOG10E_X5_INV)

#define AVG_LOG_TO_CHI2(x, y) ((x + y) * LOG10E_X5_1_2)


CONST static inline uint32_t pvalue_to_phred(double pvalue)
{
    return pvalue > 0 ? (uint32_t)(-10 * log10(pvalue) + 0.5): MAX_PV; // Add 0.5 to round up
}

// Converts a chi2 sum into a p value.
CONST static inline double igamc_pvalues(int num_pvalues, double x)
{
    return (x < 0) ? 1.0 :  igamc((double)num_pvalues, x / 2.0);
}

CONST static inline uint32_t agreed_pvalues(uint32_t pv1, uint32_t pv2)
{
    return pvalue_to_phred(igamc(2., AVG_LOG_TO_CHI2(pv1,  pv2))); // AVG divides by two while converting
}

CONST static inline uint32_t disc_pvalues(uint32_t pv_better, uint32_t pv_worse)
{
    return pv_better - pv_worse;
}


#ifdef __cplusplus
}
#else
#endif

#endif
