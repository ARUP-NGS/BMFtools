#include "mtherr.c"
//#include <quadmath.h>

typedef __float128 float128_t;


/* Polynomial evaluator:
 *  P[0] x^n  +  P[1] x^(n-1)  +  ...  +  P[n]
 */
inline float128_t polevll(float128_t x, void *p, int n)
{
register float128_t y;
register float128_t *P = (float128_t *)p;

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
inline float128_t p1evll(float128_t x, void *p, int n)
{
register float128_t y;
register float128_t *P = (float128_t *)p;

n -= 1;
y = x + *P++;
do
	{
	y = y * x + *P++;
	}
while( --n );
return( y );
}

