/* Source file for the GAPS interval class */



/* Include files */

#include "RNBasics.h"



/* Public variables */

const RNInterval RNnull_interval(FLT_MAX, -FLT_MAX);
const RNInterval RNzero_interval(0.0, 0.0);
const RNInterval RNunit_interval(-1.0, 1.0);
const RNInterval RNpositive_interval(RN_SMALL_EPSILON, FLT_MAX);
const RNInterval RNnegative_interval(-FLT_MAX, -RN_SMALL_EPSILON);
const RNInterval RNnonpositive_interval(-FLT_MAX, 0.0);
const RNInterval RNnonnegative_interval(0.0, FLT_MAX);
const RNInterval RNinfinite_interval(-FLT_MAX, FLT_MAX);



/* Public functions */

int 
RNInitInterval()
{
    /* Return success */
    return TRUE;
}



void 
RNStopInterval()
{
}



RNInterval::
RNInterval(void)
{
}



RNInterval::
RNInterval(RNScalar a, RNScalar b)
    : lo(a), hi(b) 
{
}



RNInterval::
RNInterval(const RNInterval& interval)
    : lo(interval.lo), hi(interval.hi) 
{
}



#if FALSE

RNInterval operator *(const RNInterval& a, const RNInterval& b) const
{
	register double p, q;
	RNInterval tmp;

	if (a.lo < 0.0 && b.lo >= 0.0)
		return b * a;
	if (a.lo >= 0.0) {
		if (b.lo >= 0.0)
			return RNInterval(a.lo * b.lo, a.hi * b.hi);
		if (b.hi >= 0.0)
			return RNInterval(a.hi * b.lo, a.hi * b.hi);
		else
			return RNInterval(a.hi * b.lo, a.lo * b.hi);
	}
	if (a.hi > 0.0) {
		if (b.hi > 0.0) {
			p = a.lo * b.hi;
			q = a.hi * b.lo;
			tmp.lo = (p < q) ? p : q;
			p = a.lo * b.lo;
			q = a.hi * b.hi;
			tmp.hi = (p > q) ? p : q;
			return tmp;
		}
		return RNInterval(a.hi * b.lo, a.lo * b.lo);
	}
	if (b.hi <= 0.0)
		return RNInterval(a.hi * b.hi, a.lo * b.lo);
	else
		return RNInterval(a.lo * b.hi, a.lo * b.lo);
}

static double INF()
{
	double inf;

	inf = 0.0;
	inf = 1.0/inf;
	return inf;
}

RNInterval operator /(const RNInterval& a, const RNInterval& b)
{
	RNInterval binverse;

	if (b.lo*b.hi <= 0.0) {
		if (b.lo == 0.0 && b.hi > 0.0)
			binverse = RNInterval(1.0/b.hi, INF());
		else if (b.lo < 0.0 && b.hi == 0.0)
			binverse = RNInterval(-INF(), 1.0/b.lo);
		else
			return RNInterval(-INF(), INF());
	} else
		binverse = RNInterval(1.0/b.hi, 1.0/b.lo);
	return a * binverse;
}

void RNInterval::operator *=(RNInterval a)
{

	*this = *this * a;
}

void RNInterval::operator /=(RNInterval a)
{

	*this = *this / a;
}

RNInterval fabs(const RNInterval& a)
{

	if (a.lo >= 0.0)
		return a;
	else if (a.hi <= 0.0)
		return RNInterval(-a.hi, -a.lo);
	else
		return RNInterval(0.0, (-a.lo > a.hi) ? -a.lo : a.hi);
}

RNInterval pow(RNInterval a, double y)
{
	RNInterval aa;

	aa = fabs(a);
	return RNInterval(pow(aa.lo, y), pow(aa.hi, y));
}

RNInterval sqrt(RNInterval a)
{
	RNInterval aa;

	aa = fabs(a);
	return RNInterval(sqrt(aa.lo), sqrt(aa.hi));
}

RNInterval exp(RNInterval a)
{

	return RNInterval(exp(a.lo), exp(a.hi));
}

RNInterval log(RNInterval a)
{
	return RNInterval(log(a.lo), log(a.hi));
}

RNInterval log10(RNInterval a)
{

	return RNInterval(log10(a.lo), log10(a.hi));
}

RNInterval sqr(const RNInterval& a)
{

	if (a.lo >= 0.0)
		return RNInterval(a.lo*a.lo, a.hi*a.hi);
	else if (a.hi <= 0.0)
		return RNInterval(a.hi*a.hi, a.lo*a.lo);
	else
		return RNInterval(0.0, (-a.lo > a.hi) ? a.lo*a.lo : a.hi*a.hi);
}

RNInterval cos(RNInterval a)
{
	double loCeiling, hi;
	double coslo, coshi, mincos, maxcos;

	loCeiling = ceil(a.lo*C_PIINV);
	hi = a.hi*C_PIINV;
	if (1.0 + loCeiling <= hi)
		return RNInterval(-1.0, 1.0);
	coslo = cos(a.lo);
	coshi = cos(a.hi);
	if (coslo < coshi) {
		mincos = coslo;
		maxcos = coshi;
	} else {
		mincos = coshi;
		maxcos = coslo;
	}
	if (loCeiling <= hi) {
		if ((int)loCeiling % 2 == 0)
			return RNInterval(mincos, 1.0);
		else
			return RNInterval(-1.0, maxcos);
	} else
		return RNInterval(mincos, maxcos);
}

RNInterval sin(RNInterval a)
{

	return cos(a - 0.5*C_PI);
}

RNInterval tan(RNInterval a)
{

	if (a.hi*C_PIINV + 0.5 < ceil(a.lo*C_PIINV + 0.5))
		return RNInterval(tan(a.lo), tan(a.hi));
	else
		return RNInterval(-INF(), INF());
}

RNInterval hypot(RNInterval x, RNInterval y)
{
	double mindx, mindy, maxdx, maxdy;

	if (x.lo > 0.0) {
		mindx = x.lo;
		maxdx = x.hi;
	} else if (x.hi < 0.0) {
		mindx = x.hi;
		maxdx = x.lo;
	} else {
		mindx = 0.0;
		maxdx = (x.hi > -x.lo) ? x.hi : x.lo;
	}
	if (y.lo > 0.0) {
		mindy = y.lo;
		maxdy = y.hi;
	} else if (y.hi < 0.0) {
		mindy = y.hi;
		maxdy = y.lo;
	} else {
		mindy = 0.0;
		maxdy = (y.hi > -y.lo) ? y.hi : y.lo;
	}
	return RNInterval(hypot(mindx, mindy), hypot(maxdx, maxdy));
}

RNInterval atan2(RNInterval y, RNInterval x)
{

	if (subset(0.0, x)) {
		if (subset(0.0, y))
			return RNInterval(-C_PI, C_PI);
		else if (y > 0.0)
			return RNInterval(atan2(y.lo,x.hi), atan2(y.lo,x.lo));
		else
			return RNInterval(atan2(y.hi,x.lo), atan2(y.hi,x.hi));
	} else if (subset(0.0, y)) {
		if (x > 0.0)
			return RNInterval(atan2(y.lo,x.lo), atan2(y.hi,x.lo));
		else
			/*
			 *	RNInterval angle extends beyond PI
			 */
			return RNInterval(atan2(y.hi,x.hi),
						atan2(y.lo,x.hi) + 2.0*C_PI);
	} else if (x > 0.0) {
		if (y > 0.0)
			return RNInterval(atan2(y.lo,x.hi), atan2(y.hi,x.lo));
		else
			return RNInterval(atan2(y.lo,x.lo), atan2(y.hi,x.hi));
	} else {
		if (y > 0.0)
			return RNInterval(atan2(y.hi,x.hi), atan2(y.lo,x.lo));
		else
			return RNInterval(atan2(y.hi,x.lo), atan2(y.lo,x.hi));
	}
}

RNInterval atan(RNInterval x)
{

	return RNInterval(atan(x.lo), atan(x.hi));
}

RNInterval asin(RNInterval x)
{

	return RNInterval(asin(x.lo), asin(x.hi));
}

RNInterval acos(RNInterval x)
{

	return RNInterval(acos(x.hi), acos(x.lo));
}

RNInterval sinh(RNInterval x)
{

	return RNInterval(sinh(x.lo), sinh(x.hi));
}

RNInterval tanh(RNInterval x)
{

	return RNInterval(tanh(x.lo), tanh(x.hi));
}

RNInterval cosh(RNInterval x)
{

	if (x.lo >= 0.0)
		return RNInterval(cosh(x.lo), cosh(x.hi));
	else if (x.hi <= 0.0)
		return RNInterval(cosh(x.hi), cosh(x.lo));
	else
		return RNInterval(1.0, (-x.lo > x.hi) ? cosh(x.lo) : cosh(x.hi));
}

RNInterval ceil(RNInterval x)
{

	return RNInterval(ceil(x.lo), ceil(x.hi));
}

RNInterval floor(RNInterval x)
{

	return RNInterval(floor(x.lo), floor(x.hi));
}

#endif

