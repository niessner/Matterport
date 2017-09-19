/* Include file for the GAPS interval class */



/* Initialization functions */

int RNInitInterval();
void RNStopInterval();



/* Class definition */

class RNInterval {
    public:
        // Constructor/destructor functions
	RNInterval();
	RNInterval(RNScalar a, RNScalar b);
	RNInterval(const RNInterval& interval);

	// Property functions
	const RNScalar Min(void) const;
	const RNScalar Max(void) const;
	const RNScalar Mid(void) const;
	const RNScalar Diameter(void) const;
	const RNScalar Radius(void) const;
	const RNBoolean Contains(RNScalar a) const;
	const RNBoolean Contains(const RNInterval& iv) const;
	const RNBoolean Inside(const RNInterval& iv) const;
	const RNBoolean Intersects(const RNInterval& iv) const;
	const RNBoolean Abuts(const RNInterval& iv) const;
	const RNBoolean Disjoint(const RNInterval& iv) const;
	const RNBoolean IsEmpty(void) const;

	// Arithmetic operators
	friend RNInterval operator +(const RNInterval &iv1, const RNInterval &interval2);
	friend RNInterval operator +(const RNInterval& iv, RNScalar a);
	friend RNInterval operator +(RNScalar a, const RNInterval& iv);
	friend RNInterval operator -(const RNInterval& iv1, const RNInterval& iv2);
	friend RNInterval operator -(const RNInterval& iv, RNScalar a);
	friend RNInterval operator -(RNScalar a, const RNInterval& iv);
	friend RNInterval operator -(const RNInterval& iv);
	friend RNInterval operator *(const RNInterval& iv1, const RNInterval& iv2);
	friend RNInterval operator *(const RNInterval& iv, RNScalar a);
	friend RNInterval operator *(RNScalar a, const RNInterval& iv);
	friend RNInterval operator /(const RNInterval& iv1, const RNInterval& iv2);
	friend RNInterval operator /(const RNInterval& iv, RNScalar a);

	// Relationship functions/operators
	friend RNBoolean operator ==(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator !=(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator  <(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator  >(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator <=(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator >=(const RNInterval& iv1, const RNInterval& iv2);
	friend RNBoolean operator  <(const RNInterval& iv, RNScalar a);
	friend RNBoolean operator  >(const RNInterval& iv, RNScalar a);
	friend RNBoolean operator <=(const RNInterval& iv, RNScalar a);
	friend RNBoolean operator >=(const RNInterval& iv, RNScalar a);
	friend RNBoolean operator  <(RNScalar a, const RNInterval& iv);
	friend RNBoolean operator  >(RNScalar a, const RNInterval& iv);
	friend RNBoolean operator <=(RNScalar a, const RNInterval& iv);
	friend RNBoolean operator >=(RNScalar a, const RNInterval& iv);

#if FALSE
	// ANSI math functions
	friend RNInterval acos(const RNInterval& iv);
	friend RNInterval asin(const RNInterval& iv);
	friend RNInterval atan(const RNInterval& iv);
	friend RNInterval atan2(const RNInterval& iv1, const RNInterval& iv2);
	friend RNInterval ceil(const RNInterval& iv);
	friend RNInterval cos(const RNInterval& iv);
	friend RNInterval cosh(const RNInterval& iv);
	friend RNInterval exp(const RNInterval& iv);
	friend RNInterval fabs(const RNInterval& iv);
	friend RNInterval floor(const RNInterval& iv);
	friend RNInterval hypot(const RNInterval& iv1, const RNInterval& iv2);
	friend RNInterval log(const RNInterval& iv);
	friend RNInterval log10(const RNInterval& iv);
	friend RNInterval pow(const RNInterval& iv, RNScalar a);
	friend RNInterval sin(const RNInterval& iv);
	friend RNInterval sinh(const RNInterval& iv);
	friend RNInterval sqr(const RNInterval& iv);
	friend RNInterval sqrt(const RNInterval& iv);
	friend RNInterval tan(const RNInterval& iv);
	friend RNInterval tanh(const RNInterval& iv);
#endif

	// Manipulation functions/operators
	void Empty(void);
	void SetMin(RNScalar a);
	void SetMax(RNScalar b);
	void Union(RNScalar a);
	void Union(const RNInterval& iv);
	void Intersect(const RNInterval& iv);
	void Reset(RNScalar a, RNScalar b);
	void operator +=(const RNInterval& iv);
	void operator +=(RNScalar a);
	void operator -=(const RNInterval& iv);
	void operator -=(RNScalar a);
	void operator *=(const RNInterval& iv);
	void operator *=(RNScalar a);
	void operator /=(const RNInterval& iv);
	void operator /=(RNScalar a);

    private:
	RNScalar lo;
	RNScalar hi;
};



/* Public variables */

extern const RNInterval RNnull_interval;
extern const RNInterval RNzero_interval;
extern const RNInterval RNunit_interval;
extern const RNInterval RNpositive_interval;
extern const RNInterval RNnegative_interval;
extern const RNInterval RNnonpositive_interval;
extern const RNInterval RNnonnegative_interval;
extern const RNInterval RNinfinite_interval;



/* Inline functions */

inline const RNScalar RNInterval::
Min(void) const
{
    return lo;
}



inline const RNScalar RNInterval::
Max(void) const
{
    return hi;
}



inline const RNScalar RNInterval::
Mid(void) const
{
    return 0.5 * (lo + hi);
}



inline const RNScalar RNInterval::
Diameter(void) const
{
    return (hi - lo);
}



inline const RNScalar RNInterval::
Radius(void) const
{
    return 0.5 * Diameter();
}



inline const RNBoolean RNInterval::
Contains(RNScalar a) const
{
    return (a >= lo) && (a <= hi);
}



inline const RNBoolean RNInterval::
Contains(const RNInterval& iv) const
{
    return (iv.lo >= lo) && (iv.hi <= hi);
}



inline const RNBoolean RNInterval::
Inside(const RNInterval& iv) const
{
    return iv.Contains(*this);
}



inline const RNBoolean RNInterval::
Intersects(const RNInterval& iv) const
{
    return (iv.lo <= hi) && (iv.hi >= lo);
}



inline const RNBoolean RNInterval::
Abuts(const RNInterval& iv) const
{
    return (iv.lo == hi) || (iv.hi == lo);
}



inline const RNBoolean RNInterval::
Disjoint(const RNInterval& iv) const
{
    return (iv.lo > hi) || (iv.hi < lo);
}



inline const RNBoolean RNInterval::
IsEmpty(void) const
{
    return (lo > hi);
}



inline RNInterval 
operator +(const RNInterval& iv1, const RNInterval& iv2) 
{
    return RNInterval(iv1.lo + iv2.lo , iv1.hi + iv2.hi);
}



inline RNInterval 
operator +(const RNInterval& iv, RNScalar a)
{
    return RNInterval(iv.lo + a, iv.hi + a);
}



inline RNInterval 
operator +(RNScalar a, const RNInterval& iv)
{
    return RNInterval(iv.lo + a, iv.hi + a);
}



inline RNInterval 
operator -(const RNInterval& iv1, const RNInterval& iv2)
{
    return RNInterval(iv1.lo - iv2.hi, iv1.hi - iv2.lo);
}



inline RNInterval 
operator -(const RNInterval& iv, RNScalar a)
{
    return RNInterval(iv.lo - a, iv.hi - a);
}



inline RNInterval 
operator -(RNScalar a, const RNInterval& iv)
{
    return RNInterval(a-iv.hi, a-iv.lo);
}



inline RNInterval 
operator -(const RNInterval& iv)
{
    return RNInterval(-iv.hi, -iv.lo);
}



inline RNInterval 
operator *(const RNInterval& iv, RNScalar a)
{
    return (a < 0) ? RNInterval(iv.hi*a, iv.lo*a) : RNInterval(iv.lo*a, iv.hi*a);
}



inline RNInterval 
operator *(RNScalar a, const RNInterval& iv)
{
    return (a < 0) ? RNInterval(iv.hi*a, iv.lo*a) : RNInterval(iv.lo*a, iv.hi*a);
}



inline RNInterval 
operator /(const RNInterval& iv, RNScalar a)
{
    return (a < 0) ? RNInterval(iv.hi/a, iv.lo/a) : RNInterval(iv.lo/a, iv.hi/a);
}



inline RNBoolean 
operator ==(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.lo == iv2.lo && iv1.hi == iv2.hi;
}



inline RNBoolean 
operator !=(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.lo != iv2.lo || iv1.hi != iv2.hi;
}



inline RNBoolean 
operator <(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.hi < iv2.lo;
}



inline RNBoolean 
operator >(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.lo > iv2.hi;
}



inline RNBoolean 
operator <=(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.hi <= iv2.hi;
}



inline RNBoolean 
operator >=(const RNInterval& iv1, const RNInterval& iv2)
{
    return iv1.lo >= iv2.lo;
}



inline RNBoolean 
operator <(const RNInterval& iv, RNScalar a)
{
    return iv.hi < a;
}



inline RNBoolean 
operator >(const RNInterval& iv, RNScalar a)
{
    return iv.lo > a;
}



inline RNBoolean 
operator <=(const RNInterval& iv, RNScalar a)
{
    return iv.lo <= a;
}



inline RNBoolean 
operator >=(const RNInterval& iv, RNScalar a)
{
    return iv.hi >= a;
}



inline RNBoolean 
operator <(RNScalar a, const RNInterval& iv)
{
    return a < iv.lo;
}



inline RNBoolean 
operator >(RNScalar a, const RNInterval& iv)
{
    return a > iv.hi;
}



inline RNBoolean 
operator <=(RNScalar a, const RNInterval& iv)
{
    return a <= iv.hi;
}



inline RNBoolean 
operator >=(RNScalar a, const RNInterval& iv)
{
    return a >= iv.lo;
}



inline void RNInterval::
Empty(void)
{
    lo = RN_INFINITY;
    hi = -RN_INFINITY;
}



inline void RNInterval::
SetMin(RNScalar a)
{
    lo = a;
}



inline void RNInterval::
SetMax(RNScalar b)
{
    hi = b;
}



inline void RNInterval::
Union(RNScalar a)
{
    if (a < lo) lo = a;
    if (a > hi) hi = a;
}



inline void RNInterval::
Union(const RNInterval& iv)
{
    if (iv.lo < lo) lo = iv.lo;
    if (iv.hi > hi) hi = iv.hi;
}



inline void RNInterval::
Intersect(const RNInterval& iv)
{
    if (iv.lo > lo) lo = iv.lo;
    if (iv.hi < hi) hi = iv.hi;
}



inline void RNInterval::
Reset(RNScalar a, RNScalar b)
{
    lo = a;
    hi = b;
}



inline void RNInterval::
operator +=(const RNInterval& iv)
{
    lo += iv.lo;
    hi += iv.hi;
}



inline void RNInterval::
operator +=(RNScalar a)
{
    lo += a;
    hi += a;
}



inline void RNInterval::
operator -=(const RNInterval& iv)
{
    lo -= iv.hi;
    hi -= iv.lo;
}



inline void RNInterval::
operator -=(RNScalar a)
{
    lo -= a;
    hi -= a;
}



inline void RNInterval::
operator *=(RNScalar a)
{
    lo *= a;
    hi *= a;
}



inline void RNInterval::
operator /=(RNScalar a)
{
    if (RNIsZero(a)) {
        lo = RN_INFINITY;
        hi = RN_INFINITY;
    }
    else {
        lo /= a;
        hi /= a;
    }
}



