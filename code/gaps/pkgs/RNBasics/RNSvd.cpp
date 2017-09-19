// Source file for SVD (based on code from Numerical Recipes,
// provided by Ronen Barzel and Chuck Rose


/*************************************************************
** Singular-Value Decomposition.
** 
** Had to write this by hand, since the version in Numerical
** Recipes has a bug, carried over from a bug in Golub's
** original algorithm:

        Gene H. Golub & Charles F. Van Loan
        MATRIX COMPUTATIONS
        Johns Hopkins University Press, 1983
        Fifth printing, 1987

        Page 293
        Algorithm 8.3-2: The SVD Algorithm

        ...Find the largest q and the smallest p such that if

                 +-            -+
                 | A11   0   0  |          p
                 |  0   A22  0  |        n-p-q
           A  =  |  0    0  A33 |          q
                 |  0    0   0  |         m-n
                 +-            -+

        then A33 is diagonal and A22 has a nonzero superdiagonal.
        If q = n then quit.
        If any diagonal entry in A22 is zero, then zero the superdiagonal
        element in the same row and go to Repeat...

The last sentence above is the cause of the trouble, in the following case:

                     +-   -+
                     |0 1 0|
                 A = |0 1 1|
                     |0 0 0|
                     +-   -+

In this case, q is 0, A33 is null, and A22 is the same as A.  The instruction

            "if any diagonal entry in A22 is zero, then
             zero the superdiagonal element in the same row"

cannot be applied -- A22 has a diagonal entry that is zero, but there
is no superdiagonal element in that same row.  The proper thing to
do seems to be to special-case: If A22 has 0 in its lower-right
diagonal element, zero the superdiagonal element above it.


**
** 
**
** Since the Num. Rec. code was cribbed from EISPACK or LINPACK
** or some such, which in turn cribbed from other places, I would
** consider an SVD routine suspect unless it works on the following
** test data:
**		a = 	0 1 0
**      		0 1 1
**      		0 0 0
**
** Anyway, Al Barr & I went through the references to figure out
** how it all worked.  I reimplimented it from scratch, first in
** lisp, and now in C.  Martin L. Livesey used the code and reported
** various bugs, which have been fixed.
**
**
** public routines:
**	num_svd()
**	num_svd_backsubst()
**
**	Ronen Barzel   July 1989
**       (bug fixes: Feb 1990)
**	 (ansi-C, backsubst: May 1993)
**
****************************************************************/

#include "RNBasics.h"

#define MIN(A,B)	(((A)<(B)) ? (A):(B))
#define MAX(A,B)	(((A)>(B)) ? (A):(B))
#define ALLOC2D(m,n)	(RNScalar *) malloc((unsigned)(m)*(n)*sizeof(RNScalar))
#define ALLOC1D(n)	(RNScalar *) malloc((unsigned)(n)*sizeof(RNScalar))
#define FREE(p)			free((char*)(p))
#define CLEAR2D(a,m,n)	(void) memset((char*)(a),0,(m)*(n)*sizeof(RNScalar))
#define COPY2D(d,m,n,s)	(void) memcpy((char*)(d),(char*)(s),(int)((m)*(n)*sizeof(RNScalar)))
#define REF2D(a,m,n,i,j)	(a[(n)*(i)+(j)])

RNScalar fhypot (RNScalar a, RNScalar b)
{
    return ((RNScalar) sqrt (a*a + b*b));
}

/* householder transformations
**
** looks at the submatrix of A below and to the right of the [i,j]th element
** (inclusive) (where i and j go from 0 to n-1).  Performs a householder
** transformation on the submatrix, to zero out all but the first
** elements of the first column.
**
** Matrix u (a rotation matrix) is transformed by the inverse of the householder
** transformation, so that (old u)(old a) = (new u)(new a)
**
** a is m x n,  u is m x m
**
** If the submatrix is (X1, X2, ...) the householder transformation is
** P_X1,  X1 transforms to (|X1|,0,0,0,...), and the resulting matrix
** is ((|X1|,0,0,...),P_X1 X2, P_X1 X3, ...)
*/


void householder_zero_col(RNScalar *a,RNScalar *u,int i,int j,int m,int n,RNScalar *hv)
	/* hv is a work vector, length m */
{
    RNScalar	scale,	/* a scale factor for X1 */
		sigma,		/* the norm^2 of X1 */
		snorm,		/* +/- the norm of X1 */
		f, h, s, dot;
    int k,l;

    /* we will scale X1 by its l1-norm.  squawk! */
    for (scale=0, k=i; k<m; k++) scale += ((RNScalar) fabs(REF2D(a,m,n,k,j)));

    /* if X1 is 0, no point doing anything else */
    if (!scale) return;

    /* divide out the l1-norm, and calculate sigma = |X|^2 */
    for (k=i; k<m; k++) hv[k] = REF2D(a,m,n,k,j) / scale;
    for (sigma=0, k=i; k<m; k++) sigma += hv[k]*hv[k];

    /* The householder vector is X1 +/- (|X1|,0,0,...).  We will
    ** contruct this vector in hv and use it to perform the householder
    ** transformation on the matrix.  The plus or minus on the norm
    ** is to reduce roundoff error.  */
    f = hv[i];
    snorm = ((RNScalar) ((f > 0) ? -sqrt(sigma) : sqrt(sigma)));
    h = f*snorm - sigma;	/* eqn 11.2.4 in Press et.al. */
    hv[i] = f - snorm;

    /* set the first column of a to be the (|X1|,0,...) -- this is
    ** what we would get by performing the transformation, but this
    ** way's faster */
    REF2D(a,m,n,i,j) = scale * snorm;
    for (k=i+1; k<m; k++) REF2D(a,m,n,k,j) = 0;

    /* Now perform the householder transformation on the rest of
    ** the columns.  If the householder vector is X, and -half its norm^2
    ** is h, the householder transform is P_ij : delta_ij + x_i x_j / h.
    ** We don't actually create the P matrix, we just do the calcuation.
    */
    for (k=j+1; k<n; k++) {
	for (dot=0, l=i; l<m; l++) dot += REF2D(a,m,n,l,k)*hv[l];
	s = dot/h;
	for (l=i; l<m; l++) REF2D(a,m,n,l,k) += s*hv[l];
    }

    /* Similarly, perform the householder transformation on (all)
    ** the rows of u.  Note that it's the rows of u rather than
    ** the columns, because we want U to invert what we've done to a.
    */
    for (k=0; k<m; k++) {
	for (dot=0, l=i; l<m; l++) dot += REF2D(u,m,m,k,l)*hv[l];
	s = dot/h;
	for (l=i; l<m; l++) REF2D(u,m,m,k,l) += s*hv[l];
    }
}

/* this is the same as householder_zero_col, but with rows
** and cols swapped.  
*/

void householder_zero_row(RNScalar *a,RNScalar *v,int i,int j,int m,int n,RNScalar *hv)
{
    RNScalar scale, sigma,snorm, f, h, s, dot;
    int k,l;

    for (scale=0, k=j; k<n; k++) scale += ((RNScalar) fabs(REF2D(a,m,n,i,k)));
    if (!scale) return;

    for (k=j; k<n; k++) hv[k] = REF2D(a,m,n,i,k) / scale;
    for (sigma=0, k=j; k<n; k++) sigma += hv[k]*hv[k];

    f = hv[j];
    snorm = ((RNScalar) ((f > 0) ? -sqrt(sigma) : sqrt(sigma)));
    h = f*snorm - sigma;
    hv[j] = f - snorm;

    REF2D(a,m,n,i,j) = scale * snorm;
    for (k=j+1; k<n; k++) REF2D(a,m,n,i,k) = 0;

    for (k=i+1; k<m; k++) {
	for (dot=0, l=j; l<n; l++) dot += REF2D(a,m,n,k,l)*hv[l];
	s = dot/h;
	for (l=j; l<n; l++) REF2D(a,m,n,k,l) += s*hv[l];
    }

    for (k=0; k<n; k++) {
	for (dot=0, l=j; l<n; l++) dot += REF2D(v,n,n,l,k)*hv[l];
	s = dot/h;
	for (l=j; l<n; l++) REF2D(v,n,n,l,k) += s*hv[l];
    }
}

/*
** performs a Givens rotation on the ith and jth columns of a, looking
** at the n elements starting at start.  a is mm x nn.
*/
void rotate_cols(int i,int j,RNScalar cos,RNScalar sin,RNScalar* a,int start,int n,int mm,int nn)
{
    int end = start+n, k;
    RNScalar x,y;
    for (k=start; k<end; k++) {
	x = REF2D(a,mm,nn,k,i);
	y = REF2D(a,mm,nn,k,j);
	REF2D(a,mm,nn,k,i) =  cos*x + sin*y;
	REF2D(a,mm,nn,k,j) = -sin*x + cos*y;
    }
}


/*
** performs a Givens rotation on the ith and jth rows of a, looking
** at the n elements starting at start.  a is mm x nn.
*/
void rotate_rows(int i,int j, RNScalar cos,RNScalar sin,RNScalar* a,int start,int n,int mm,int nn)
{
    int end = start+n, k;
    RNScalar x,y;
    for (k=start; k<end; k++) {
	x = REF2D(a,mm,nn,i,k);
	y = REF2D(a,mm,nn,j,k);
	REF2D(a,mm,nn,i,k) =  cos*x + sin*y;
	REF2D(a,mm,nn,j,k) = -sin*x + cos*y;
    }
}

/*
** This takes a submatrix of b (from p through z inclusive).  The submatrix
** must be bidiagonal, with a zero in the upper-left corner element.  Mucks
** with the sumatrix so that the top row is completely 0, accumulating
** the rotations into u.  b is m x n.  u is m x min.
**
** Suppose the matrix looks like
**   0 R 0 0 0...
**   0 X X 0 0...
**   0 0 X X 0...
**   ...
** Where R & X's are different & non-zero.  We can rotate the first
** and second rows, giving
**   0 0 R 0 0...
**   0 X X 0 0...
**   0 0 X X 0...
**   ...
** with new R's and X's.  We rotate the first and third rows, moving
** R over again, etc.  till R falls off the end.
*/

void clr_top_supdiag_elt(RNScalar * b,int p,int z,RNScalar * u,int m,int n,int min)
{
    int i;
    RNScalar r, x, hypot, cos, sin;

    for (i=p+1; i<=z; i++) {
	r = REF2D(b,m,n,p,i);
	x = REF2D(b,m,n,i,i);
	if (r==0) break;
	hypot = ((RNScalar) sqrt (r*r + x*x)); //hypot = pythag(r,x);
	cos = x/hypot;
	sin = r/hypot;
	/* update the diagonal and superdiagonal elements */
	REF2D(b,m,n,p,i) = 0;
	REF2D(b,m,n,i,i) = hypot;
	/* Rotate the remainder of rows p and i (only need to
	** rotate one more element, since the rest are zero) */
	if (i != z) rotate_rows(p,i,cos,sin,b,i+1,1,m,n);
	/* accumulate the rotation into u */
	rotate_cols(i,p,cos,sin,u,0,m,m,min);
    }
}

/*
** This takes a submatrix of b (from p through z inclusive).  The submatrix
** must be bidiagonal, with a zero in the lower-right corner element.  Mucks
** with the sumatrix so that the right column is completely 0, accumulating
** the rotations into v.  b is m x n.  v is min x n
**
** Suppose the matrix looks like
**   X X 0 0 
**   0 X X 0
**   0 0 X R
**   0 0 0 0
** Where R & X's are different & non-zero.  We can rotate the last
** and second-to-last columns, yielding
**   X X 0 0
**   0 X X R
**   0 0 X 0
**   0 0 0 0
** with new R's and X's.  We rotate the last and third-to-last cols, moving
** R over again, etc.  till R falls off the end.
*/

void clr_bot_supdiag_elt(RNScalar * b,int p,int z,RNScalar * v,int m,int n,int min)
{
    int i;
    RNScalar r, x, hypot, cos, sin;

    for (i=z-1; i>=p; i--) {
	r = REF2D(b,m,n,i,z);
	x = REF2D(b,m,n,i,i);
	if (r==0) break;
	hypot = ((RNScalar) sqrt (r*r + x*x)); //hypot = pythag(r,x);
	cos = x/hypot;
	sin = r/hypot;
	/* update the diagonal and superdiagonal elements */
	REF2D(b,m,n,i,z) = 0;
	REF2D(b,m,n,i,i) = hypot;
	/* Rotate the remainder of cols z and i (only need to
	** rotate one more element, since the rest are zero) */
	if (i != p) rotate_cols(i,z,cos,sin,b,i-1,1,m,n);
	/* accumulate the rotation into v */
	rotate_rows(i,z,cos,sin,v,0,n,min,n);
    }
}

/*
** This takes a submatrix of b (from p through z inclusive).  The submatrix
** must be bidiagonal except that the topmost subdiagonal element is non-zero.
** Mucks with the submatrix to make it bidiagonal, accumulating the rotations
** into u and v.  b is m x n  u is m x min   v is min x n
**
** Suppose the matrix looks like
**   X X 0 0 0...
**   R X X 0 0...
**   0 0 X X 0...
**   ...
** Where R & X's are different & non-zero.  We can rotate the first and
** second rows, giving
**   X X R 0 0...
**   0 X X 0 0...
**   0 0 X X 0...
**   ...
** with new R &X's.  Now rotate the second and third columns, getting
**   X X 0 0 0...
**   0 X X 0 0...
**   0 R X X 0...
**   ...
** which is the same as the initial problem, but decreased in
** size by 1.  Eventually, we'll reach the situation where we have
**      ...
**   ... X X
**   ... R X
** and the row rotation will eliminate R.
*/

void clr_top_subdiag_elt(RNScalar * b,int p,int z,RNScalar * u,RNScalar * v,int m,int n,int min)
{
    int i;
    RNScalar x, r, hypot, cos, sin;

    for (i=p; ; i++) {
	/* figure out row rotation to zero out the subdiagonal element */
	x = REF2D(b,m,n,i,i);
	r = REF2D(b,m,n,i+1,i);
	hypot = fhypot(x,r);
	cos = x/hypot;
	sin = r/hypot;
	/* rotate the leading elements of the row */
	REF2D(b,m,n,i,i) = hypot;
	REF2D(b,m,n,i+1,i) = 0;
	/* rotate out the rest of the row */
	rotate_rows(i,i+1,cos,sin,b,i+1,(i+1==z)?1:2,m,n);
	/* accumulate transformation into columns of u */
	rotate_cols(i,i+1,cos,sin,u,0,m,m,min);

	/* end with a row rotation */
	if (i+1==z) break;

	/* figure out column rotation */
	x = REF2D(b,m,n,i,i+1);
	r = REF2D(b,m,n,i,i+2);
	hypot = fhypot(x,r);
	cos = x/hypot;
	sin = r/hypot;
	/* rotate the leading elements of the column */
	REF2D(b,m,n,i,i+1) = hypot;
	REF2D(b,m,n,i,i+2) = 0;
	/* rotate the rest of the column */
	rotate_cols(i+1,i+2,cos,sin,b,i+1,2,m,n);
	/* accumulate transformation into columns of v */
	rotate_rows(i+1,i+2,cos,sin,v,0,n,min,n);
    }
}

/*
** This is the first part of an implicit-shift QR step.  We do some
** magic eigenvalue calculation, to calculate a Givens rotation for
** the top-left corner.
**
** This rotation is described as part of Golub's Algorithm 3.3-1
**
** This is also described as the implicit QL algorithm for symmetric
** tridiagonal matrices, in Numerical Recipes.  Here's what's going
** on, as far as I can tell:
**
** Hypothetically, one could do a series of QR decompositions of a symmetric
** tridiagonal matrix A.
**   - Ignoring the R's, one serially computes An+1 = Qn^t An Q
**   - Eventually, An will approach diagonal
**   - The rate of convergence goes like A[i,j] = u_i/u_j, where
**     u_i and u_j are the eigenvalues of A
**   - To make it converge faster, we shift the eigenvalues by some amount
**     ("uu" in the code) by decomposing the matrix A - uu I.
**   - uu is computed so as to attempt to minimize (u_i-uu)/(u_j-uu), which
**     is the convergence rate of the shifted matrix.  Ideally uu would
**     be equal to the smallest eigenvalue.  Since we don't know the
**     smallest eigenvalue, we look at the eigenvalues of the
**     trailing 2x2 submatrix.
**   - Rather than actually computing the matrix A - uu I, we just keep
**     track of the shift in the calculation of the coefficients for the
**     rotations that make up the Q's.  Hence the "implicit" in the
**     name of the algorithm.
**
** For SVD, we are looking at a bidiagonal matrix, rather than a tridiagonal
** one.  Thus we will do one more level of implicitness, by computing the
** coefficients we WOULD get were we to consider the tridiagonal matrix
** T = B^t B.
**
** This particular routine just performs the initial rotation on the
** bidiagonal matrix, based on the eigenvalue-shift stuff.  This
** makes the matrix no loger bidiagonal.  Other routines will have to
** clean up the non-bidiagonal terms.  The net rotation that is performed
** by this routine and the cleanup routines is the Q at one step of
** the above iteration.  We don't ever explicitly compute Q, we
** just keep updating the U and V matrices.
**
** b is m x n,   v is min x n
*/

void golub_kahn_svd_rot(RNScalar * b,int p,int q,RNScalar * v,int m,int n,int min)
{
    RNScalar t11, t12, t22;
    RNScalar uu;
    RNScalar b11, b12;
    RNScalar dm,fm,dn,fn;
    RNScalar hypot,cos,sin;
    RNScalar d, s, y, z;

    /* grab the last diagonal and superdiagonal elements of b22 */
    fm = (q-2) < 0 ? 0 : REF2D(b,m,n,q-2,q-1);
    dm = REF2D(b,m,n,q-1,q-1);	fn = REF2D(b,m,n,q-1,q);
				dn = REF2D(b,m,n,q,q);

    /* create the trailing 2x2 submatrix of T = b22^t b22 */
    t11 = dm*dm + fm * fm;
    t12 = dm * fn;		t22 = dn * dn + fn * fn;

    /* find the eigenvalues of the T matrix.
    **
    ** The quadratic equation for the eigenvalues has coefficients:
    ** a = 1
    ** b = -(t11+t22)
    ** c = t11 t22 - t12 t12
    **
    ** b^2 - 4ac = (t11^2 + t22^2 + 2 t11 t12) - 4 t11 t22 + 4 t12 t12
    **           = (t11 - t22)^2 + 4 t12 t12
    **
    ** sqrt(b^2-4ac) = sqrt((t11 - t22)^2 + 4 t12 t12) -- use "pythag()"
    **
    ** using quadratic formula, we have the eigenvalues:
    ** (u1,u2) = .5 (t11 + t22 +/- pythag(...))
    **         = t22 + .5 (t11 - t22 +/- pythag(...))
    **
    ** We propogate the .5 into the pythag term, yielding golub's equation 8.2-2
    ** for the "wilkinson shift".
    ** [Note:  Golub says to use (t22 + d - signum(d) s) to find the eigenvalue
    ** closer to t22.  He's right.  ]
    **/
    d = ((RNScalar) 0.5)*(t11 - t22);
    s = fhypot (d,t12);
    uu = t22 + d + ((d > 0) ? -s : s);

    /* grab the leading 2x1 of b */
    b11 = REF2D(b,m,n,p,p);	b12 = REF2D(b,m,n,p,p+1);

    /* make the leading 2x1 submatrix of T */
    t11 = b11 * b11;
    t12 = b11 * b12;

    /* calculate the rotation that would zero the off-diagonal term of the
    ** shifted matrix T - uu I
    */
    y = t11 - uu;
    z = t12;
    hypot = fhypot (y,z);
    cos = y / hypot;
    sin = z / hypot;

    /* perform the rotation on B.  This sprays some glop into the upper-left
    ** subdiagonal element.
    */
    rotate_cols(p,p+1,cos,sin,b,p,2,m,n);
    /* accumulate the rotation into the rows of v */
    rotate_rows(p,p+1,cos,sin,v,0,n,min,n);
}


/* bidiagonalize
**
** Given a (m x n)
** computes u (m x min)	orthogonal cols
**	    b (m x n)	bidiagonal
**	    v (min x n)	orthogonal rows
** such that a = u b v  (looking at only the min x min leading part of b)
**
** Works by starting with B = A, then performing a series of Householder
** transformations to zero-out the off-diagonal elements, while
** accumulating the transformations into U and V.
**
** b may point to the same array as a, in which case the result
** overwrites the original data.
*/

void bidiagonalize(const RNScalar *a,int m,int n,RNScalar *u,RNScalar *b,RNScalar *v)
{
    int i,j, min=MIN(m,n);
    RNScalar *usave = u, *vsave = v, *h = ALLOC1D(MAX(m,n));

    /* we need square matrices to accumulate the transformations */
    if (min != m) u = ALLOC2D(m,m);
    if (min != n) v = ALLOC2D(n,n);

    /* start off with u and v equal to the identity, and b equal to a */
    CLEAR2D(u,m,m);
    CLEAR2D(v,n,n);
    for (i=0; i<m; i++) REF2D(u,m,m,i,i) = 1;
    for (i=0; i<n; i++) REF2D(v,n,n,i,i) = 1;
    if (b != a) COPY2D(b,m,n,a);

    /* walk down the diagonal */
    for (i = 0; i<min; i++) {
	/* zero the entries below b[i,i] */
	householder_zero_col(b,u,i,i,m,n,h);
	/* zero the entries to the right of b[i,i+1] */
	householder_zero_row(b,v,i,i+1,m,n,h);
    }
    /* when m < n (matrix wider than tall), the above
       leaves a non-0 element in the [m-1,m]th spot.
       This is the bottom superdiagonal element of an (m+1)x(m+1);
       use the clear routine to get rid of it. */
    if (min != n)
	clr_bot_supdiag_elt(b,0,min,v,m,n,n);

    /* For non-square arrays, lop off the parts we don't need */
    if (min!=m) {
	for (i=0; i<m; i++) for (j=0; j<min; j++)
	    REF2D(usave,m,min,i,j)=REF2D(u,m,m,i,j);
	FREE(u);
    }
    if (min!=n) {
	for (i=0; i<n; i++) for (j=0; j<min; j++)
	    REF2D(vsave,min,n,j,i)=REF2D(v,n,n,j,i);
	FREE(v);
    }
    FREE(h);
}

/*
** Finds the SVD of a bidiagonal matrix.
**
** Zeroes the superdiagonal of a bidiagonal matrix, accumulating
** left- and right-hand transforms into u and v.  The matrix
** is modified in place.
**
** That is, given u, b, v  where b is bidiagonal and u & v are
** rotations, modify the matrices so that (old) u b v = (new) u b v,
** where the new b is diagonal, and u & v are still rotations.
**
** b is m x n,  u is m x min ,  and v is min x n
**
** This is Golub's SVD algorithm (Algorithm 8.3-2)
*/

void bidiagonal_svd(RNScalar *b,int m,int n,RNScalar *u,RNScalar *v)
{
    RNScalar anorm, t;
    int p, q, i, j, z, iter, min=MIN(m,n);

    /* use the max of the sum of each diagonal element and the
    ** superdiagonal element above it as a norm for the array.
    ** The use of this norm comes from the Numerical Recipes (Press et al)
    ** code.  squawk!
    */

    // CHUCKR: this is where the array access if fouling up
    //for (anorm=REF2D(b,m,n,0,0), i=1; i<n; i++) {
    for (anorm = REF2D(b,m,n,0,0), i=1 ; i < min ; i++) 
    {
	    t = ((RNScalar) (fabs(REF2D(b,m,n,i,i)) + fabs(REF2D(b,m,n,i-1,i))));

	    if (t > anorm) 
            anorm = t;
    }

    /* Each iteration of this loop will zero out the superdiagonal
    ** element above b[q][q] -- i.e. b[q][q] will be a singular value.
    */
    for (q = min-1; q>=0; q--) {

	/* Each iteration will make the superdiagonal elements smaller,
	** till they drop below numerical noise.  It ought to converge
	** by 30 tries.  squawk!  (Increased to 100 tries by funk)
	*/
        const int max_iter = 100;
	for (iter=0; iter<max_iter; iter++) {
	    
	    /* Find a block that has a zero on the diagonal or superdiagonal.
	    ** That is, we are breaking up b into three submatrices,
	    ** b11, b22, b33, such that b33 is diagonal and b22 has no zeros 
	    ** on the superdiagonal.
	    ** b33 goes from q+1 to n-1 -- it's the part we already did
	    ** b22 goes from p through q -- it's a bidiagonal block
	    */

	    /* sweep back till we reach a numerically 0 superdiagonal element */
	    for (p = q; p>0; p--) if (REF2D(b,m,n,p-1,p) + anorm == anorm) {
		REF2D(b,m,n,p-1,p) = 0;
		break;
	    }

	    /* if b22 is 1x1, i.e. there is a 0 above b[q,q], then
	    ** b[q,q] is the singular value.  Go on to the next
	    ** singular value.  (But first, make sure it's
	    ** positive.  If it's negative, negate it and the
	    ** corresponding row of v)
	    */

	    if (p==q) {
		if (REF2D(b,m,n,q,q) < 0) {
		    REF2D(b,m,n,q,q) *= -1;
		    for (j=0; j<n; j++) REF2D(v,min,n,q,j) *= -1;
		}
		break;
	    }

	    /* check for zero on the diagonal */
	    for (z= -1, i=q; i>=p; i--) if (REF2D(b,m,n,i,i)+anorm==anorm) {
		z = i;
		break;
	    }
	    if (z >= 0) {
		/* get rid of zero on the diagonl.  where is it? */
		if (z == q)
		    /* lower-right corner.  clr element above it */
		    clr_bot_supdiag_elt(b,0,z,v,m,n,min);
		else
		    /* in the middle.  clr elt to its right */
		    clr_top_supdiag_elt(b,z,q,u,m,n,min);
	    }
	    else {
		/* b22 has no zeroes on the diagonal or superdiagonal.
		** Do magic rotation, leaving glop in the uppermost
		** subdiagonal element
		*/
		golub_kahn_svd_rot(b,p,q,v,m,n,min);
		/* get rid of the glop in the uppermost subdiagonal */
		clr_top_subdiag_elt(b,p,q,u,v,m,n,min);
	    }
	}
#if 0
	if (iter>=max_iter)
	    (void) fprintf(stderr,"svd: didn't converge after %d iterations\n", max_iter);
#endif
    }
}


/**************************************************************************
**
** SVD
**
** Given a (m x n)
** computes u (m x min)	column-orthonormal
**          w (min x min)	diagonal
**          vt (min x n)	column-orthonormal
** where min=min(m,n),
** such that a = u w vt
** 
** w is returned as a vector of length min
**
**                                                t
** NOTE:  the SVD is commonly represented as U W V
**        where U is m x min and V is n x min.  This routine
**        computes the min x n transpose of V, placing the result
**        in (the memory pointed to by) parameter "vt"
**
****************************************************************************/

void RNSvdDecompose(int m, int n, 
    const RNScalar *a,
    RNScalar *u, RNScalar *w, RNScalar *vt)
{
    int i, j, k, min = MIN(m,n);
    RNScalar *g, *p, wi;

    g = ALLOC2D(m,n);
    bidiagonalize(a,m,n,u,g,vt);


    bidiagonal_svd(g,m,n,u,vt);
    for (i=0; i<min; i++) w[i] = REF2D(g,m,n,i,i);

    /* insertion sort singular values.  Note that since
    ** the svd algorithm produces values sorted or nearly sorted,
    ** an insertion sort is efficient.
    */
    for (i=1; i<min; i++) {
	if (w[i] > w[i-1]) {
	    /* save w[i], and col i of u and row i of vt.  use "g" as buffer */
	    wi = w[i];
	    p=g;
	    for (k=0; k<m; k++) *p++ = REF2D(u,m,min,k,i);
	    for (k=0; k<n; k++) *p++ = REF2D(vt,min,n,i,k);
	    /* slide columns over */
	    for (j=i; j>0 && wi > w[j-1]; j--) {
		w[j] = w[j-1];
		for (k=0; k<m; k++) REF2D(u,m,min,k,j)=REF2D(u,m,min,k,j-1);
		for (k=0; k<n; k++) REF2D(vt,min,n,j,k)=REF2D(vt,min,n,j-1,k);
	    }
	    /* put original i stuff where we ended up */
	    w[j] = wi;
	    p=g;
	    for (k=0; k<m; k++) REF2D(u,m,min,k,j) = *p++;
	    for (k=0; k<n; k++) REF2D(vt,min,n,j,k) = *p++;
	}
    }
	    
    FREE(g);
}




/**************************************************************************
**
** BACKSUBSTITUTE
**
** Given the svd of a matrix A, and a vector B, solves A X = B for
** a vector X.  Takes a condition-threshold factor "eps", singular
** values less than "eps" times the largest value are considered to be
** zero.
**
** input:
**	    u (m x min)			column-orthonormal
**          w (min-element vector)	diagonal elements
**          vt (min x n)		column-orthonormal
**    where min=min(m,n),
**    such that a = u w vt, as decomposed by num_svd()
**	    b (m-element vector)	
**
** output:
** 	    x (n-element vector)
**  
****************************************************************************/

void RNSvdBacksubstitute(int m, int n, 
    const RNScalar *u, const RNScalar *w, const RNScalar *vt, const RNScalar *b, 
    RNScalar *x, RNScalar eps)
{
    const int min = MIN(m,n);
    const RNScalar thresh = eps * w[0];	/* singular vals are sorted, w[0] is max */
    RNScalar *tmp = ALLOC1D(min);
    int i,j,k;

    for (j=0; j<min; j++) {
	RNScalar s = 0;
	if (w[j] >= thresh) {
	    for (i=0; i<m; i++) s += REF2D(u,m,min,i,j) * b[i];
	    s /= w[j];
	}
	tmp[j] = s;
    }
    for (k=0; k<n; k++) {
	RNScalar s=0;
	for (j=0; j<min; j++) s += REF2D(vt,min,n,j,k) * tmp[j];
	x[k] = s;
    }
    FREE(tmp);
}



/**************************************************************************
**
** SOLVE
**
** Given    a (m x n) 
**	    b (m-element vector)	
**
** output:
** 	    x (n-element vector)
**
** Given a matrix A, and a vector B, solves A X = B for
** a vector X.  
**
****************************************************************************/


void RNSvdSolve(int m, int n, 
    const RNScalar *a, const RNScalar *b,
    RNScalar *x, RNScalar eps)
{
  // Allocate memory
  int min = (m < n) ? m : n;
  RNScalar *u = new RNScalar [ m * min ];
  RNScalar *w = new RNScalar [ min ];
  RNScalar *vt = new RNScalar [ min * n ];
  assert(u && w && vt);

  // Decompose matrix
  RNSvdDecompose(m, n, a, u, w, vt);

  // Backsubstitute
  RNSvdBacksubstitute(m, n, u, w, vt, b, x, eps);

  // Delete memory
  delete [] u;
  delete [] w;
  delete [] vt;
}



#undef MIN
#undef MAX
#undef ALLOC2D
#undef ALLOC1D
#undef FREE
#undef CLEAR2D
#undef COPY2D
#undef REF2D
