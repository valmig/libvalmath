
#include <LA.h>

namespace val
{


#ifdef INTEGER_H

template <>
void changesign(integer& wert)
{
    wert.changesign();
}

template <>
integer exactdiv(const integer& a,const integer& b)
{
    return val::EDIV(a,b);
}


#endif



#ifdef MODQ_H

template <>
int gauss(matrix<modq>& A,vector<int>& q,modq& det)
{
 int m=A.numberofrows(),n=A.numberofcolumns(),i,j,s,l,k=0,r=0;
 modq h,zero(0),one(1),inverse;

 q=vector<int>(0,n);

 if (m<=0 || n<=0) {det=zero;return 0;}

 det = one;

 for (s=0;s<n;s++) {                     // Create left triangle matrix 
     h=zero;
     for (i=r;i<m;i++)   				 // Find pivot
	 if ((h=A(i,s))!=zero) {
	    k=i;
	    break;
	 }
     if (h!=zero) {
        r++;
        q[r-1]=s;
        if (k!=(r-1)) {                 // swap rows
            A.swaprows(k,r-1);
            det=-det;
        }

        h=A(r-1,s);
        if (h!=one) {
            inverse=inv(h);
            for (j=s;j<n;j++) A(r-1,j)*=inverse;   // divide row by pivot
            det*=h;
        }

        for (i=r;i<m;i++)
            for (j=n-1;j>=s;j--) {
            A(i,j)-=A(i,s)*A(r-1,j);
        }
     }
 }

 for (i=r-1;i>=1;i--) {             // Create right triangle matrix 
     s=q[i];
     for (k=i-1;k>=0;k--) {
        for (l=s+1;l<n;l++) A(k,l)-=A(k,s)*A(i,l);
        A(k,s)=zero;
     }
 }

 if (n<m) k=n; else k=m;
 for (i=0;i<k;i++) det*=A(i,i);

 return r;
}


#endif // MODQ_H


#ifdef RATION_H

template <>
void changesign(rational& wert)
{
    wert.changesign();
}


template <>
const vector<rational>& makeprimitiv(vector<rational>& v)
{
    if (v.isempty()) return v;
    int i=0,n=v.dimension();
    rational zero(0),wert;
    integer eins(1),minuseins(-1),ggt(1),kgv(1);

    while (v(i) == zero) {
        i++;
        if (i==n) break;
    }
    if (i<n) ggt=nominator(v(i));
    i++;

    for (;i<n;i++) {
        if (v(i)==zero) continue;
        if ((ggt==eins || ggt == minuseins)) break; 
	    else {
            ggt = gcd(ggt,nominator(v(i)));
        }
    }

    i=0;
    while (denominator(v[i])==eins) {
            i++;
            if (i==n) break;
    }

    if (i<n) kgv=denominator(v[i]);
    i++;
    for (;i<n;i++) {
         if (denominator(v[i])==eins) continue;
         kgv=EDIV((kgv*denominator(v[i])),gcd(kgv,denominator(v[i])));
    }
    wert=rational(kgv,ggt);

    for (i=0;i<n;i++) v(i)*=wert;

    return v;
}


#endif // RATION_H


#ifdef COMPLEX_TYPE_H_INCLUDED


template <>
complex innerproduct(const vector<complex>& x,const vector<complex>& y,const matrix<complex>& A)
{
 int i,n=x.dimension();
 complex wert(0);


 if (x.dimension()==0 || y.dimension()==0) return wert;

 if (A.numberofrows()==0) {  // standard-product:
	 for (i=0;i<n;i++) wert += x(i).conjugate() * y(i);
 }
 else {
	 int j;
	 complex help;

	 for (i=0;i<n;i++) {
		 help = complex(0);
		 for (j=0;j<n;j++)
			 help += x(j).conjugate() * A(j,i);
		 wert += help * y(i);
	 }
 }

 return wert;
}



// Orthogonal space to  <V(0),...,V(n-1)>
template <>
vector<vector<complex> > orthogonalspace(const vector<vector<complex> >& V,const matrix<complex> &A)
{
    if (V.isempty()) return vector<vector<complex> >();

    int m=V.dimension(),n,i,j;
    complex zero(0),det;

    // n
    n=V(0).dimension();
    for (i=1;i<m;i++)
        if (V(i).dimension()>n) n = V(i).dimension();
    if (n==0) return vector<vector<complex> >();
    //
    // Create extended matrix to solve the homogeneous LES:
    matrix<complex> B(m,n+1),X;

    if (A.isempty()) {
        for (i=0;i<m;i++) {
            for (j=0;j<n;j++) B(i,j)=(V(i)(j)).conjugate();
            B(i,n) = zero;
        }
    }
    else {
        int k;
        for (i=0;i<m;i++) {
            for (j=0;j<n;j++) {
                B(i,j)=zero;
                for (k=0;k<n;k++) B(i,j)+=((V(i)(k)).conjugate())*A(k,j);
            }
            B(i,n) = zero;
        }

    }
    // Solve homogeneous LES;
    if (les(B,X,det)==0) return vector<vector<complex> >();
    m=X.numberofrows()-1;

    if (m==0) {
        vector<vector<complex> > W(vector<complex>(zero,n),1);
        return W;
    }

    vector<vector<complex> > W(vector<complex>(n),m);

    for (i=0;i<m;i++)
        for (j=0;j<n;j++) W(i)(j) = std::move(X(i+1,j));

    return W;
}


#endif // COMPLEX_TYPE_H_INCLUDED

val::matrix<double> rotationmatrix(const double& alpha,int mu,int nu,int dim)
{
	if (dim<=0 || mu<0 || mu >=dim || nu<0 || nu>=dim) return val::matrix<double>();
	val::matrix<double> A(dim);

	A.make_identity();
	if (mu==nu) return A;
	A(mu,mu) = val::cos(alpha); A(mu,nu) = -val::sin(alpha);
	A(nu,mu) = val::sin(alpha); A(nu,nu) = val::cos(alpha);
	return A;
}



} //end namespace val
