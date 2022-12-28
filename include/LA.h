#ifndef LA_H_INCLUDED
#define LA_H_INCLUDED

#include <val_basics.h>
#include <matrix.h>
#include <analysis.h>
#include <vector.h>
#include <pol.h>
#include <fraction.h>


namespace val
{


template <class T>
void changesign(T& wert);

template <class T>
T exactdiv(const T& a,const T& b);

template <class T>
T euclid(T,T,T&,T&);


// Exact Gauss-algorithm. Computes also the Determinant of A.
// Returns the rank of A.
template <class T>
int gauss(matrix<T>& A,vector<int>& q,T& det);


// Gauss-Algorithm for float-types. (rounding to eps)
template <class T>
int gauss_double(matrix<T>& A,vector<int>& q,T& det,const double& eps=1e-9);


// Computes determinant of A by Gauss. Afterwards A is in diagonal form.
template <class T>
T det_by_gauss(matrix<T>& A);

template <class T>
val::fraction<val::pol<T>> det_by_gauss(val::matrix<val::fraction<val::pol<T>>>& A);

template <class T>
T (*det)(matrix<T>& A) = det_by_gauss;

// Characteristic-polynomial of A.
template <class T>
val::pol<T> char_pol(const val::matrix<T>& A);

// Computes exact solutions of a linear equation system given by the extended matrix of A.
// First row of X is a special solution, the rest of the rows form a basis for the solution of 
// the homogeneous LES. Returns the number of rows of X, e.g. 0 if LES is unsolvable. 
template <class T>
int les(matrix<T>& A,matrix<T>& X,T& det);

// Same as les but for float types. (rounding to eps)
template <class T>
int les_double(matrix<T>& A,matrix<T>& X,T& det,const double& eps=1e-9);


// i-th row = x*i-th row + y*j-th row
template <class T>
void addtorow(matrix<T> &A,int i,int j,const T& x,const T& y);


// i-th column = x*i-th column + y*j-th column
template <class T>
void addtocolumn(matrix<T> &A,int i,int j,const T& x,const T& y);

// LES for PIDs via elementary divisors. X is solution space. A is extended matrix of the LES
// Returns 1 if solvable 0 otherwise.
template <class T>
int les_integer(matrix<T> &A,matrix<T>& X);



template <class T>
vector<T> operator *(const matrix<T>& A,const vector<T>& x);


template <class T>
T innerproduct(const vector<T>& x,const vector<T>& y,const matrix<T>& A=matrix<T>());


// Standard-product.
template <class T>
T operator *(const vector<T>& x,const vector<T>& y);

template <class T>
const vector<T>& makeprimitiv(vector<T>& v) {return v;}

// Returns the vector space spanned by V via Gauss-algorithm.
template <class T>
const vector<vector<T> >& spanvspace(vector<vector<T> >& V);

// Erhard Schmidt orthogonalization  of the vectors  V(0),...,V(m).
// Returns 0 if function fails (Division by 0).
template <class T>
int orthogonalize(vector<vector<T> > &V,const matrix<T> &A=matrix<T>());


// Computes the orthogonal space of  <V(0),...,V(n-1)>
template <class T>
vector<vector<T> > orthogonalspace(const vector<vector<T> >& V,const matrix<T> &A=matrix<T>());

// Computes the rotation matrix by the angle alpha of the variables mu,nu (dim = dimension, angle in RAD)
DLL_PUBLIC val::matrix<double> rotationmatrix(const double& alpha,int mu=0,int nu=1,int dim=2); 




// ---------------------------------------------------------------------------------------------------------------------------

// friend-functions of classes  vector, matrix:
// ======================================

// Only the first r rows.
template <class T>
const vector<vector<T> >& movefrommatrixtovector(matrix<T>& A,vector<vector<T> >& V,int r)
{
    V=vector<vector<T> >();

    if (A.isempty()) return V;
    int i,m=A.numberofrows(),n=A.numberofcolumns();

    if (r<=0 || r>=m) r=m;

    V=vector<vector<T> >(r);

    for (i=0;i<r;i++) {
        V(i).coeff=A.coeff[i];
        V(i).dim=n;
        A.coeff[i]=NULL;
    }
    for (;i<m;i++) {
        delete[] A.coeff[i];
    }
    delete[] A.coeff;
    A.coeff=NULL; A.rows=A.columns=0;
    return V;
}



template <class T>
const matrix<T>& movefromvectortomatrix(vector<vector<T> >& V,matrix<T>& A)
{
    A=matrix<T>();
    if (V.isempty()) return A;

    int i,m=V.dimension(),n=V(0).dimension(),isequaldim=1,maxdim;

    maxdim=n;

    for (i=1;i<m;i++) {
        if ((n=V(i).dimension())!=maxdim) {
            isequaldim=0;
            if (n>maxdim) maxdim=n;
        }
    }
    if (isequaldim) {
        if (maxdim==0) {
            V=vector<vector<T> >();
            return A;
        }
        A.coeff= new T*[A.rows=m];
        A.columns=maxdim;
        for (i=0;i<m;i++) {
            A.coeff[i]=V(i).coeff;
            V(i).coeff=NULL;
            V(i).dim=0;
        }
        V=vector<vector<T> >();
    }
    else {
        int j,k;
        T zero(0);
        A=matrix<T>(m,maxdim);
        for (i=0;i<m;i++) {
            for (j=0;j<V(i).dim;j++) A(i,j) = std::move(V(i)(j));
            for (k=j;k<maxdim;k++) A(i,k)=zero;
        }
        V=vector<vector<T> >();
    }
    return A;
}


// ---------------------------------------------------------------------------------------------------------------------------

template <class T>
void changesign(T& wert)
{
    wert=-wert;
}


template <class T>
int gauss(matrix<T>& A,vector<int>& q,T& det)
{
 int m=A.numberofrows(),n=A.numberofcolumns(),i,j,s,l,k=0,r=0;
 T h,zero(0),one(1);

 q=vector<int>(0,n);

 if (m<=0 || n<=0) {det=zero;return 0;}

 det = one;

 for (s=0;s<n;s++) {                     // Create left triangle matrix !!
     h=zero;
     for (i=r;i<m;i++)   // Find pivot
	 if ((h=A(i,s))!=zero) {
	    k=i;
	    break;
	 }
     if (h!=zero) {
        r++;
        q[r-1]=s;
        if (k!=(r-1)) {        // swap rows
            A.swaprows(k,r-1);
            changesign(det);
        }

        h=A(r-1,s);
        if (h!=one) {
            for (j=s;j<n;j++) A(r-1,j)/=h;   // divide row by pivot
            det*=h;
        }

        for (i=r;i<m;i++)
            for (j=n-1;j>=s;j--) {
            A(i,j)-=A(i,s)*A(r-1,j);
        }
     }
 }

 for (i=r-1;i>=1;i--) {            // Create right triangle matrix
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



template <class T>
int gauss_double(matrix<T>& A,vector<int>& q,T& det,const double& eps)
{
 int m=A.numberofrows(),n=A.numberofcolumns(),i,j,s,l,k=0,r=0;
 double h,g;
 T zero(0),hi;
 vector<double> b(m);   // vector for scaling;

 q=vector<int>(0,n);
 //q.make_zero();

 if (m<=0 || n<=0) return 0;

 det = T(1);
 for (i=0;i<m;i++) {
     g=0;
     for (j=0;j<n;j++) {
        h= double(abs(A(i,j)));
        if (h>g) g= h;
     }
     b(i)= (g==0)? 1.0 : 1.0/g;
 }

 for (s=0;s<n;s++) {                     // Create left triangle matrix
     h=0;
     for (i=r;i<m;i++) {                 // Find pivot
        g= double(abs(A(i,s))) * b(i);
        if (g<eps) A(i,s)=zero;
        if (g>h) { h=g;k=i; }
     }
     if (h>=eps) {
        r++;
        q[r-1]=s;
        if (k!=(r-1)) {       // Swap rows
            swap(b(k),b(r-1));
            A.swaprows(k,r-1);
            changesign(det);
        }
        hi=A[r-1][s];
        for (j=s;j<n;j++) A(r-1,j)/=hi;   // Divide row by pivot
        det*=hi;

        for (i=r;i<m;i++)
            for (j=n-1;j>=s;j--) {
                A(i,j)-=A(i,s)*A(r-1,j);
	    }
     }
 }

 for (i=r-1;i>=1;i--) {            // Create right triangle matrix
     s=q(i);
     for (k=i-1;k>=0;k--) {
        for (l=s+1;l<n;l++) A(k,l)-=A(k,s)*A(i,l);
        A(k,s)=zero;
     }
 }
 if (n<m) k=n; else k=m;
 for (i=0;i<k;i++) det*=A(i,i);
 return r;
}


template <class T>
T det_by_gauss(matrix<T>& A)
{
 int m=A.numberofrows(),n=A.numberofcolumns(),i,j,s,k=0,r=0;
 T h,zero(0),one(1),det(1);

 if (m<=0 || n<=0) return zero;

 for (s=0;s<n;s++) {                     // Create left triangle matrix
     h=zero;
     for (i=r;i<m;i++)   // Find pivot
        if ((h=A(i,s))!=zero) {
            k=i;
            break;
        }
     if (h!=zero) {
        r++;
        if (k!=(r-1)) {        // Swap rows 
            A.swaprows(k,r-1);
            changesign(det);
        }

        h=A(r-1,s);
        if (h!=one) {
            for (j=s;j<n;j++) A(r-1,j)/=h;   // Divide row by pivot
            det*=h;
        }

        for (i=r;i<m;i++)
            for (j=n-1;j>=s;j--) {
            A(i,j)-=A(i,s)*A(r-1,j);
        }
     }
 }


 if (n<m) k=n; else k=m;
 for (i=0;i<k;i++) det*=A(i,i);
 return det;
}


template <class T>
val::fraction<val::pol<T>> det_by_gauss(val::matrix<val::fraction<val::pol<T>>>& A)
{
 using namespace val;
 int m=A.numberofrows(),n=A.numberofcolumns(),i,j,s,k=0,r=0;
 fraction<pol<T>> h,zero,one,det,minusone;
 pol<T> onepol;

 set_unity_element(onepol);
 one = onepol;

 det=one;
 minusone=-one;

 if (m<=0 || n<=0) return zero;

 for (s=0;s<n;s++) {                     // Create left triangle matrix
     h=zero;
     for (i=r;i<m;i++)   // Find pivot
        if ((h=A(i,s))!=zero) {
            k=i;
            break;
        }
     if (h!=zero) {
        r++;
        if (k!=(r-1)) {        // swap rows
            A.swaprows(k,r-1);
            det*=minusone;
        }

        h=A(r-1,s);
        if (h!=one) {
            for (j=s;j<n;j++) A(r-1,j)/=h;   // Divide row by pivot
            det*=h;
        }
        for (i=r;i<m;i++)
            for (j=n-1;j>=s;j--) {
            A(i,j)-=A(i,s)*A(r-1,j);
        }
     }
 }

 if (n<m) k=n; else k=m;
 for (i=0;i<k;i++) det*=A(i,i);
 return det;
}


template <class T>
val::pol<T> char_pol(const val::matrix<T>& A)
{
    using namespace val;
    fraction<pol<T>> determinant;
    int i,j,n=val::Min(A.numberofrows(),A.numberofcolumns());
    matrix<fraction<pol<T>>> B(n);
    pol<T> f,x(T(1),1);

    //
    for (i=0;i<n;++i) {
        for (j=0;j<n;++j) {
            if (i==j) f= pol<T>(A(i,i)) - x;
            else f= pol<T>(A(i,j));
            B(i,j)=f;
        }
    }
    determinant=det_by_gauss(B);
    if (determinant.denominator().leader()==T(-1)) return -determinant.nominator();
    else return determinant.nominator();
}



template <class T>
int les(matrix<T>& A,matrix<T>& X,T& det)
{
 int n=A.numberofcolumns(),r,i,j,s,dim;

 vector<int> q(n);

 X=val::matrix<T>();

 r=gauss(A,q,det);
 if (r>0)
    if (q[r-1]==n-1) return 0;      // LES has no solutions;

 X=matrix<T>(T(0),dim=n-r,n-1);

 vector<int> p(0,dim-1);            // not-scale-indexes
 //p.make_zero();

 for (i=j=s=0;i<n-1;i++)
     if (q(j)==i) j++;
     else { p(s)=i;s++;}
 //cout<<"\n p= "<<p;

 for (i=0;i<r;i++)                // Get X[0]
      X(0,q(i)) = A(i,n-1);

 for (i=1;i<dim;i++) X(i,p(i-1)) = T(-1);
  for (i=0;i<r;i++)
	 for (j=0;j<dim-1;j++) X(j+1,q(i))  = A(i,p(j));

 return dim;
}



template <class T>
int les_double(matrix<T>& A,matrix<T>& X,T& det,const double& eps)
{
 int n=A.numberofcolumns(),r,i,j,s,dim;

 vector<int> q(n);

 X=val::matrix<T>();

 r=gauss_double(A,q,det,eps);
 //cout<<endl<<a;
 //cout<<endl<<q;
 if (r>0)
    if (q[r-1]==n-1) return 0;      // LES has no solutions;

 X=matrix<T>(T(0),dim=n-r,n-1);

 vector<int> p(0,dim-1);            // not-scale-indexes

 for (i=j=s=0;i<n-1;i++)
     if (q(j)==i) j++;
     else { p(s)=i;s++;}
 //cout<<"\n p= "<<p;

 for (i=0;i<r;i++)                // Get X[0]
      X(0,q(i)) = A(i,n-1);

 for (i=1;i<dim;i++) X(i,p(i-1)) = T(-1);
  for (i=0;i<r;i++)
	 for (j=0;j<dim-1;j++) X(j+1,q(i))  = A(i,p(j));

 return dim;
}



template <class T>
void rowsoperation(val::matrix<T> &A,int i,int j,const T& u, const T& v, const T& x, const T& w)
{
    int n=A.numberofcolumns(),k;

    if (n==0) return;
    val::vector<T> i_row(n),j_row(n);

    for (k=0;k<n;k++) {
        i_row(k) = u*A(i,k) + v*A(j,k);
        j_row(k) = x*A(i,k) + w*A(j,k);
    }
    for (k=0;k<n;k++) {
        A(i,k) = std::move(i_row(k));
        A(j,k) = std::move(j_row(k));
    }
}

// i-column = u*i-column + v*j-column, j-row= x*i-column + w*j-column
template <class T>
void columnsoperation(val::matrix<T> &A,int i,int j,const T& u, const T& v, const T& x, const T& w)
{
    int n=A.numberofrows(),k;

    if (n==0) return;
    val::vector<T> i_column(n),j_column(n);

    for (k=0;k<n;k++) {
        i_column(k) = u*A(k,i) + v*A(k,j);
        j_column(k) = x*A(k,i) + w*A(k,j);
    }
    for (k=0;k<n;k++) {
        A(k,i) = std::move(i_column(k));
        A(k,j) = std::move(j_column(k));
    }
}



template <class T>
void addtorow(matrix<T> &A,int i,int j,const T& x,const T& y)
{
    int k,n=A.numberofcolumns();

    for (k=0;k<n;k++) {
        A(i,k) = x*A(i,k) + y*A(j,k);
    }
}


template <class T>
void addtocolumn(matrix<T> &A,int i,int j,const T& x,const T& y)
{
    int k,n=A.numberofrows();
    for (k=0;k<n;k++) A(k,i) = x*A(k,i) + y*A(k,j);
}



template <class T>
int les_integer(matrix<T> &A,matrix<T>& X)
{
    int r=0,m=A.numberofrows(),n=A.numberofcolumns()-1,actz,i,j,loesbar,ready;
    T zero(0),eins(1),x,u,v,w,d;

    if (n<=0) return 0;

    matrix<T> R(n);
    R.make_identity();

    for (actz=0;actz<m && actz<n;actz++) {
        // Find A(actz,actz) !=0,
        if (A(actz,actz) == zero) {
            // first in act-column
            for (i=actz+1;i<m;i++) if (A(i,actz) !=  zero) break;
            if (i<m) {
                A.swaprows(actz,i);
                r++;
            }
            else {
                // then in  act-column
                 for (i=actz+1;i<n;i++) if (A(actz,i)!= zero) break;
                 if (i<n) {
                    A.swapcolumns(actz,i);
                    R.swaprows(actz,i);
                    r++;
                 }
                 else continue;
            }
        }
        else r++;
        ready=0;
        while (!ready) {
            ready=1;
            for (i=actz+1;i<m;i++) {
                if (A(i,actz) != zero) {
                    ready=0;
                    if (A(i,actz)%A(actz,actz)==zero) {
                        d=exactdiv(A(i,actz),A(actz,actz));
                        changesign(d);
                        addtorow(A,i,actz,eins,d);
                    }
                    else {
                        d=euclid(A(actz,actz),A(i,actz),u,v);
                        x=exactdiv(A(actz,actz),d);
                        w=exactdiv(A(i,actz),d);
                        changesign(w);
                        rowsoperation(A,actz,i,u,v,w,x);
                    }
                }
            }
            // Create zeros in actz-row:
            for (i=actz+1;i<n;i++) {
                if (A(actz,i) != zero) {
                    ready=0;
                    if (A(actz,i)%A(actz,actz)==zero) {
                        d=exactdiv(A(actz,i),A(actz,actz));
                        changesign(d);
                        addtocolumn(A,i,actz,eins,d);
                        addtorow(R,i,actz,eins,d);
                    }
                    else {
                        d=euclid(A(actz,actz),A(actz,i),u,v);
                        x=exactdiv(A(actz,actz),d);
                        w=exactdiv(A(actz,i),d);
                        changesign(w);
                        columnsoperation(A,actz,i,u,v,w,x);
                        rowsoperation(R,actz,i,u,v,w,x);
                    }
                }
            }
        }
    }
    // diagonal-form reached.
    for (actz=0;actz<m && actz<n;actz++) {
        for (i=actz+1;i<m && i<n;i++) {
            if (A(i,i)==zero) continue;
            if (A(actz,actz)!=zero && ((A(i,i)%A(actz,actz))==zero)) continue;
            else {
                addtocolumn(A,actz,i,eins,eins);
                addtorow(R,actz,i,eins,eins);
                d=euclid(A(actz,actz),A(i,actz),u,v);
                x=exactdiv(A(actz,actz),d);
                w=exactdiv(A(i,actz),d);
                changesign(w);
                rowsoperation(A,actz,i,u,v,w,x);
                d=exactdiv(A(actz,i),A(actz,actz));
                changesign(d);
                addtocolumn(A,i,actz,eins,d);
                addtorow(R,i,actz,eins,d);
            }
        }
    }
    // Check if solvable:
    loesbar=1;
    for (i=0;i<r;i++) {
        if (A(i,n)%A(i,i)!=zero) {
            loesbar=0;
            break;
        }
    }

    if (loesbar) {
        for (;i<m;i++) {
            if (A(i,n)!=zero) {
                loesbar=0;
                break;
            }
        }
    }

    if (r==0) {
        X=matrix<T>();
        return 0;
    }

    X=matrix<T>(zero,n-r+1,n);

    // if solvable:
    if (loesbar) {
        vector<T> y(n),z(n);
        for (i=0;i<r;i++) y(i) = exactdiv(A(i,n),A(i,i));

        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) z(i) += y(j)*R(j,i);
        }

        // special solution is 1. row in X.
        for (i=0;i<n;i++) X(0,i) = std::move(z(i));
    }

    //Set the rest of the rows in X as solution space of the  hom. LES
    for (i=1;i<=n-r;i++)
        for (j=0;j<n;j++) X(i,j) = std::move(R(r+i-1,j));

    return loesbar;
}


template <class T>
vector<T> operator *(const matrix<T>& A,const vector<T>& x)
{
    vector<T> y;

    if (A.numberofrows()==0 || x.dimension()==0) return y;
    int n,m=A.numberofrows(),i,j;

    y=vector<T>(m);

    if (A.numberofcolumns()<x.dimension()) n=A.numberofcolumns();
    else n=x.dimension();
    for (i=0;i<m;i++) {
        y(i)=T(0);
        for (j=0;j<n;j++) y(i) += A(i,j)*x(j);
    }
    return y;
}


template <class T>
T innerproduct(const vector<T>& x,const vector<T>& y,const matrix<T>& A)
{
 int i,n=x.dimension();
 T wert(0);

 if (x.dimension()==0 || y.dimension()==0) return wert;

 if (A.numberofrows()==0) {  // Standardprodukt:
	 for (i=0;i<n;i++) wert += x(i) * y(i);
 }
 else {
	 int j;
	 T help;

	 for (i=0;i<n;i++) {
		 help = T(0);
		 for (j=0;j<n;j++)
			 help += x(j) * A(j,i);
		 wert += help * y(i);
	 }
 }

 return wert;
}



template <class T>
T operator *(const vector<T>& x,const vector<T>& y)
{
    return innerproduct(x,y);
}


template <class T>
const vector<vector<T> >& spanvspace(vector<vector<T> >& V)
{
    int rang;
    matrix<T> A;
    vector<int> q;
    T det;

    movefromvectortomatrix(V,A);
    rang=gauss(A,q,det);
    if (rang==0) return V;
    movefrommatrixtovector(A,V,rang);
    return V;
}


template <class T>
int orthogonalize(vector<vector<T> > &V,const matrix<T> &A)
{
 int i,j,k,m=V.dimension(),n=0;
 T zero(0);

 if (m==0 || m==1) return 1;
 
 vector<T> v(m-1),w(m-1);

 n=V(0).dimension();

 for (i=1;i<m;i++) {
	 v(i-1) = innerproduct(V(i-1),V(i-1),A);
	 if (v(i-1)==zero) return 0;
	 for (k=0;k<i;k++) {w(k)= innerproduct(V(i),V(k),A);}
	 if (V(i).dimension()!=n) return 0;
	 for (j=0;j<n;j++) {
		 for (k=0;k<i;k++) V(i)(j) -= (w(k)/v(k))*(V(k)(j));
	 }
 }
 return 1;
}


template <class T>
vector<vector<T> > orthogonalspace(const vector<vector<T> >& V,const matrix<T> &A)
{
    if (V.isempty()) return vector<vector<T> >();

    int m=V.dimension(),n,i,j;
    T zero(0),det;

    // n
    n=V(0).dimension();
    for (i=1;i<m;i++)
        if (V(i).dimension()>n) n = V(i).dimension();
    if (n==0) return vector<vector<T> >();
    //
    // Create extended matrix to solve the homog. LES:
    matrix<T> B(m,n+1),X;

    if (A.isempty()) {
        for (i=0;i<m;i++) {
            for (j=0;j<n;j++) B(i,j)=V(i)(j);
            B(i,n) = zero;
        }
    }
    else {
        int k;
        for (i=0;i<m;i++) {
            for (j=0;j<n;j++) {
                B(i,j)=zero;
                for (k=0;k<n;k++) B(i,j)+=V(i)(k)*A(k,j);
            }
            B(i,n) = zero;
        }

    }
    // Solve homogeneous LES;
    if (les(B,X,det)==0) return vector<vector<T> >();
    m=X.numberofrows()-1;

    if (m==0) {
        vector<vector<T> > W(vector<T>(zero,n),1);
        return W;
    }

    vector<vector<T> > W(vector<T>(n),m);

    for (i=0;i<m;i++)
        for (j=0;j<n;j++) W(i)(j) = std::move(X(i+1,j));

    return W;
}

} // end namespace val

#endif // LA_H_INCLUDED
