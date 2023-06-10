#ifndef MODINT_CPP_INCLUDED
#define MODINT_CPP_INCLUDED

#include <modint.h>
#include <error.h>

#define __int64 long long
namespace val
{

 int modint::default_q=32003;

 // => gcd (a,b) = euklid (a,b,x0,y0)  ( = x0*a + y0*b )
 int modint::euklid(int a0,int a1,int& x0,int& y0)
 {
     int h,quot,x1=0,y1=1;
     x0=1;y0= 0;
     while (a1) {
        quot=a0/a1;
        h=a1;
        a1=a0%a1;
        a0=h;
        h= x1;
        x1=x0 - quot*x1;
        x0= h;
        h= y1;
        y1= y0 - quot*y1;
        y0= h;
    }
    return a0;
 }

 int modint::ggT(int a,int b)
 {
     int vorrest=a;

     if (a<b) {
        a=b;
        b=vorrest;
        vorrest=a;
	 }
     while (b != 0) {
       vorrest=b;
       b= a % b;
       a= vorrest;
     }

     return vorrest;
 }


 void modint::error_q()
 {
     Error::error("\nError: modint:: different q's!");
 }


 void modint::error_div()
 {
     Error::error("\nERROR: modint:: element is not a unit!");
 }


 modint modint:: operator +(const modint &a) const
 {
     if (!x) return a;
     if (!a.x) return *this;
     if (q!=a.q) error_q();
     
     if (x<0 && a.x>0) return modint(x+a.x,q);
     if (x>0 && a.x<0) return modint(x+a.x,q);
     if (x<0 && a.x<0) return modint(x+q+a.x,q);
     return modint(x-q+a.x,q);
 }


 modint modint::operator *(const modint &a) const
 {
     if (!x) return modint(0,a.q);
     if (!a.x) return modint(0,q);
     if (a.x==1) return *this;
     if (x==1) return a;
     __int64 h;

     h=(__int64) (x) * (__int64) (a.x);
     return modint(int(h % (__int64) (a.q)),a.q);
 }

const modint& modint::operator =(const modint &a)
 {
     q=a.q;x=a.x;
     return *this;
 }

 modint inv(const modint &a)
 {
    int x0,y0,x1;

    x1=modint::euklid(a.q,a.x,x0,y0);
    if  ((x1 != 1) || (a.x==0)) modint::error_div();

    return modint(y0,a.q);
 }

 modint modint::operator /(const modint &a) const
 {
     if (!x) return modint(0,a.q);
     if (a.x==1) return *this;
     int x1,y1,d;
     d=euklid(a.x,a.q,x1,y1);
     if (x%d) error_div();
     return modint(x1,a.q)*modint(x/d,a.q);

 }

 const modint& modint::operator +=(const modint& a)
 {
     if (!x) {x=a.x;q=a.q;}
     else if (!a.x) return *this;
     else if (q!=a.q) error_q();
     else if (x<0 && a.x>0) x+=a.x;
     else if (x>0 && a.x<0) x+=a.x;
     else if (x<0 && a.x<0) x+=q+a;
     else x+=a.x-q;     

     return *this;
 }


 const modint& modint::operator *=(const modint& a)
 {
     if (!x) q=a.q;
     else if (!a.x) x=0;
     else if (a.x==1) return *this;
     else if (x==1) {x=a.x;q=a.q;}
     else {
		__int64 h;

		h=(__int64) (x) * (__int64) (a.x);
		x = int(h % (__int64) (a.q)); 
		q = a.q;
	}
    return *this;
 }

 const modint& modint::operator /=(const modint& a)
 {
     *this = *this / a;
     return *this;
 }

 int modint::operator ==(const modint& a) const
 {
     if (x==a.x) return 1;
     return !((x-a.x)%q);
 }

 int modint::operator !=(const modint& a) const
 {
     return ((x-a.x)%q);
 }


 std::istream& operator >>(std::istream& is, modint& a)
 {
     int b;

     is>>a.x;
     b= a.x;
     b%=a.q;
     a.x=b;
     return is;
 }

 std::ostream& operator <<(std::ostream& os,const modint &a)
 {
     if (a.x<0) os<<a.x+a.q;
     else os<<a.x;
     return os;
 }

} //end namespace val

#endif // MODINT_CPP_INCLUDED
