#ifndef MODQ_CPP
#define MODQ_CPP

#include <modq.h>
#include <error.h>

#define __int64 long long


namespace val
{


// ================ Hilfsfktionen =======================================================


// => gcd (a,b) = euklid (a,b,x0,y0)  ( = x0*a + y0*b )
int modqhilfsfktn::euklid(int a0,int a1,int& x0,int& y0)
{
 int h,q,x1=0,y1=1;

 x0=1;y0= 0;
 while (a1) {
       q=a0/a1;
       h=a1;
       a1=a0%a1;
       a0=h;
       h= x1;
       x1=x0 - q*x1;
       x0= h;
       h= y1;
       y1= y0 - q*y1;
       y0= h;
 }
 return a0;
}


int modqhilfsfktn::ggT(int a,int b)
{
 int vorrest=a;

 if (a<b) {
	 a=b;
	 b=vorrest;
	 vorrest=a;
 }
 while (b != 0)
 {
       vorrest=b;
       b= a % b;
       a= vorrest;
 }

 return vorrest;
}


// class modq functions :

int modq::q=1;


modq:: modq(int a)
{
   x = ((x=a%q)<0) ? (q+x) : (x); // hier gilt: 0<= x < q
}

modq modq:: operator+(modq a) const
{
 modq b;
 int m;
 m=q-x;
 if (a.x<m)
    b.x = x+a.x;
 else
    b.x=a.x-m;
 return b;
}

modq modq:: operator -() const
{
 modq b;
 b.x = q-x;
 return b;
}

modq modq:: operator -(modq a) const
{
 modq b;
 b = (*this) + (-a);
 return b;
}


modq modq::operator *(modq a) const
{
 modq b;
 __int64 h ;

 h=(__int64) (x) * (__int64) (a.x);
 b.x = int(h % (__int64) (modq::q));
 return b;
}


modq inv(modq a)
{
 int x0,y0,x1,q;
 modq b;
 q=modq::q;
 x1=modqhilfsfktn::euklid(q,a.x,x0,y0);
 if  ((x1 != 1) || (a.x==0)){
     Error::error("\nERROR modq: inv(a): a is not a unit!!");
 }
 b.x = (y0<0) ? q+y0 : y0;
 return b;
}


int modq:: operator |(modq b) const
{
 int d=modqhilfsfktn::ggT(x,q);
 return ( (b.x%d)? 0:1 );
}



modq modq:: operator /(modq a) const
{
 int x1,y1,d;
 d=modqhilfsfktn::euklid(a.x,q,x1,y1);
 if (x%d) {
    Error::error("\nERROR: modq::operator/ : division by zero!");
 }
 return modq(x1)*modq(x/d);
}


modq& modq:: operator =(modq a)
{
 x = a.x;
 return (*this);
}


modq& modq:: operator =(int a)
{
 a%=q;
 x= (a<0) ? q+a : a;
 return (*this);
}


modq operator *(int x,modq a)
{
 modq b(0);
 int y,q=modq::q;
 x %= q;
 x = (x<0) ? x+q : x;
 if (x==0) return b;
 if (x==1) return a;
 y=x/2;
 b=y*a;
 return ((x%2) ? (b+b+a) : (b+b) );
}

modq operator /(int x,modq a)
{
 return (x*inv(a));
}


int inU(modq a)
{
 int q=modq::q;
 return ( (modqhilfsfktn::ggT(a.x,q)==1) ? 1:0 );
}


int ord(modq a)
{
 int i=1;
 modq b=a;

 if (modqhilfsfktn::ggT(modq::q,b.x)!=1) return 0;

 while ((b!=1) && (b!=-1)) {
       b*=a;
       i++;
 }
 if (b==1) return i;
 else return 2*i;
}


std::istream& operator >>(std::istream& is, modq& a)
{
 int b,q=modq::q;
 is>>a.x;
 b= a.x;
 b%=q;
 if (b<0) b+=q;
 a.x=b;
 return is;
}

std::ostream& operator <<(std::ostream& os,modq a)
{
 os<<a.x;
 return os;
}

} // end namespace val

#endif  //MODQ_H
