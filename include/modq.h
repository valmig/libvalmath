
#ifndef MODQ_H
#define MODQ_H

// class for Z/qZ, q static. Not suitable for multithreading

#include <iostream>
#include <val_basics.h>

namespace val
{


namespace modqhilfsfktn
{
	int euklid(int,int,int&,int&);
	int ggT(int a1,int a2);
}

class DLL_PUBLIC modq            
{
  protected:
	  int x;
  public:
	   static int q;
	   modq(const modq&) = default;
	   modq(int a=0);      
	   operator int()const {return x;}
	   modq operator +(modq) const;
	   modq operator -() const;
	   modq operator -(modq) const;
	   modq operator *(modq) const;
	   friend DLL_PUBLIC modq inv(modq);  // computes inverse element
	   int operator |(modq) const;
	   modq operator /(modq a) const;
	   modq& operator =(modq);
	   modq& operator =(int);
	   void operator +=(modq a) {*this = *this + a;}
	   void operator -=(modq a) {*this = *this - a;}
	   void operator *=(modq a) {*this = *this * a;}
	   void operator /=(modq a) { *this = *this / a;}
	   int operator ==(modq a) const {return ( x==(a.x) );}
	   int operator !=(modq a) const {return ( x!=(a.x) );}
	   int operator ==(int a) const {return (*this == modq(a));}
	   int operator !=(int a) const {return (*this != modq(a));}
	   friend DLL_PUBLIC modq operator *(int,modq);
	   friend DLL_PUBLIC modq operator /(int,modq);
	   friend DLL_PUBLIC int inU(modq);             // Check if element is unit
	   friend DLL_PUBLIC int ord(modq);             // order of element
	   friend DLL_PUBLIC std::istream& operator >>(std::istream&,modq&);
	   friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,modq);
};

} // end namespace val
#endif // MODQ_H

