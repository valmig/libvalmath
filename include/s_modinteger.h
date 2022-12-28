#ifndef S_MODINTEGER_H_INCLUDED
#define S_MODINTEGER_H_INCLUDED

// class for Z/qZ, where q is a static integer-type. Analog to class modq.

#include <integer.h>

namespace val
{

class DLL_PUBLIC s_modinteger            // Struktur fuer Z/qZ
{
  protected:
	  integer x;
  public:
	   static integer q;
	   s_modinteger() = default;
	   s_modinteger(const s_modinteger&) = default;
	   s_modinteger(s_modinteger&& a) : x(std::move(a.x)) {}
       s_modinteger(const integer &a);        
	   operator integer() const {return x;}
	   s_modinteger operator +(const s_modinteger&) const;
	   s_modinteger operator -() const;
	   s_modinteger operator -(const s_modinteger&) const;
	   s_modinteger operator *(const s_modinteger&) const;
	   friend s_modinteger inv(const s_modinteger&);  // Berechnet inverses Elem.
	   s_modinteger operator /(const s_modinteger&) const;
	   s_modinteger& operator =(const s_modinteger &a) {x=a.x;return *this;};
	   s_modinteger& operator =(s_modinteger&& a) {x=std::move(a.x); return *this;}
	   s_modinteger& operator =(const integer& a) {x=a%q; if (x.signum()==-1) x+=q; return *this;}
	   void operator +=(const s_modinteger& a) {*this = *this + a;}
	   void operator -=(const s_modinteger& a) {*this = *this - a;}
	   void operator *=(const s_modinteger& a) {*this = *this * a;}
	   void operator /=(const s_modinteger& a) {*this = *this / a;}
	   int operator ==(const s_modinteger& a) const {return x==a.x;}
	   int operator !=(const s_modinteger& a) const {return x!=a.x;}
	   int operator ==(const integer& a) const {return *this == s_modinteger(a);}
	   int operator !=(const integer& a) const {return *this != s_modinteger(a);}
	   friend DLL_PUBLIC std::istream& operator >>(std::istream&,s_modinteger&);
	   friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,const s_modinteger&);
};


} // end namespace val

#endif // S_MODINTEGER_H_INCLUDED
