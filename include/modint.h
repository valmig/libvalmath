#ifndef MODINT_H_INCLUDED
#define MODINT_H_INCLUDED

// class for Z/qZ, q not static.

#include <iostream>
#include <val_basics.h>


namespace val
{

 class DLL_PUBLIC modint
 {
 private:
     int x;
     int q;

	 static int default_q;
	 static int euklid(int,int,int&,int&);
	 static int ggT(int a1,int a2);
	 static void error_q();
	 static void error_div();

 public:
     int get_q() const {return q;}
     modint() : x(0), q(default_q) {}
     modint(int wert) :x(wert),q(default_q) {}
     modint(int a,int q1) : x(a), q(q1) {}
     modint (const modint& a) :x(a.x),q(a.q) {}
     operator int() const {return x;}
     modint operator +(const modint&) const;
     modint operator -() const {return modint(-x,q);} 
     modint operator -(const modint &a) const {return operator+(-a);}

     modint operator *(const modint&) const;
     const modint& operator =(const modint &);
     friend DLL_PUBLIC modint inv(const modint &);  // computes the inverse element
     modint operator /(const modint&) const;
     const modint& operator +=(const modint&);
     const modint& operator -=(const modint &a) {return operator+=(-a);}
     const modint& operator *=(const modint&);
     const modint& operator /=(const modint&);
     int operator ==(const modint&) const;
     int operator !=(const modint&) const;

     friend DLL_PUBLIC std::istream& operator >>(std::istream&,modint&);
     friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,const modint&);
 };

}

#endif // MODINT_H_INCLUDED
