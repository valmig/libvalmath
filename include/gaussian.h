#ifndef GAUSSIAN_H_INCLUDED
#define GAUSSIAN_H_INCLUDED


#include <iostream>
#include <rational.h>


namespace val
{




class gaussian
{
   private:
	     rational re,im;
   public:
	   gaussian(): re(0),im(0) {}
	   gaussian(const rational& a) : re(a) {}
	   gaussian(rational && a) : re(std::move(a)) {}
	   gaussian(const rational& a,const rational& b) : re(a),im(b) {}
       gaussian(rational &&a,rational &&b) : re(std::move(a)), im(std::move(b)) {}
       gaussian(const gaussian &a) : re(a.re) , im(a.im) {}
       gaussian(gaussian &&a) : re(std::move(a.re)), im(std::move(a.im)) {}
	   gaussian operator +(const gaussian&) const;
	   gaussian operator -() const;
	   gaussian operator -(const gaussian&) const;
	   gaussian operator *(const gaussian&) const;
	   gaussian operator /(const gaussian&) const;
	   const gaussian& operator =(const gaussian&);
	   const gaussian& operator =(gaussian&&);
	   const gaussian& operator +=(const gaussian&);
	   const gaussian& operator -=(const gaussian&);
	   const gaussian& changesign() {re.changesign();im.changesign();return *this;}
	   const gaussian& operator *=(const gaussian&);
	   const gaussian& operator /=(const gaussian&);
	   const rational& real() const {return re;}
	   const rational& imaginary() const {return im;}
	   gaussian conjugate() const;
	   //gaussian operator ~();         // Konjugiert gaussian
	   const gaussian& operator =(const rational&);   // gaussian = rational
	   const gaussian& operator =(rational&&);
	   int iszero() const {return (re.iszero() && im.iszero());}
	   int operator ==(const gaussian&) const;
	   int operator !=(const gaussian&) const;
	   int operator ==(const rational&) const;
	   int operator !=(const rational&) const;
	   friend gaussian operator *(const rational&,const gaussian&);  // Skalarmultiplikationm
	   friend std::istream& operator >>(std::istream&,gaussian&);
	   friend std::ostream& operator <<(std::ostream&,const gaussian&);
};




} // end namespace val




#endif // GAUSSIAN_H_INCLUDED
