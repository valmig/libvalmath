#ifndef MULTIFLOAT_H_INCLUDED
#define MULTIFLOAT_H_INCLUDED

#include <iostream>
#include <val_basics.h>
//#include <cstdlib>


namespace val
{

class multifloat;

namespace hilfmultifloat
{
int highestbit(unsigned x); // Stelle 1<= i <= 32 des höchsten bit in x
int abs(int);
}


multifloat abs(const multifloat &x);

class multifloat
{
public:
    // Konstruktoren:
    multifloat() = default;//:exp(0),laenge(0) {mantissa=nullptr;}
    multifloat(const multifloat&);
    multifloat(multifloat&&);
    multifloat(int);
    multifloat(const double&);
    // Destruktor
    ~multifloat() {if (mantissa!=nullptr) delete[] mantissa;}
    //
    const multifloat& operator =(const multifloat&);
    const multifloat& operator =(multifloat&&);
    //
    int operator ==(const multifloat& x) const;
    int operator !=(const multifloat& x) const {return !(*this==x);}
    int operator <(const multifloat&) const;
    int operator <=(const multifloat&) const;
    int operator >(const multifloat&) const;
    int operator >=(const multifloat&) const;
    //
    const multifloat& changesign() {laenge*=-1;return *this;}
    friend multifloat abs(const multifloat &x) {multifloat y(x);if (x.laenge<0) y.laenge*=-1; return y;}
    //
    static void setprecision(int);
    static int getprecision(int) {return prec;}
    // nan, inf , -inf
    int iszero() const {return mantissa==nullptr;}
    int isNaN() const;          // laenge=0, mantissa[0]=0;
    int isInfPos() const;       // laenge=1, mantissa[0]=0;
    int isInfNeg() const;       // laenge=-1, mantissa[0]=0;
    //
    int signum() const;
    int abslength() const {if (laenge<0) return -laenge; else return laenge;}
    int lowestexp() const {return exp;}
    int highestexp() const;
    //
    multifloat round(int precisionbits) const;    // x = *round(*zhis) rundet, so dass x.mantissa höchstens precisionbits bits  enthält
    //
    operator double() const;
    //
    friend std::istream& operator >>(std::istream&,multifloat&);
    friend std::ostream& operator <<(std::ostream&,const multifloat&);
    void write() const;

private:
    unsigned *mantissa=nullptr;
    int exp=0;
    int laenge=0;
    static int prec;
    static const int numberofbits;
    static multifloat add(const multifloat& a,const multifloat& b);  // |a| + |b|
    static multifloat sub(const multifloat& a,const multifloat& b);  // |a| - |b| |a| >= |b|
};







}



#endif // MULTIFLOAT_H_INCLUDED
