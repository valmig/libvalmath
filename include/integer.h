#ifndef INTEGER_H
#define INTEGER_H


#include <val_basics.h>
#include <iostream>

namespace val {

class integer;
class rational;

namespace hilfinteger {
 DLL_PUBLIC inline int abs(int n ) {return ( (n<0)? -n:n );}
 inline int max(int n,int m )  {return (n>m)? n:m;}
 inline int min(int n ,int m) { return (n<m)? n:m;}
 inline int signum(int n) {return (n<0)? -1:1;}

 void multunsigned(unsigned a,unsigned b,unsigned &high,unsigned &low);
 int multbyunsigned(unsigned*,int,unsigned,unsigned*);
 unsigned divrest64(unsigned,unsigned,unsigned);
 unsigned trialdiv(unsigned*,unsigned*);
 //unsigned *trialproduct=new unsigned[3];
 unsigned modularinvers(unsigned);
 integer exactdivision(integer,integer);
 unsigned ggT(unsigned,unsigned);
 integer restsingle(const integer&,const integer&);
/*
#ifdef STDTHREAD
 void karatsuba(const integer&,const integer&,integer&);
#endif
*/
}

// Declaration of friend functions, so they belong to namespace val 
DLL_PUBLIC integer ggTspez(const integer &a,const integer &b);
DLL_PUBLIC integer shift(const integer&,int);  // shift(x,n) = x * 2^(n*numberofbits)
DLL_PUBLIC integer trunc(const integer&,int);
DLL_PUBLIC integer rest(const integer&,int); // Regel: x = shift(trunc(x,n),n) + rest(x,n)

DLL_PUBLIC void divrem(const integer&,const integer&,integer&,integer&);

// Exact division
DLL_PUBLIC integer EDIV(const integer&,const integer&);    

DLL_PUBLIC integer gcdeuk(integer,integer);
DLL_PUBLIC integer lcm(const integer&,const integer&);
DLL_PUBLIC integer gcd(const integer &a,const integer &b);
DLL_PUBLIC integer euclid(integer,integer,integer&,integer&);
DLL_PUBLIC integer abs(const integer& a);
// -------------------------------------------------------------------------

// Binary representation of integers.
// abs(laenge) is the length of the array, sign(laenge) the sign of the integer.
class DLL_PUBLIC integer
{
   //private:
public:
       static int mem_optimized;
       // constructors and destructor
	   integer() = default;//{dat=NULL;laenge=0;}
	   integer(int);
	   integer(unsigned);
	   integer(const integer&);
	   integer(const rational&);
	   integer(integer&&);// Move-Konstruktor:
	   integer(const integer&,int); // entspr. integer(const integer&)*2^(n*numberofbits)
	   ~integer() {delete[] dat;}
       // cast-operators 
	   operator int() const;
	   operator unsigned() const;
	   operator double() const;    //nutzt ieee 754 - Eigenschaften.
       // assignment-operators
	   integer& operator =(const integer&);
	   integer& operator =(integer&&); //move-Zuweisung.
       integer& operator =(int);
	   // misc
	   int signum() const {return (laenge < 0) ? -1:1;}         // x.signum=1 <=> x>=0; x.signum=-1 <=> x<0
	   int length() const {return laenge;}
	   int abslength() const {return hilfinteger::abs(laenge);}
	   int iseven() const;
	   int isodd() const;
	   int isNull() const {return (dat==NULL);}
	   int iszero() const {return (dat==NULL);}
	   int isunit() const {return ((laenge==1 || laenge==-1) && dat[0]==1);}
	   static int GetMaxlength(){return Maxlength;}
	   unsigned operator[] (int i) const {if (i<0 || i>=hilfinteger::abs(laenge)) return 0; else return dat[i];}
       // compare-operators 
	   int operator ==(const integer&) const;
	   int operator !=(const integer&) const;
	   int operator ==(int) const;
	   int operator !=(int) const;
	   int operator <(const integer&) const;
	   int operator <=(const integer&) const;
	   int operator >(const integer&) const;
	   int operator >=(const integer&) const;
	   int operator <(int) const;
	   int operator <=(int) const;
	   int operator >(int) const;
	   int operator >=(int) const;
	   int abslower(const integer&) const;
	   // bit-operators
	   integer operator <<(int) const;
	   integer operator >>(int) const;
	   //additive operators
	   integer operator +(const integer&) const;
	   integer operator -(const integer&) const;
	   integer operator -() const;
	   integer& operator +=(const integer&);
	   integer& operator -=(const integer&);
	   integer& operator++();     // praefix ++ Gibt neuen Wert zurück
	   integer& operator++(int);  // postfix ++ Gibt neuen Wert zurück
	   integer& operator--();     // praefix --
	   integer& operator--(int);  // postfix --
	   const integer& changesign() {laenge*=-1;return *this;}
	   //multiplicative operators
	   friend integer shift(const integer&,int);  // shift(x,n) = x * 2^(n*numberofbits)
	   friend integer trunc(const integer&,int);
	   friend integer rest(const integer&,int); // Regel: x = shift(trunc(x,n),n) + rest(x,n)
	   integer operator *(const integer&) const;
	   integer& operator *=(const integer&);
	   friend void divrem(const integer&,const integer&,integer&,integer&);
	   integer operator /(const integer&) const;
	   integer& operator /=(const integer&);
	   integer operator %(const integer&) const;
	   integer& operator %=(const integer&);
	   integer& EDIVBY(const integer&);           
	   friend integer ggTspez(const integer &a,const integer &b);
	   friend integer EDIV(const integer&,const integer&);
	   friend integer gcdeuk(integer,integer);
	   friend DLL_PUBLIC integer abs(const integer& a);
	   friend DLL_PUBLIC integer hilfinteger::exactdivision(integer,integer);
	   friend DLL_PUBLIC integer hilfinteger::restsingle(const integer&,const integer&);
	   // random generator
	   friend integer integerrand(int);
	   // ----------------------------------------
       // friend class 
       friend class rational;
       // input, output:
       enum Output_Type {INTEGER,FLOAT};
       static Output_Type GetOutput_Style() {return Output_Style;}
       static void SetOutput_Style(Output_Type t) {Output_Style=t;}
       static void SetOutput_Digits(int n) {if (n<2) n=2;Output_Digits=n;}
       static integer char_to_integer(char* buf,int l=-1);
	   friend DLL_PUBLIC std::istream& operator >>(std::istream&,integer&);
	   friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,integer);
private:
       // -----------------------------------------------------------------------------
       unsigned int *dat=nullptr;
	   int laenge=0;                             // |laenge| = length of dat. laenge < 0 <=> *this < 0. laenge = 0  <=> *this = 0
	   // -----------------------------------------------------------------------------
	   static int Maxlength;
	   static Output_Type Output_Style;
	   static int Output_Digits;
	   //int abslower(const integer&) const;
	   void add(const integer&,integer&) const;
	   integer& addto(const integer&);
	   void sub(const integer&,integer&) const;
	   void subto(const integer&);
	   void truncby(int);
	   void shiftrightby(int);
	   void classicmult(const integer&,integer&) const;
	   static integer gcdbin(const integer&,const integer&);
       static void divrestsingle(const integer&,const integer&,integer&,integer&);
       static void divsingle(const integer&,const integer&,integer&);
       static void restsingle(const integer&,const integer&,integer&);
       static void divrest(const integer&,const integer&,integer&,integer&);
};


} // end namespace val

#endif //INTEGER_H

