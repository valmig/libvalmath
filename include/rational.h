#ifndef RATION_H
#define RATION_H

// class for rational numbers

#include <integer.h>


namespace val
{


namespace hilfratiofktn
{
 DLL_PUBLIC void errornenner();
}

//===================================================
// RULES:      1/0 = inf , -1/0 = -inf  , 0/0 = nan
// ==================================================


DLL_PUBLIC const integer& nominator(const rational& a);
DLL_PUBLIC const integer& numerator(const rational& a);
DLL_PUBLIC const integer& denominator(const rational& a);
DLL_PUBLIC rational abs(const rational&);


class DLL_PUBLIC rational
{
public:
		rational() : zaehler(0), nenner(1) {}
		rational(int a) : zaehler(a),nenner(1) {}
	    rational(int, int);
		rational(const integer& a) : zaehler(a), nenner(1) {}
		rational(integer&& a) : zaehler(std::move(a)),nenner(1) {}
		rational(const integer&,const integer&);
		rational(integer&& a,integer&& b) :zaehler(std::move(a)),nenner(std::move(b)) {kuerze();}
		rational(const rational& r) : zaehler(r.zaehler),nenner(r.nenner) {}
		rational(rational&& r) : zaehler(std::move(r.zaehler)),nenner(std::move(r.nenner)) {}
	    rational (const double&);
	    void kuerze();                                               // reduce
	    void reduce() {kuerze();} 
	    // Input/Output:
	    enum Output_Type {RATIONAL,FLOAT};
        static Output_Type GetOutput_Style() {return Output_Style;}
        static void SetOutput_Style(Output_Type t) {Output_Style=t;}
        static void SetOutput_Digits(int n) {if (n<2) n=2;Output_Digits=n;}
	    friend DLL_PUBLIC std::istream& operator >>(std::istream&,rational&);
	    friend DLL_PUBLIC std::ostream& operator <<(std::ostream&,const rational&);
	    // cast-operators
		operator double () const;
        // const-getters:
        const integer& nominator() const {return zaehler;}
        const integer& numerator() const {return zaehler;}
        const integer& denominator() const {return nenner;}
		friend const integer& nominator(const rational& a) {return a.zaehler;}
		friend const integer& numerator(const rational& a) {return a.zaehler;}
		friend const integer& denominator(const rational& a) {return a.nenner;}
        //
	    rational operator +(const rational&) const;
	    rational operator -(const rational&) const;
		rational operator -() const;
	    rational operator *(const rational&) const;
	    rational operator /(const rational&) const;
	    rational operator /(const integer&) const;
	    friend DLL_PUBLIC rational operator *(const integer&,const rational&);
	    rational& operator =(const rational&);
	    rational& operator =(rational&&);
	    rational& operator =(const integer&);
	    rational& operator =(integer&&);
	    rational& operator +=(const rational&);
	    rational& operator -=(const rational&);
	    rational& operator *=(const rational&);
	    rational& operator /=(const rational&);
	    int operator ==(const rational&) const;
	    int operator !=(const rational&) const;
	    int operator ==(const integer&) const;
	    int operator !=(const integer&) const;
		int operator >(const rational&) const;
		int operator <(const rational&) const;
		int operator >=(const rational&) const;
		int operator <=(const rational&) const;
		friend rational abs(const rational&);        
		int signum() const {return zaehler.signum();}
		int abslength() const {return zaehler.abslength()+nenner.abslength();}
        const rational& changesign() {zaehler.changesign();return *this;}
        int isNull() const {return zaehler.isNull();}
        int iszero() const {return (zaehler.dat==NULL && nenner==integer(1));}
private:
        //
	    integer zaehler,nenner;
        //
        static char gekuerzt;
	    static Output_Type Output_Style;
	    static int Output_Digits;
        static rational char_to_rational(char* buf,int l=-1);
};

} // End namespace val


#endif
