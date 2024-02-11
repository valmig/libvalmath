#ifndef RATION_CPP
#define RATION_CPP

#include <rational.h>
#include "integer.cpp"

namespace val
{


integer::integer(const rational& r) :integer(nominator(r)/denominator(r))
{

}



// ------------------------------ Hilfsfktn -------------------------------------------



void hilfratiofktn::errornenner()
{
 std::cout<<"\nDIVISION BY 0 (in RATIO.H)!!!!!                  ";
 std::cout<<"\nABORT PROGRAM!!!!!                  ";
 exit(-1);
}

// ---------------------------------------------------------------------------------------




char rational::gekuerzt = 1;
rational::Output_Type rational::Output_Style=RATIONAL;
int rational::Output_Digits=16;

void rational:: kuerze()
{
    integer h(ggTspez(zaehler,nenner));
    if (h.dat == NULL) return;
    if (h.laenge < 0) h.laenge *= -1;
    zaehler = EDIV(zaehler, h);
    nenner=EDIV(nenner,h);
    if (nenner.laenge<0) {zaehler.laenge*=-1;nenner.laenge*=-1;}
}


// ------------------------------ Constructors ------------------------------------------
rational::rational(int a,int b) : zaehler(a),nenner(b)
{
   if (gekuerzt) kuerze();
}


rational:: rational(const integer& a,const integer& b) : zaehler(a) , nenner(b)
{
   if (gekuerzt) kuerze();
}


rational::rational(const double &a)
{
    int sign=1,exp=0,i;
    unsigned bit=1;
    __int64 b,eins(1);
    b =  *((__int64*) &a);

    if ((eins<<63) & b) sign=-1;

    // Exponent:
    for (i=62;i>51;i--) {
        if ((eins<<i) &b) exp+=1;
        exp<<=1;
    }
    exp>>=1;
    
    // Special cases:
    if (exp==0) {nenner=integer(1);return;} // sub-normals will be treated as 0

    if (exp==2047) {   //+ or - inf or nan
        int nan=0;
        for (i=51;i>=0;i--) {
            if ((eins<<i) & b) {nan=1;break;}
        }
        if (nan) return;
        zaehler=integer(1);
        zaehler.laenge*=sign;
        return;
    }
    exp-=1023;
    exp-=52;
    
    // Get mantissa;
    zaehler=integer(1)<<52;   // = 2^52.
    for (i=51;i>=32;i--) {
         if ((eins<<i) & b) zaehler.dat[1] |= (bit << (i%32));
    }
    for (;i>=0;i--) {
         if ((eins<<i) & b) zaehler.dat[0] |= (bit << (i%32));
    }
    //
    if (exp<0) nenner=integer(1)<<hilfinteger::abs(exp);
    else {zaehler = zaehler << exp;nenner=integer(1);}

    if (sign==-1) zaehler.laenge*=-1;
    kuerze();
    return;

}

// -----------------------------------------------------------------------------------------

// -------------------Cast operators ---------------------------------------------------
rational::operator double() const
{
 int lz=hilfinteger::abs(zaehler.laenge),ln=hilfinteger::abs(nenner.laenge);

 if (nenner.dat==NULL && zaehler.dat==NULL) return val::NaN;

 if (zaehler.dat==NULL) return 0.0;

 if (lz-ln>32) {
    if (zaehler.laenge>0) return val::Inf;
    else return -val::Inf;
 }
 if (ln-lz>33) {
    if (zaehler.laenge>0) return double(0.0/1.0);
    else return double (0.0/-1.0);
 }

 if (lz <=32 && ln <=32) {
        return (double(zaehler) / double(nenner));
 }
 else if (lz>ln) { // not ready yet
        integer z1=trunc(zaehler,lz-32),n1=trunc(nenner,lz-32);
        return (double(z1) / double(n1));
 }
 else { //not ready yet
        integer z1=trunc(zaehler,ln-32),n1=trunc(nenner,ln-32);
        return (double(z1) / double(n1));
 }
}


// ------------------------------ Input/output ------------------------------------------



rational rational::char_to_rational(char* buf,int l)
{
  int i,sign=1,exp=0,kommastellen=0,signexp=1;
  rational ten(10),x;

  if (buf==NULL) return x;

  if (l==-1) for (l=0;buf[l]!='\0';l++) ;

  if (l==0) return x;


  //delete all 0 and - at the beginning
  for (i=0;i<l;) {
	  if (buf[i]=='-') {
		  sign*=-1;
		  i++;
	  }
	  else if (buf[i]=='0') {
		  i++;
	  }
	  else break;
  }
  if (i==l) {
     return x;
  }
  // Get digits:
  for (;i<l;i++) {
      if (buf[i]=='.' || buf[i]==',') {
        if (kommastellen) continue;
        kommastellen=1;
        continue;
      }
	  if (buf[i]=='e') {i++;break;}
	  if (kommastellen) kommastellen++;
	  x*=ten;
	  x+=integer(int(buf[i]-48));
  }

  if (kommastellen) kommastellen--;

  for (;i<l;i++) {
      if (buf[i]=='+' || buf[i]=='0') continue;
      else if (buf[i]=='-') {signexp*=-1;continue;}
      else break;
  }
  for (;i<l;i++) {
      exp*=10;
      exp+= int(buf[i]-48);
  }

  if (signexp==-1) exp*=-1;

  exp-=kommastellen;

  if (exp>=0) {
    for (i=0;i<exp;i++) x*=ten;
  }
  else {
    exp*=-1;
    for (i=0;i<exp;i++) x/=ten;
  }

  if (sign==-1) x.zaehler.laenge*=-1;
  return x;
}



std::istream& operator >>(std::istream& is,rational& x)
{

  char *buf,*hbuf;
  int i,l,k=1,wert=0;//sign=1;
  integer ten(10);

  buf = new char[1000];



  l=0;
  while ((wert!=10) && (wert!=32) && is) {
	wert=is.get();
	if ((wert==10) || (wert==32)){
		if (l!=0) buf[l]='\0';
		else wert=0;
	}
	else buf[l++]=wert;
	if (l>=(k*1000)) {
        k++;
        hbuf= new char[k*1000];
        for (i=0;i<((k-1)*1000);i++) hbuf[i]=buf[i];
        delete[] buf;
        buf=hbuf;
        hbuf=NULL;
	}

	if (!is) l--;
  }

  // Search '/'
  for (i=0;i<l && buf[i]!='/';i++);
  x=rational::char_to_rational(buf,i);
  if (i<l-1) {
    hbuf=&buf[i+1];
    x/=rational::char_to_rational(hbuf,l-1-i);
    hbuf=NULL;
  }

  delete[] buf;
  return is;
}


std::ostream& operator <<(std::ostream& os,const rational& y)
{
    if (rational::Output_Style==rational::RATIONAL) {
        integer::Output_Type style=integer::GetOutput_Style();
        integer::SetOutput_Style(integer::INTEGER);
        os<<y.zaehler;
        if (y.nenner!=1) os<<"/"<<y.nenner;
        integer::SetOutput_Style(style);
    }
    else {

        if (y.nenner.isNull()) {os<<double(y);return os;}
        if (y.zaehler.isNull()) {os<<'0'; return os;}

        integer x,q,r,faktor(10),ten(10),eins(1),z;
        int i=1,exp=0,osleer=1,komma=0,numberofzeros=0,sign=1;

        if (y.zaehler.signum()<0) sign=-1;

        divrem(y.zaehler,y.nenner,q,r);
        x=std::move(q);
        z=std::move(r);
        //digits left of comma:
        if (!x.isNull()) {
            if (sign==-1) os<<'-';
            osleer=0;
            while (faktor.abslower(x)) {faktor*=ten;exp++;}

            while (!x.isNull() && i<=rational::Output_Digits) {
                faktor/=ten;exp--;
                divrem(x,faktor,q,r);
                if (q.isNull()) numberofzeros++;
                else if (unsigned(q)==10) {
                    os<<'1';i++;numberofzeros++;
                }
                else {
                     if (!komma && i>1) {komma=1;os<<'.';}
                     for(;numberofzeros>0;numberofzeros--) os<<'0';
                     os<<unsigned(q);
                }
                i++;
                x=std::move(r);
                q=0;
            }
            while (eins.abslower(faktor) && i<=rational::Output_Digits) {
                faktor/=ten;
                i++;
                exp--;
                numberofzeros++;
            }
        }
        exp+=i-1;
        if (i>rational::Output_Digits || z.isNull()) {
                if (exp>0) os<<"e+0"<<exp;
                return os;
        }
        //digits right from comma
        z*=ten;
        if (osleer) exp--;
        while (z.abslower(y.nenner) && i<=rational::Output_Digits) {
            z*=ten;
            if (osleer) exp--;
            else {i++;numberofzeros++;}
        }
        if (osleer && i>rational::Output_Digits) {os<<'0';return os;}

        while (!z.isNull() && i<=rational::Output_Digits) {
            divrem(z,y.nenner,q,r);
            if (q.isNull() && !osleer) {numberofzeros++;i++;}
            else if (!q.isNull()) {
                    if (!komma && !osleer) {os<<'.';komma=1;}
                    if (osleer && sign==-1) os<<'-';
                    for (;numberofzeros>0;numberofzeros--) os<<'0';
                    os<<unsigned(q);i++;osleer=0;
            }
            z=std::move(r);
            z*=ten;
            r=0;
        }

        if (exp!=0) {
            if (exp>0) os<<"e+0"<<exp;
            else os<<"e-0"<<hilfinteger::abs(exp);
        }

        if (osleer) os<<'0';

	}
    return os;
}
// -----------------------------------------------------------------------------------------


// --------------------------- Additive Operatoren -----------------------------------------

rational rational:: operator +(const rational& x) const
{
  // Cases inf and NaN
  if (nenner.dat==NULL && x.nenner.dat==NULL) {
     if (zaehler.laenge==1 && x.zaehler.laenge==1) return rational(1,0);
     else if (zaehler.laenge==-1 &&x.zaehler.laenge==-1) return rational(-1,0);
     else return rational(0,0);
  }

  rational y(zaehler*(x.nenner)+(x.zaehler)*nenner,nenner*(x.nenner));
  return y;
}

rational rational:: operator -(const rational& x) const
{
  // Cases inf and NaN
  if (nenner.dat==NULL && x.nenner.dat==NULL) {
     if (zaehler.laenge==1 && x.zaehler.laenge==-1) return rational(1,0);
     else if (zaehler.laenge==-1 &&x.zaehler.laenge==1) return rational(-1,0);
     else return rational(0,0);
  }
  rational y(zaehler*x.nenner-x.zaehler*nenner,nenner*x.nenner);
  return y;
}

rational rational::operator -() const
{
 rational y;

 y.nenner=nenner;
 y.zaehler=zaehler;
 (y.zaehler).laenge*=-1;
 return y;
}

// -----------------------------------------------------------------------------------------


// ------------------------ Multiplicative operators --------------------------------------

rational rational:: operator *(const rational& x) const
{
  rational y(zaehler*x.zaehler,nenner*x.nenner);
  return y;
}


rational rational:: operator /(const rational& x) const
{
 rational y(zaehler*x.nenner,nenner*x.zaehler);
 return y;
}


rational rational:: operator /(const integer& n) const
{
   rational y(zaehler,nenner*n);
   return y;
}


rational operator *(const integer& z,const rational& x)
{
 rational y(z*x.zaehler,x.nenner);
 return y;
}

// ---------------------------------------------------------------------------------------


// --------------------------- Zuweisungsoperatoren --------------------------------------

rational& rational:: operator =(const rational& x)
{
 zaehler=x.zaehler;
 nenner=x.nenner;
 return *this;
}


rational& rational::operator =(rational&& r)
{
    zaehler=std::move(r.zaehler);
    nenner=std::move(r.nenner);
    return *this;
}


rational& rational:: operator=(const integer& n)
{
  zaehler=n;nenner=integer(1);
  return *this;
}

rational& rational::operator =(integer&& r)
{
    zaehler=std::move(r);
    nenner=integer(1);
    return *this;
}


rational& rational:: operator +=(const rational& x)
{
  *this = *this + x;
   return *this;
}


rational& rational:: operator -=(const rational& x)
{
 *this = *this - x;
 return *this;
}

rational& rational:: operator *=(const rational& x)
{
 *this = *this * x;
 return *this;
}


rational& rational:: operator /=(const rational& x)
{
 *this = *this / x;
 return *this;
}

// --------------------------------------------------------------------------------------


// -------------------------- Comparisons --------------------------------------

int rational :: operator ==(const rational& x) const
{
	// Special cases
	if (nenner.dat==NULL || x.nenner.dat==NULL) return 0;
	//
	if (gekuerzt) {
		return (zaehler==x.zaehler && nenner==x.nenner);
	}
	else return ((zaehler*x.nenner)==(nenner*x.zaehler));
}

int rational :: operator !=(const rational& x) const
{
   return !(*this==x);
}

int rational :: operator ==(const integer& n) const
{
  if (gekuerzt) return ((nenner==1) && (zaehler==n));
  else return (zaehler==(nenner*n));
}

int rational :: operator !=(const integer& n) const
{
 	return !(*this==n);
}

int rational::operator >(const rational& x) const
{
	return ((zaehler*x.nenner)>(nenner*x.zaehler));
}

int rational::operator <(const rational& x) const
{
	return ((zaehler*x.nenner)<(nenner*x.zaehler));
}

int rational::operator >=(const rational& x) const
{
	return ((zaehler*x.nenner)>=(nenner*x.zaehler));
}

int rational::operator <=(const rational& x) const
{
	return ((zaehler*x.nenner)<=(nenner*x.zaehler));
}


// ---------------------------------------------------------------------------------------

rational abs(const rational& x) 
{
 rational y(x);
 if (y.zaehler.signum()<0) y.zaehler.changesign();
 return y;
}

} // End namespace val

#endif // RATION_CPP
