#ifndef COMPLEX_TYPE_H_INCLUDED
#define COMPLEX_TYPE_H_INCLUDED

#include <val_basics.h>
#include <string>
#include <sstream>


namespace val
{

class rational;
template <class T> class complex_type;
template <class T> std::istream& operator >>(std::istream&,complex_type<T>&);
template <class T> complex_type<T> operator *(const T&,const complex_type<T>&);
template <class T> double abs(const complex_type<T>&);
template <class T>
void set_unity_element(complex_type<T> &one);
//template <class T> std::ostream& operator <<(std::ostream&,const complex_type<T>&);

template <class T>
class complex_type
{
private:
    static const T T_zero;
    T re = T_zero, im = T_zero;
    static T To_number(char*,int);
public:
    complex_type() = default;
    complex_type(const T& x) : re(x) , im(T_zero) {}
    complex_type(T&& x) : re(std::move(x)) , im(T_zero) {}
    complex_type(const T &x,const T &y) : re(x) , im(y) {}
    complex_type(T&& x, T&& y) : re(std::move(x)), im(std::move(y)) {}
    complex_type(const complex_type<T>&) = default;
    complex_type(complex_type &&z) : re(std::move(z.re)) , im(std::move(z.im)) {}
    //
    const T& real() const {return re;}
    const T& imaginary() const {return im;}
    int is_zero() const {return (re==T_zero && im==T_zero);}
    //
    const complex_type<T>& operator =(const complex_type<T>& z) {re=z.re;im=z.im;return *this;}
    const complex_type<T>& operator =(complex_type<T>&& z) {re=std::move(z.re); im = std::move(z.im); return *this;}
    const complex_type<T>& operator =(const T& x) {re=x;im=T_zero;return *this;}
    const complex_type<T>& operator =(T&& x) {re=std::move(x);im=T_zero;return *this;}
    //
    int operator ==(const complex_type<T> &z) const {return (re==z.re && im == z.im);}
    int operator !=(const complex_type<T> &z) const {return !(*this==z);}
    //
    complex_type<T> conjugate() const {return complex_type<T>(re,-im);}
    //
    complex_type<T> operator +(const complex_type<T>& z) const {return complex_type<T>(re+z.re,im+z.im);}
    const complex_type<T>& operator +=(const complex_type<T> &z) {*this = *this + z; return*this;}
    complex_type<T> operator -(const complex_type<T>& z) const {return complex_type<T>(re-z.re,im-z.im);}
    const complex_type<T>& operator -=(const complex_type<T> &z) {*this = *this - z; return*this;}
    complex_type<T> operator -() const {return complex_type<T>(-re,-im);}
    //
    complex_type<T> operator *(const complex_type<T>& x) const {return complex_type<T>(re*x.re - im*x.im,re*x.im + im*x.re);}
    const complex_type<T>& operator *=(const complex_type<T> &z) {*this = *this * z; return*this;}
    complex_type<T> operator /(const complex_type<T>& x) const {return complex_type<T>((re*x.re + im*x.im)/(x.re*x.re + x.im*x.im),(im*x.re - re*x.im)/(x.re*x.re + x.im*x.im));}
    const complex_type<T>& operator /=(const complex_type<T> &z) {*this = *this / z; return*this;}
    //
    friend std::istream& operator >><T>(std::istream&,complex_type<T>&);
    //
    static double outputerror;
};


//
template <class T>
const T complex_type<T>:: T_zero = zero_element<T>();

template <class T>
double complex_type<T>::outputerror=1e-09;

//

template <class T>
complex_type<T> operator *(const T& x,const complex_type<T>& z)
{
    return (x*z.real(),x*z.imaginary());
}

template <class T>
double abs(const complex_type<T> &z)
{
    return sqrt(double(z.real()*z.real() + z.imaginary()*z.imaginary()));
}
//
template <typename T, std::enable_if_t<!std::is_floating_point<T>::value,int> = 0>//,typename std::enable_if_t<!std::is_floating_point<T>>,int> = 0>
std::ostream& operator <<(std::ostream& os,const complex_type<T>& z)
{
    T h=z.imaginary(),zero=zero_element<T>();

    if (z.is_zero()) {
        os<<'0';
        return os;
    }
    if (z.real()!=zero) os<<z.real();
    if (h!=zero) {
        if (h>zero && z.real()!=zero) os<<"+";
        if (z.imaginary()!=T(1)) os<<h;
        os<<"i";
    }
    return os;
}


template <typename T, std::enable_if_t<std::is_floating_point<T>::value,int> = 0>//,typename std::enable_if_t<!std::is_floating_point<T>>,int> = 0>
std::ostream& operator <<(std::ostream& os,const complex_type<T>& x)
{
    T zero(0);

    if (x.real()==zero && abs(x.imaginary())<=complex_type<T>::outputerror) {
        os<<'0';
        return os;
    }
    if (x.real()!=zero) os<<x.real();
    if (abs(x.imaginary())>complex_type<T>::outputerror) {
        if (x.imaginary()>0 && x.real()!=zero) os<<'+';
        os<<x.imaginary()<<'i';
    }
    return os;
}


template <class T>
T complex_type<T>::To_number(char *buf,int l)
{
    T x=T_zero;
    if (l==0) return x;

    std::stringstream ss;
    std::string s;


    for (int i=0;i<l;++i) s+=buf[i];
    //std::cout<<std::endl<<s;
    ss.str(s);
    ss>>x;
    return x;
}


template <class T>
std::istream& operator >>(std::istream& is,complex_type<T>& x)
{
  x=complex_type<T>();
  T zero=complex_type<T>::T_zero;
  char *buf,*hbuf;
  int i,l,k=1,wert=0,sign=1;
  //integer ten(10);

  buf = new char[1000];



  l=0;
  //is.seekg(1l);
  while ((wert!=10) && (wert!=32) && is) {
	wert=is.get();
	//cout<<wert;
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

  //suche i: nach i werden alle Zeichen ignoriert:
  for (i=0;i<l;i++) {
       if (buf[i]=='i') break;
  }
  if (i==l) {// Imaginär-Teil = 0
      x.re=complex_type<T>::To_number(buf,i);
      delete[] buf;
      //std::cout<<"\n is (re) ready: "<<x<<std::endl;;
      return is;
  }
  // Nun ist Imginär-Teil !=0:
  l=i;
  // suche nach Trennendem + oder -
  i=0;
  while (i<l && (buf[i]=='+' || buf[i]=='-')) {
        if (buf[i]=='-') sign*=-1;
        i++; //+ und Minus werden zunächst ignoriert.
  }

  if (i==l) {
      x.im=T(sign)*T(1);
      delete[] buf;
      //std::cout<<"\n is (im) ready";
      return is;
  }
  //Nun ist i<l und buf[i]!=+ -
  for (;i<l;i++) {
       if (buf[i]=='+' || buf[i]=='-') {
           if (buf[i-1]!='e') break;
       }
  }
  if (i==l) {  // Zahl ist rein -imaginär
      x.im=complex_type<T>::To_number(buf,i);
      delete[] buf;
      return is;
  }
  //x.zaehler=integer::char_to_integer(buf,i);
  //Zahl hat im und re -Teil != 0
  x.re=complex_type<T>::To_number(buf,i);
  if (i<l-1) {
    hbuf=&buf[i+1];
    //x.nenner=integer::char_to_integer(hbuf,l-1-i);
    x.im=complex_type<T>::To_number(hbuf,l-1-i);
    if (buf[i]=='-') x.im*=-1;
    if (x.im==zero) {
        sign=1;
        x.im=T(1);
        for (;i<l;i++) {
            if (buf[i]=='-') sign*=-1;
            else {
                x.im=zero;
                break;
            }
        }
        x.im*=sign;
    }
    hbuf=NULL;
  }
  else {
    if (buf[i]=='-') x.im =T(-1);
    else x.im=T(1);
  }

  //if (x.nenner==0) hilfratiofktn::errornenner();

  //if (rational::gekuerzt) x.kuerze();

  delete[] buf;
  return is;
}

//
template <class T>
void set_unity_element(complex_type<T> &one)
{
    T one_T;
    set_unity_element(one_T);
    one = complex_type<T>(one_T);
}

//
typedef complex_type<rational> gauss_number;
typedef complex_type<double> complex;

} // end namespace val




#endif // COMPLEX_TYPE_H_INCLUDED
