#ifndef FRACTION_H_INCLUDED
#define FRACTION_H_INCLUDED

#include <error.h>
#include <iostream>
#include <val_basics.h>


namespace val
{

template <class T> class pol;
template <class T> class fraction;
template<typename T> struct is_pol_class : std::false_type {};
template<typename T> struct is_pol_class<val::pol<T>> : std::true_type {};
template <class T> void set_unity_element(fraction<T> &one);
class rational;
typedef fraction<pol<rational>> rationalfunction;

template <class T> class fraction;
template <class T> std::ostream& operator <<(std::ostream&,const fraction<T>&);

template <class T>
class fraction
{
private:
    T zaehler=zero_element<T>();      // zaehler = numerator
    T nenner=unity_element<T>();  	  // nenner = denominator
    static int reduced;

    template<class U=T,typename std::enable_if_t<!is_pol_class<U>::value,int> = 0>
    void reduce() {}
    template<class U=T,typename std::enable_if_t<is_pol_class<U>::value,int> = 0>
    void reduce() {if (!reduced) return;T div(gcd(zaehler,nenner));zaehler/=div;nenner/=div;}
    //void reduce(typename std::true_type);
public:
    fraction() = default;
    //fraction() : zaehler(zero_element<T>()) {set_unity_element(nenner);}
    fraction(const T &z) : zaehler(z) { }
    fraction(T &&z) : zaehler(std::move(z)) { }
    fraction(const T &z,const T &n) : zaehler(z) , nenner(n) {reduce();}//reduce(typename is_pol_class<T>::type());}
    fraction(T&& z,T&& n) : zaehler(std::move(z)) , nenner(std::move(n)) {reduce();}//reduce(typename is_pol_class<T>::type());}
    fraction(const fraction<T> &F) : zaehler(F.zaehler) , nenner(F.nenner) {}
    fraction(fraction<T>&& F) : zaehler(std::move(F.zaehler)) , nenner(std::move(F.nenner)) {}
    const fraction<T>& operator =(const fraction<T>& F) {zaehler=F.zaehler;nenner=F.nenner;return *this;}
    const fraction<T>& operator =(fraction<T>&& F) {zaehler=std::move(F.zaehler);nenner=std::move(F.nenner); return *this;}
    const fraction<T>& operator =(const T& f) {zaehler=f;set_unity_element(nenner);; return *this;}
    const fraction<T>& operator =(T&& f) {zaehler=std::move(f);set_unity_element(nenner);; return *this;}
    const T& nominator() const {return zaehler;}
    const T& numerator() const {return zaehler;}
    const T& denominator() const {return nenner;}
    fraction<T> operator +(const fraction<T> &x) const {fraction<T> y(zaehler*(x.nenner)+(x.zaehler)*nenner,nenner*(x.nenner));return y;}
    fraction<T> operator -(const fraction<T> &x) const {fraction<T> y(zaehler*(x.nenner)-(x.zaehler)*nenner,nenner*(x.nenner)); return y;}
    fraction<T> operator -() const {fraction<T> r;r.zaehler=-zaehler;r.nenner=nenner;return r;}
    const fraction<T>& operator +=(const fraction<T> &F) {*this=*this+F;return *this;}
    const fraction<T>& operator -=(const fraction<T> &F) {*this=*this-F;return *this;}
    fraction<T> operator *(const fraction<T> &F) const {return fraction<T>(zaehler*F.zaehler,nenner*F.nenner);}
    fraction<T> operator /(const fraction<T> &F) const;
    const fraction<T>& operator *=(const fraction<T> &F) {*this=*this*F;return *this;}
    const fraction<T>& operator /=(const fraction<T> &F) {*this=*this/F;return *this;}
    int operator ==(const fraction<T> &F) const {return  (zaehler==F.zaehler && nenner==F.nenner);}
    int operator !=(const fraction<T> &F) const {return !(*this==F);}
    static int isreduced() {return reduced;}
    static void setreduced(int a) {reduced=a;}
    void reducethis() {reduce();}
    friend std::ostream& operator <<<T>(std::ostream&,const fraction<T>&);
};

/*
template <class T>
//void fraction<T>::reduce(std::true_type)
typename std::enable_if<is_pol_class<T>::value,void> reduce()
{
    T div(gcd(zaehler,nenner));
    zaehler/=div;
    nenner/=div;
}
*/

template <class T>
int fraction<T>::reduced=1;

template <class T>
fraction<T> fraction<T>::operator /(const fraction<T> &F) const
{
    T zero = zero_element<T>();
    if (F.zaehler==zero) Error::error("Error: fraction<T>::operator /: division by zero!;");
    return fraction<T>(zaehler*F.nenner,nenner*F.zaehler);
}



template <class T>
std::ostream& operator <<(std::ostream &os,const fraction<T> &F)
{
    T one;
    set_unity_element(one);
    os<<F.zaehler;
    if (F.nenner!=one) std::cout<<"\n/\n"<<F.nenner;
    return os;
}


template <class T>
void set_unity_element(fraction<T> &one)
{
    T one_T;
    set_unity_element(one_T);
    one = fraction<T>(one_T);
}



}



#endif // FRACTION_H_INCLUDED
