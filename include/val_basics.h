#ifndef VAL_BASICS_H_INCLUDED
#define VAL_BASICS_H_INCLUDED

#include <iostream>
#include <type_traits>



#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  //#define DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    //#define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_PUBLIC
    //#define DLL_LOCAL
  #endif
#endif



namespace val
{

class integer;
class rational;
class modq;
class s_modinteger;
class modint;

// Change Path to your location of primlis1.txt.
// Example if primlis1.txt is in /home/user/include, so set primlispath[] = "/home/user/include/primlis1.txt"
extern DLL_PUBLIC char primlistpath[];
// -----------------------------------------------------------------------------------------------------------

template <class T> struct DLL_PUBLIC is_nongeneric_ring : std::false_type {};
template<> struct DLL_PUBLIC is_nongeneric_ring<rational> : std::true_type {};
template<> struct DLL_PUBLIC is_nongeneric_ring<integer> : std::true_type {};
template<> struct DLL_PUBLIC is_nongeneric_ring<modq> : std::true_type {};
template<> struct DLL_PUBLIC is_nongeneric_ring<s_modinteger> : std::true_type {};
template<> struct DLL_PUBLIC is_nongeneric_ring<modint> : std::true_type {};


template <class T> struct DLL_PUBLIC is_archimedian_ring : std::false_type {};
template <> struct DLL_PUBLIC is_archimedian_ring<integer> : std::true_type {};
template <> struct DLL_PUBLIC is_archimedian_ring<rational> : std::true_type {};



template <class T,class S = T>
struct GPair
{
    T x;
    S y;
    GPair() = default;
    GPair(const T &a,const S &b) : x(a),y(b) {}
    GPair(T&& a,S&& b) : x(std::move(a)), y(std::move(b)) {}
};


DLL_PUBLIC int gcd(int a,int b);

template <class T>
void swap(T& a,T& b)
{
    T h(std::move(b));
    b=std::move(a);
    a=std::move(h);
}


template <class T>
inline const T& Max(const T& a,const T& b)
{
    if (a>b) return a;
    else return b;
}

template <class T>
inline const T& Min(const T& a, const T& b)
{
    if (b<a) return b;
    else return a;
}




template <typename T,typename std::enable_if_t<std::is_pointer<T>::value || std::is_arithmetic<T>::value,int> = 0 >
T zero_element(void)
{
    return T(0);
}


template <typename T, typename std::enable_if_t<!(std::is_pointer<T>::value || std::is_arithmetic<T>::value),int> = 0 >
T zero_element(void)
{
    T zero;
    return zero;
}



template <typename T,typename std::enable_if_t<std::is_arithmetic<T>::value || is_nongeneric_ring<T>::value,int> = 0 >
void set_unity_element(T &one)
{
    one = T(1);
}


template <typename T,typename std::enable_if_t<std::is_arithmetic<T>::value || is_nongeneric_ring<T>::value,int> = 0 >
T unity_element()
{
    return T(1);
}


template <typename T,typename std::enable_if_t<!std::is_arithmetic<T>::value && !is_nongeneric_ring<T>::value,int> = 0 >
T unity_element()
{
    T one;
    set_unity_element(one);
    return one;
}


template <typename T,typename std::enable_if_t<std::is_arithmetic<T>::value,int> = 0 >
T abs(const T& a)
{
    if (a<T(0)) return -a;
    else return a;
}

template <class T>
const T& as_const(const T& W)
{
    return const_cast<const T&> (W);
}


template <class T,class U>
const T& MemCopy(T& dest,const U& src)
{
    int n=Min(sizeof(dest),sizeof(src));
    const char* s = (char*) &src;
    char* d = (char*) &dest;
    for (int i=0;i<n;++i) d[i] = s[i];
    return dest;
}

extern DLL_PUBLIC const double Inf,NaN;

DLL_PUBLIC double sqrt(const double& a);

DLL_PUBLIC double round(const double& x,int k=4); // round to e-k if k>=0.


template <typename T, typename std::enable_if_t<!std::is_floating_point<T>::value,int> = 0 >
T power(const T& a,int z)
{
    int n=abs(z);
    T y,x=unity_element<T>(),zero=zero_element<T>();

    if (a==zero) {
        if (z==0) return x;   // = 1
        else if (z>0) return zero;
        else return x/a;
    }

    if (z==0) return x; // = 1
    //if (z>0) y=a;
    //else y=1.0/a;
    y=a;

    while (n!=0) {
       if (n%2!=0) {
           x*=y;
       }
       y*=y;
       n=n/2;
    }
    if (z<0) return (unity_element<T>()/x);
    else return x;
}


} // end namespace val

#endif // VAL_BASICS_H_INCLUDED
