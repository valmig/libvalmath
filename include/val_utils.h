#ifndef VAL_UTILS_H_INCLUDED
#define VAL_UTILS_H_INCLUDED

//#include <val_utils.cpp>
#include <val_basics.h>
#include <string>
#include <sstream>
#include <typeinfo>
#ifdef __GNUC__
#include <cxxabi.h>
#endif

namespace val
{

class integer;
class rational;
template <class T> class complex_type;
typedef complex_type<double> complex;

template <class T> struct DLL_PUBLIC is_string_type : std::false_type {};
template <typename T> struct DLL_PUBLIC is_string_type<std::basic_string<T>> : std::true_type {};
//template <> struct DLL_PUBLIC is_string_type<std::wstring> : std::true_type {};



template <class T> std::string ToString(const T& a,unsigned prec=6);

template <class T> std::string gettypename(const T& a);

/*
template <typename S, typename std::enable_if_t<val::is_string_type<S>::value,int> = 0>
S headofstring(const S &value, int m  = 1);

template <typename S, typename std::enable_if_t<val::is_string_type<S>::value,int> = 0>
S tailofstring(const S &value, int m  = 1);
*/

//template <template<typename> class C, typename S, typename std::enable_if_t<val::is_string_type<S>::value,int> = 0>
//S getfirstwordofstring(const S& value, const C<char> &separators = C<char>{'\n', ' ', ';', ',', '.', ':'});
template <template<typename> class C, typename T  = char>
std::basic_string<T> getfirstwordofstring(const std::basic_string<T>& value, const C<T> &separators = C<T>{'\n', ' ', ';', ',', '.', ':'});


DLL_PUBLIC char* StringToChar(const std::string& s);


template <typename T,template <typename> class C>
int isinContainer(const T& value, const C<T>& G);

//DLL_PUBLIC rational char_to_rational(const char*);
//DLL_PUBLIC rational string_to_rational(const std::string&);
//DLL_PUBLIC val::complex string_to_complex(const std::string&);

template <class T>
std::string ToString(const T& a,unsigned prec)
{
   std::stringstream ss;//create a stringstream
   ss.precision(prec);
   ss <<a;//add number to the stream
   return ss.str();
}

//template <class T>
template <typename T, typename std::enable_if_t<!std::is_floating_point<T>::value,int> = 0 >
T FromString(const std::string& s)
{
   std::stringstream ss;
   T x;

   ss.str(s);
   ss>>x;
   return x;
}


template <typename T, typename std::enable_if_t<std::is_floating_point<T>::value,int> = 0 >
T FromString(const std::string& s)
{
   std::stringstream ss;
   T x;

   if (s=="inf") return T(val::Inf);
   if (s=="-inf") return T(-val::Inf);

   ss.str(s);
   ss>>x;
   return x;
}


template <typename S, typename std::enable_if_t<val::is_string_type<S>::value,int> = 0> 
S headofstring(const S &value, int m)
{
	int n = value.length();
	S head;
	
	m = Min(n,m);
	
	if (m == 0) return head;
	
	for (int i = 0; i < m; ++i) {
		head += value[i];
	}
	return head;
}


template <typename S, typename std::enable_if_t<val::is_string_type<S>::value,int> = 0> 
S tailofstring(const S &value,int m)
{
	int n = value.length();
	S tail;
	
	m = Min(m,n);
	
	for (int i = n-m; i < n; ++i) {
		tail += value[i];
	}
	return tail;	
}

//template <template<typename> class C, typename S,typename std::enable_if_t<val::is_string_type<S>::value,int> = 0>
//S getfirstwordofstring(const S& value, const C<char> &separators)
template <template<typename> class C, typename T>
std::basic_string<T> getfirstwordofstring(const std::basic_string<T>& value, const C<T> &separators)
{
	std::basic_string<T> first, empty;
	int n = value.length();
	for (int i = 0; i < n; ++i) {
		if (isinContainer(value[i],separators)) {
			if (first != empty) return first;
		}
		else first += value[i];
	}
	return first;
}



template <class T> std::string gettypename(const T& a)
{
    std::string name(typeid(a).name());

#ifdef __GNUC__
	int status=0;
    char *cname = abi::__cxa_demangle(name.c_str(), nullptr, nullptr,&status);
    if (status==0) {
        name = std::string(cname);
    }
    free(cname);
#endif

    return name;
}


// Check if s is an integer.
//DLL_PUBLIC int isinteger(const std::string &s);
template <class STRG_TYPE>
int isinteger(const STRG_TYPE &s)
{
    int l = s.length();
    if (!l) return 0;

    for (int i = 0; i < l; ++i) {
        if (i==0 && int(s[i]) == int('-') ) continue;
        else if (s[i]<48 || s[i] > 57) return 0;
    }
    return 1;
}

template <class STRG_TYPE>
int isrationalnumber(const STRG_TYPE &s)
{
	int l = s.length();
	if (!l) return 0;
	
	
	for (int i = 0; i < l; ++i) {
		if (s[i] == '-' || s[i]  == ',' || s[i] =='.' || s[i] == 'e' || s[i] == '/' ) continue;
		else if (s[i]<48 || s[i] > 57) return 0;
	}
	return 1;
}



// ---

template <class STRG_TYPE>
STRG_TYPE capitalize(const STRG_TYPE &word)
{
    int l = word.length(), c;
    STRG_TYPE n_word;

    for (int i = 0; i < l; ++i) {
        c = int(word[i]);
        if (c < 0 && i < l-1) {
            n_word += word[i];
            ++i;
            n_word += word[i] - 32;
            continue;
        }
        if ((c > 96 && c < 123) || (c > 223 && c < 256)) n_word += (word[i] -32);
        else n_word += word[i];
    }
    return n_word;
}


template <class STRG_TYPE>
STRG_TYPE capitalize_first(const STRG_TYPE &word)
{
    int l = word.length(), c, i = 0;
    STRG_TYPE n_word;

    if (!l) return n_word;
    c = int(word[0]);
    if ( c < 0 && l > 1) {
        n_word += word[0];
        n_word += word[1] -32;
        ++i;
    }
    else if ((c > 96 && c < 123) || (c > 223 && c < 256)) n_word += word[0] -32;

    ++i;
    for (; i < l; ++i) {
        n_word += word[i];
    }
    return n_word;
}


template <class STRG_TYPE>
STRG_TYPE decapitalize(const STRG_TYPE &word)
{
    int l = word.length(), c;
    STRG_TYPE n_word;

    for (int i = 0; i < l; ++i) {
        c = int(word[i]);
        if (c < 0 && i < l-1) {
            n_word += word[i];
            ++i;
            n_word += word[i] + 32;
            continue;
        }
        if ((c > 64 && c < 91) || (c > 191 && c < 222)) n_word += (word[i] + 32);
        else n_word += word[i];
    }
    return n_word;
}


template <typename T,template <typename> class C>
int isinContainer(const T& value, const C<T>& G)
{
	for (const auto& v : G) {
		if (v == value) return 1;
	}
	return 0;
}

}

#endif // VAL_UTILS_H_INCLUDED
