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

template <class T>
int replace(std::basic_string<T> &s, const std::basic_string<T> &from, const std::basic_string<T> &to);


DLL_PUBLIC char* StringToChar(const std::string& s);

#ifdef __WIN32
DLL_PUBLIC int system(std::string &command);
#else
// wrapper for std::system
DLL_PUBLIC int system(const std::string &command, int silent = 1);
#endif


// Checks if valu is in container G
template <typename T,template <typename> class C>
int isinContainer(const T& value, const C<T>& G);

// Crates a string Container of type G, of the single words in sf, separated by a list of separators (e.g {',', ';'})
// if emptywords = 0, empty words are ignored, container ignore, consists of char, that will be ignored.
template <class T,template <typename> class C, template <typename> class G>
G<std::basic_string<T>> getwordsfromstring(const std::basic_string<T> &sf,const C<T>& separators,int emptywords = 0, const C<T> &ignore = C<T>());


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

    m = Min(m, n);

    for (int i = n - m; i < n; ++i) {
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


// replace every appearende of 'from' to 'to' in string s
template <class T>
int replace(std::basic_string<T> &s, const std::basic_string<T> &from, const std::basic_string<T> &to)
{
    int found = 0;
    if (from.empty()) return found;

    size_t pos = 0, n = to.length(), m = from.length();
    while ((pos = s.find(from,pos)) != std::basic_string<T>::npos) {
        s.replace(pos,m,to);
        pos += n;
        found = 1;
    }
    return found;
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
        if (s[i] == '-' || s[i]  == ',' || s[i] =='.' || s[i] == 'e' || s[i] == '/' || s[i] == ')' || s[i] == '(' ) continue;
        else if (s[i]<48 || s[i] > 57) return 0;
    }
    return 1;
}


template <class STRG_TYPE>
int isfloatnumber(const STRG_TYPE &s)
{
    int l = s.length();
    if (!l) return 0;

    for (int i = 0; i < l; ++i) {
        if (s[i] == '-' || s[i] =='.' || s[i] == 'e') continue;
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


// returns natural number found from s[i]. shifts i to the next char that is not 0,...,9.
template <class STRG_TYPE>
STRG_TYPE findnumber(const STRG_TYPE &s, int &i)
{
    int n=s.size(),e_set=0,p_set=0;
    STRG_TYPE out;

    for (;i<n;i++) {
        if (s[i]>='0'  && s[i]<='9') out+=s[i];
        else if (s[i]=='e' && i<n-1 && s[i+1]!='x') {
            if (e_set) {i++; break;}
            if (i==n-1) {i++;break;}
            if (s[i+1] == '-' || s[i+1]=='+') {
                out+=s[i];i++;out+=s[i];e_set=1;
            }
            else {out+=s[i];e_set=1;}
        }
        else if (s[i]==',' || s[i]=='.') {
            if (p_set) {i++;break;}
            out+=s[i];
            p_set=1;
        }
        else {break;}
    }

    return out;
}



template <typename T,template <typename> class C>
int isinContainer(const T& value, const C<T>& G)
{
    for (const auto& v : G) {
        if (v == value) return 1;
    }
    return 0;
}


template <class T,template <typename> class C, template <typename> class G>
G<std::basic_string<T>> getwordsfromstring(const std::basic_string<T> &sf,const C<T>& separators,int emptywords, const C<T> &ignore)
{
    G<std::basic_string<T>> values;
    std::basic_string<T> s="";
    int n = sf.length();

    for (int i = 0; i < n ; ++i) {
        if (val::isinContainer(sf[i],ignore)) continue;
        if (val::isinContainer(sf[i],separators)) {
            if (emptywords || s != "") values.push_back(s);
            //else if (s != "") values.push_back(s);
            s = "";
        }
        else s += sf[i];
    }

    if (emptywords || s!= "") values.push_back(s);

    return values;
}

}

#endif // VAL_UTILS_H_INCLUDED
