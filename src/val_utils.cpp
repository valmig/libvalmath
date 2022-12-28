#ifndef VAL_UTILS_CPP_INCLUDED
#define VAL_UTILS_CPP_INCLUDED

#include <val_utils.h>
//#include <rational.h>
//#include <complex.h>
//#include <sstream>

namespace val
{

/*
std::string ToString(const int& a)
{
   std::stringstream ss;//create a stringstream
   ss << a;//add number to the stream
   return ss.str();
}


std::string ToString(const double& a)
{
   std::stringstream ss;//create a stringstream
   ss <<a;//add number to the stream
   return ss.str();
}
*/

char* StringToChar(const std::string& s)
{
   int n=0,i;
   char *c=NULL;


    n=s.length();
    c = new char[n+1];
    for (i=0;i<n;i++) c[i]=s[i];
    if (n>0) c[n]='\0';
    return c;
}


/*
int isinteger(const std::string &s)
{
    int i;
    if (s=="") return 0;
    for (i=0;s[i]!='\0';i++) {
        if (s[i]=='-' && i==0 ) continue;
        else if (s[i]<48 || s[i] > 57) return 0;
    }
    return 1;
}
*/

/*
#ifdef RATION_H

rational char_to_rational(const char* s)
{
   std::stringstream ss;
   rational x;

   ss.str(s);
   ss>>x;
   return x;
}

rational string_to_rational(const std::string& s)
{
   std::stringstream ss;
   rational x;

   ss.str(s);
   ss>>x;
   return x;
}

#endif

#ifdef COMPLEX_H_INCLUDED

val::complex string_to_complex(const std::string& s)
{
   std::stringstream ss;
   val::complex x;

   ss.str(s);
   ss>>x;
   return x;
}

#endif
*/


} // End namespace val
#endif // VAL_UTILS_CPP_INCLUDED
