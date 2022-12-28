#ifndef NUMBERS_H_INCLUDED
#define NUMBERS_H_INCLUDED

#include <d_array.h>
#include <val_basics.h>

namespace val
{

class modq;
class integer;
class rational;
template <class T> class pol;
template <class T> class Glist;



DLL_PUBLIC val::Glist<int> getPrimelist(const std::string &path=primlistpath);

//modq power(modq,int);

// naive :
DLL_PUBLIC int isprime(int n);


// Probabilistic function to check if a number is pseudo-prime.
// If 1 is returned, so the probability of n being prime is 1-(3/4)^steps.
// If 0 is returned, n is definitely not prime.  
DLL_PUBLIC int ispseudoprime(int n,int steps=5);

DLL_PUBLIC int pseudoprim(const val::integer &n1,int schritte=5);

// Computes least prime p > n.
DLL_PUBLIC int nextprime(int n);               // needs  isprime

DLL_PUBLIC val::integer nextprime(const val::integer &n1,int schritte=5);

// simultaneous congruence: z1= z[i] mod m[i]  (Chinese Remainder Theorem)
// M[i] = m[0]*...*m[i] for 0<=i<n
DLL_PUBLIC integer simultcong(const d_array<integer> &z,const d_array<integer> &m,const d_array<integer> &M = d_array<integer>());

// Rational Construction: Computes to a in Z/mZ, r,t with |r| < k, 0 <= t < k, so that r = t*a mod m, where k in {1,...,m}.
// If r*t^(-1) = a in Z/mZ, 1 is returned; 0 otherwise.
DLL_PUBLIC int ratconstruction(const val::integer &m,const val::integer &a,const val::integer &k,val::integer &r,val::integer &t);

// Computes to a in Z/mZ r,t with |r|<=N, 0<=t<=D, so that r=t*a mod m, N,D in {1,...,m}
// If r*t^(-1) = a in Z/mZ, 1 is returned; 0 otherwise.
DLL_PUBLIC int ratconstruction(const val::integer &m,const val::integer &a,const val::integer &N,const val::integer &D,val::integer &r,val::integer &t);


// Lift root x of f mod p up to root  s of f mod p^n.
// Return: 1 if success; 0 otherwise.
DLL_PUBLIC int hensel_lift(const pol<integer> &f,const modq &x,int n,integer &s);

DLL_PUBLIC val::pol<val::integer> simultanpolynomial(const val::d_array<val::pol<val::integer> > &Fint,const val::d_array<val::integer> &m,const val::d_array<val::integer> &M);

DLL_PUBLIC val::pol<val::rational> rationalconstruction(const val::pol<val::integer> &fint,const val::integer& M);

DLL_PUBLIC val::pol<val::rational> modular_gcd(const val::pol<val::rational> &f, const val::pol<val::rational> &g,
                                    const val::Glist<int> &Prim,int n=3,int multithread=1);
} // end namespace val


#endif // NUMBERS_H_INCLUDED
