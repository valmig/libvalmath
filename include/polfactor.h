#ifndef POLFACTOR_H_INCLUDED
#define POLFACTOR_H_INCLUDED

// Factorization of univariate polynomials 

#include <s_modinteger.h>
#include <modq.h>
#include <d_array.h>
#include <Glist.h>
#include <rational.h>
#include <pol.h>


namespace val
{

extern DLL_PUBLIC Glist<int> Primes;
//template <class T> class pol;

DLL_PUBLIC pol<modq> findirreducible(int n);

DLL_PUBLIC int polfactor(const pol<modq> &f,d_array<pol<modq>> &afaktor,d_array<int> &amult);

// computes primitive irreducible factors of g (therefore the factors have integer coefficients) in Q[x]
DLL_PUBLIC d_array<pol<rational>> polfactor(const pol<rational> &g,int comment=0);

DLL_PUBLIC d_array<rational> rational_roots(const pol<rational> &g,const Glist<int> &Primelist);

DLL_PUBLIC d_array<rational> rational_roots(const pol<rational> &g);

} // end namespace val




#endif // POLFACTOR_H_INCLUDED
