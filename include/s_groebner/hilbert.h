#ifndef HILBERT_H_INCLUDED
#define HILBERT_H_INCLUDED

#include <pol.h>
#include <rational.h>
#include <s_polynom.h>
#include <vector.h>
#include <string>


namespace hilbert
{


int isdisjunct(const val::s_expo& a,const val::s_expo& b);

int notdisjunct(const val::s_expo& a,const val::s_expo& b);

// c = lcm(a,b)/b, Prem.: c=0
void kgVdiv(const val::s_expo &a,const val::s_expo &b,val::s_expo &c);

// vector feld is sorted in increasing order wrt. ordtype = -1
int reduce(const val::vector<val::s_expo> &feld,val::vector<char> &iszero);

// vector feld is sorted in increasing order wrt. ordtype = -1
val::vector<val::s_expo> reduce(const val::vector<val::s_expo> &feld);


// Prem.: I is sorted in increasing order and already reduced.
// Computes  Hilbert-numerator (not reduced) of K[X]/I
val::pol<val::integer> Hilbertnum(const val::vector<val::s_expo> &I,int anzI,int MCIcrit=0);


// Computes (c,d) with: f(z)-q(z) = c*z^d + a1*z^(d+1)+...
void getlowestmon(const val::pol<val::integer> &f,const val::pol<val::integer> &g,int &c,int &d);

//Prem.: I already reduced and ordered ( -1 )
val::pol<val::rational> Hilbertpolynomial(const val::vector<val::s_expo> &I, int affin=1);

val::pol<val::rational> Hilbertpolynomial(const std::string &filename,int affin=1);

} //end namespace

#endif // HILBERT_H_INCLUDED
