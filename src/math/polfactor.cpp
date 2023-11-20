
#include <polfactor.h>
#include <rand.h>
#include <s_modinteger.h>
#include <val_basics.h>
#include <pol_arithmetic.h>
#include <fstream>
#include <thread>
#include <modint.h>
#include <numbers.h>



std::string primlistpath = val::primlistpath;

template <class T>
void subsets(const val::Glist<T>& G,val::Glist<val::Glist<T>>& S,int m,int i0=0)
{
	int i,n=G.length();
	if ((m==0) || m>n-i0) {
		val::Glist<T> A;
		S.push_back(A);
		return;
	}
	else if (m==n-i0) {
		val::Glist<T> A;
		for (i=i0;i<n;++i) A.push_back(G[i]);
		S.push_back(A);
		return;
	}
	else if (m==1) {
		for (i=i0;i<n;++i) {
			S.push_back(val::Glist<T>{G[i]});
		}
	}
	else {
		val::Glist<val::Glist<T>> R;
		const T &a=G[i0];
		subsets(G,R,m-1,i0+1);
		for (auto& A : R) A.push(a);
		S.append(std::move(R));
		subsets(G,S,m,i0+1);
		return;
	}
}

template <class T>
int inSet(const val::Glist<T> &G,const T &n)
{
    for (const auto& i : G) {
        if (i==n) return 1;
    }
    return 0;
}

// modular gcd - functions ---------------------------------------------------------------------------------------------

template <class T>
val::pol<val::integer> GetFint(const val::pol<T> &fmodq)
{
    using namespace val;
    pol<integer> f;
    polIterator<T> itfmodq;

    for (itfmodq=fmodq;itfmodq;itfmodq++) f.insert(integer(itfmodq.actualcoef()),itfmodq.actualdegree());
    return f;
}


// ---------------------------------------------------------------------------------------------------------------------


val::pol<val::modq> euclid(val::pol<val::modq> a0,val::pol<val::modq> a1,val::pol<val::modq>& x0,val::pol<val::modq>& y0);

struct henselnode
{
    val::pol<val::s_modinteger> f,s,t;
};


int log2(int z)
{
    if (!z) return 0;
    int i,l,numberofbits=8*sizeof(int);
    unsigned  m = 1<<(numberofbits -1); // Highest bit 1 the rest 0;

    l=numberofbits-1;
    for (i=l;i>=0;--i,m>>=1,--l) {
        if (m&z) break;
    }
    for (--i,m>>=1;i>=0;--i,m>>=1) {
        if (m&z) return l+1;
    }
    return l;
}

int logi(int p,const val::integer &B)
{
    val::integer q(1),P(p),l;
    if (B<=q) return 0;
    l=1;
    do {
        q*=P;
        if (B<=q) return l;
        ++l;
    }
    while (1);
    return l;
}


// pol<integer> help functions ---------------------------------------------------------------------------------------------------------------

val::integer MaxNorm(const val::pol<val::integer> &f)
{
    val::integer m;
    for (const auto& pf : f) m = val::Max<val::integer>(m,val::abs(pf.actualcoef()));
    return m;
}


val::integer Norm1(const val::pol<val::integer> &f)
{
    val::integer m;
    for (const auto& pf : f) m+=val::abs(pf.actualcoef());
    return m;
}

// ------------------------------------------------------------------------------------------------------------------------------

// pol<s_modinteger> - help functions  ------------------------------------------------------------------------------------------------

val::pol<val::s_modinteger> TomodintegerPol(const val::pol<val::integer> &f)
{
    const val::integer &q = val::s_modinteger::q;
    val::pol<val::s_modinteger> g;

    for (const auto& pf : f) g.insert(pf.actualcoef()%q,pf.actualterm());

    return g;
}


val::pol<val::s_modinteger> TomodintegerPol(const val::pol<val::modq> &f)
{
    val::pol<val::s_modinteger> g;

    for (const auto& pf : f) g.insert(val::s_modinteger(int(pf.actualcoef())),pf.actualterm());

    return g;
}



val::pol<val::integer> TointegerPol(const val::pol<val::s_modinteger> &f)
{
    val::pol<val::integer> g;
    val::integer m = val::s_modinteger::q/val::integer(2),c;

    for (const auto& pf : f) {
        c = pf.actualcoef();
        if (c>m) c-=val::s_modinteger::q;
        g.insert(std::move(c),pf.actualterm());
    }

    return g;
}

val::pol<val::modq> TomodqPol(const val::pol<val::s_modinteger> &f)
{
    val::pol<val::modq> g;

    for (const auto& pf : f) g.insert(val::modq(int(val::integer(pf.actualcoef()))),pf.actualterm());

    return g;
}


void hensellift(const val::pol<val::s_modinteger>& f,val::pol<val::s_modinteger> &g,val::pol<val::s_modinteger> &h,
                val::pol<val::s_modinteger> &s,val::pol<val::s_modinteger> &t)
{
    val::pol<val::s_modinteger> e,q,r,b,c,d;

    e = f-g*h;
    val::divrest(s*e,h,q,r);
    if (h*q+r != s*e) {std::cout<<"\n Error on divrest!\n"<<h*q+r -s*e;}
    g+= t*e + q*g; h+=r;
    b = s*g+t*h - val::pol<val::s_modinteger>(val::s_modinteger(1));
    val::divrest(s*b,h,c,d);
    if (h*c + d != s*b) {std::cout<<"\n Error on divrest!\n"<<h*c+d - s*b;}
    s-=d; t-=(t*b + c*g);
}


val::d_array<val::d_array<henselnode>> createhenseltree(const val::d_array<val::pol<val::modq>> &factor)
{
    int m = log2(factor.length()) + 1,r,i,j,even,k;
    val::d_array<val::d_array<henselnode>> F(m);
    val::pol<val::modq> s,t,h;
    val::modq a;

    F[m-1].reserve(r = factor.length());
    for (i=0;i<r;++i) F[m-1][i].f = TomodintegerPol(factor[i]);

    for (i=m-2;i>=0;--i) {
        even=1;
        k = F[i+1].length();
        r = k/2;
        if (k%2) {++r; even=0;}
        F[i].reserve(r);
        for (j=0;j<k-1;j+=2) {
            F[i][j/2].f = F[i+1][j].f * F[i+1][j+1].f;
            h=euclid(TomodqPol(F[i+1][j].f),TomodqPol(F[i+1][j+1].f),s,t);
            a = inv(h.LC());
            s*=a; t *=a;
            F[i][j/2].s= TomodintegerPol(s);
            F[i][j/2].t= TomodintegerPol(t);
        }
        if (!even) {
            F[i][r-1].f = F[i+1][k-1].f;
        }
    }

    return F;
}


void hensellift(val::d_array<val::d_array<henselnode>> &F,int l)
{
    int d = log2(l),k,i,j,m=F.length(),n,even;

    std::cout<<"\n d = "<<d<<std::flush;

    for (k=0;k<d;++k) {
        val::s_modinteger::q*=val::s_modinteger::q;
        for (i=1;i<m;++i) {
            n=F[i].length();
            even=1;
            if (n%2) even=0;
            for (j=0;j<n-1;j+=2) {
                hensellift(F[i-1][j/2].f,F[i][j].f,F[i][j+1].f,F[i-1][j/2].s,F[i-1][j/2].t);
            }
            if (!even) F[i][n-1] = F[i-1][F[i-1].length() -1];
        }
        std::cout<<"\n factors modulo "<<val::s_modinteger::q;
        for (i=0;i<F[m-1].length();++i) {
            std::cout<<std::endl<<F[m-1][i].f;
            if (!(F[0][0].f%F[m-1][i].f).iszero()) {
                std::cout<<"\n Not a factor!";
            }
            else std::cout<<"\n Is factor!";
        }
    }
}

val::d_array<val::pol<val::s_modinteger>> factorlift(const val::pol<val::integer>& f,const val::d_array<val::pol<val::modq>> &factormodq,int l)
{
    using namespace val;
    s_modinteger::q=integer(modq::q);
    int i,r=factormodq.length(),m,d,k,j,n,even;
    d_array<pol<s_modinteger>> factor(r);
    modq a1 = inv(modq(int(f.LC()%integer(modq::q))));
    s_modinteger a = s_modinteger(int(a1));

    if (l==1) {
        for (i=0;i<r;++i) factor[i] = TomodintegerPol(factormodq[i]);
        return factor;
    }

    d_array<d_array<henselnode>> F=createhenseltree(factormodq);

    m=F.length();
    d=log2(l);

    // Hensellift:

    for (k=0;k<d;++k) {
        val::s_modinteger::q*=val::s_modinteger::q;
        a *= (s_modinteger(2) - s_modinteger(f.LC())*a);
        F[0][0].f = TomodintegerPol(f);
        F[0][0].f*=a;
        for (i=1;i<m;++i) {
            n=F[i].length();
            even=1;
            if (n%2) even=0;
            for (j=0;j<n-1;j+=2) {
                hensellift(F[i-1][j/2].f,F[i][j].f,F[i][j+1].f,F[i-1][j/2].s,F[i-1][j/2].t);
            }
            if (!even) F[i][n-1] = F[i-1][F[i-1].length() -1];
        }
    }

    for (i=0;i<r;++i) factor[i] = std::move(F[m-1][i].f);

    return factor;
}



//

// pol<modq> - help functions ---------------------------------------------------------------------------------------------------------------

// => gcd (a,b) = euclid (a,b,x0,y0)  ( = x0*a + y0*b )
val::pol<val::modq> euclid(val::pol<val::modq> a0,val::pol<val::modq> a1,val::pol<val::modq>& x0,val::pol<val::modq>& y0)
{
 val::pol<val::modq> h,q,x1,y1(val::modq(1));

 x0=val::pol<val::modq>(val::modq(1));y0 = h; // => y0 = 0;
 while (!a1.iszero()) {
       q=a0/a1;
       h=a1;
       a1=a0%a1;
       a0=std::move(h);
       h= x1;
       x1=x0 - q*x1;
       x0= std::move(h);
       h= y1;
       y1= y0 - q*y1;
       y0= std::move(h);
 }
 return a0;
}


// Creates random polynomial of degree n
val::pol<val::modq> randpoly(int n)
{
 int i;
 val::pol<val::modq> h;
 val::RandomClass rand;


 h.insert(val::modq(1),n);            // leading coef = 1

 for (i=0;i<n;i++) h.insert(val::modq(rand.random(val::modq::q)),i);

 return h;
}


val::pol<val::modq> spur2(int m,val::pol<val::modq>& g,const val::pol<val::modq>& f)
{
 int i;
 val::pol<val::modq> x,h;

 x=g;
 h=x;
 for (i=1;i<m;i++) {
     x*=x;
     if (x.degree()>=f.degree()) x%=f;
     h+=x;
 }
 return h;
}


// Computes (a^n mod basis) non-recursive.
template <class T>
val::pol<T> potenzmod(val::pol<T> a,val::integer n,const val::pol<T>& basis)
{
 val::pol<T> x(T(1),0);
 val::integer zwei(2);

 while (n!=0) {
      if (!n.iseven()) {
         x*=a;
         if (x.degree()>=basis.degree()) x%=basis;
       }
       a*=a;
       if (a.degree()>=basis.degree()) a%=basis;
       n/=zwei;
 }
 return x;
}



// Computes non-trivial divisors of f.
// Prem.: f=f1*...*fr where fi different, irreducible of deg. d<n
val::pol<val::modq> splitpol(const val::pol<val::modq>& f,int d)
{
 val::integer exponent;
 int k,n=f.degree();
 val::pol<val::modq> g,one(val::modq(1)),h;
 val::RandomClass rand;

 exponent=(val::power(val::integer(val::modq::q),d)-val::integer(1))/val::integer(2);
 do {
     k=rand.random(n-1)+1;
     g=randpoly(k);
     h = val::gcd(f,g);
     if (h!=one) return h;
     if (val::modq::q!=2) g=potenzmod(g,exponent,f);
     else g=spur2(d,g,f);
     g=val::gcd(g-one,f);
 }
 while((g==one) || (g==f));

 return g;
}


int equaldegree(const val::pol<val::modq>& f,int d,val::d_array<val::pol<val::modq>> &apol,int i=0)
{
  val::pol<val::modq> g;

  if (f.degree()==d) {
     apol[i++]=f;
     return i;
  }
  g=splitpol(f,d);
  i=equaldegree(f/g,d,apol,i);
  i=equaldegree(g,d,apol,i);
  return i;
}

// --------------------------------------------------------------------------------------------------------------------------------------


//  -----------------------------------------------------------------------------------------------------------------------------------
namespace rationalroots
{
// Least  nat. number l >= 1 with p^l >=B
int lowestpowergreater(const val::integer &B,const val::integer &p,val::integer &q)
{
    int l=1;
    q=val::integer(p);

    while (q<B) {
        q*=p;
        ++l;
    }
    return l;
}


// Computes all roots of f in modq.
// Return : 0 if double roots exists, otherwise 1.
int polynomzeros(const val::pol<val::modq> &f,val::d_array<val::modq> &zeros)
{
    val::pol<val::modq> g,h;
    val::modq zero(0);

    zeros.del();

    h.insert(val::modq(1),1);
    h = potenzmod(h,val::modq::q,f); // h = x^p mod f
    h.insert(val::modq(-1),1);

    g=val::gcd(f,h);

    if (g.degree()==0) return 1;
    if (val::deg(val::gcd(g,g.derive()))!=0) return 0;
    zeros.reserve(g.degree());
    for (int i =0 ; i<val::modq::q;++i) {
        if (g(val::modq(i)) == zero) zeros.push_back(val::modq(i));
    }
    return 1;
}

}  // end namespace rationalroots

namespace val
{

val::Glist<int> Primes = val::getPrimelist();

pol<modq> findirreducible(int n)
{
 int i,m,fertig;
 pol<modq> f,g,x(modq(1),1),one(modq(1));
 integer q(modq::q);

 if (n==1) return randpoly(1);
 m=n/2;

 fertig=0;
 while (!fertig) {
       f=randpoly(n);
       if (f[0]==val::modq(0)) f.insert(val::modq(1),0);
       g=potenzmod(x,q,f);
       for (i=1;i<=m;i++)
	   if (val::gcd(g-x,f)!=one) {
	      fertig=0;
	      break;
	   }
	   else {
	      fertig=1;
	      g=potenzmod(g,q,f);
	   }
 }
 return f;
}


int polfactor(const pol<modq> &f,d_array<pol<modq>> &afaktor,d_array<int> &amult)
{
 int e,i,j,n=f.degree(),anzahl,s;
 modq lead;
 pol<modq> x,h(modq(1),1),f0,g,onepol(modq(1)),zeropol;
 d_array<pol<modq>> agi;
 integer q(modq::q);

 afaktor.del(); amult.del();
 // => h = x

 if (n<=1) {
    afaktor= val::d_array<val::pol<val::modq>>(1);
    amult= val::d_array<int>(1);
    afaktor[0]=f;
    amult[0]=1;
    return 1;
 }

 f0 = f;
 f0.normalize();
 x = h;

 afaktor.reserve(n);
 amult.reserve(n);
 i=0;
 anzahl=0;

 do {
    i++;
    h=potenzmod(h,q,f0);
    g=val::gcd(h-x,f0);
    if (g!=onepol) {
       s=g.degree()/i;
       agi=val::d_array<val::pol<val::modq>>(s);
       equaldegree(g,i,agi);
       for (j=0;j<s;j++) {
			e=0;
			while ((f0%agi[j])==zeropol) {
				f0/=agi[j];
				e++;
			}
			afaktor[anzahl]=agi[j];
			amult[anzahl++]=e;
       }
    }
 }
 while (f0.degree()>1);

 if (f0.degree()>=1) {
    afaktor[anzahl]=f0;
    amult[anzahl++]=1;
 }
 return anzahl;
}


d_array<pol<rational>> polfactor(const pol<rational> &g,int comment)
{
    d_array<pol<rational>> afaktor;
    if (g.degree()<=1) {
        afaktor.reserve(1);
        afaktor[0] = toRationalPolynom(primitivpart(g));
        if (afaktor[0].LC().signum()==-1) afaktor[0]*=rational(-1);
        return afaktor;
    }
    pol<integer> f = primitivpart(g/modular_gcd(g,g.derive(),Primes)),f1,f2;
    pol<rational> h,zeroratpol;

    if (f.degree()==1) {
        afaktor.reserve(1);
        afaktor[0] = toRationalPolynom(f);
        return afaktor;
    }
    if (f.LC().signum()==-1) f*=integer(-1);
    
    int divbyx = 0;
    pol<rational> X(rational(1),1);
    
    if (f[0] == integer(0)) {
		divbyx = 1;
		f.getdivbypower(1);
	}
    

    int n=f.degree(),p=5,l,i,j,r,found;// ,m
    integer b= f.LC(), A = MaxNorm(f),izero,B;
    s_modinteger b1;
    pol<modq> fmodq;
    rational BR = rational(sqrt(double(n)+1)) * rational(power(integer(2),n)) * rational(A) * rational(b);
    d_array<pol<modq>> factormodq;
    pol<s_modinteger> g1,g2,onepol(s_modinteger(1));
    Glist<int> T;
    Glist<Glist<int>> S;
    d_array<int> amult;

    for (const auto& q : val::Primes) {
        p=q;
        if (b%integer(p)==izero) continue;
        modq::q=p;
        fmodq = val::toModqPolynom(f);
        if (gcd(fmodq,fmodq.derive()).degree()==0) break;
    }

    B = nominator(BR)/denominator(BR);
    l = logi(p,integer(2)*B+integer(1));

    fmodq.normalize();
    polfactor(fmodq,factormodq,amult);
    amult.del();
    r = factormodq.length();
    s_modinteger::q=integer(p);
    for (i=0;i<r;++i) T.push_back(i);
    afaktor.reserve(r + divbyx);

    if (comment) {
        integer::SetOutput_Style(integer::FLOAT);
        std::cout<<"\n r = "<<r;
        std::cout<<"\n B = "<<B<<" , l = "<<l<<std::flush;
        integer::SetOutput_Style(integer::INTEGER);
    }

    d_array<pol<s_modinteger>> F = factorlift(f,factormodq,l);

    b1=s_modinteger(b);

    if (divbyx) afaktor.push_back(X);
    
    for (i=1;2*i<=r;) {
        S.dellist();
        subsets(T,S,i);
        found=0;
        for (const auto& H : S) {
            g1=g2=pol<s_modinteger>(b1);
            for (j=0;j<r;++j) {
                if (inSet(H,T[j])) g1*=F[T[j]];
                else g2*=F[T[j]];
            }
            f1=TointegerPol(g1);
            f2=TointegerPol(g2);
            if (Norm1(f1)*Norm1(f2)<B) {
                found=1;
                for (const auto& k : H) {
                    for (j=0;j<r;j++) {
                        if (T[j]==k) {
                            T.delelement(j);
                            --r;
                            break;
                        }
                    }
                }
                afaktor.push_back(toRationalPolynom(primitivpart(f1)));
                f=primitivpart(f2);
            }
            if (found) break;
        }
        if (!found) ++i;
    }
    afaktor.push_back(toRationalPolynom(f));

    return afaktor;
}


d_array<rational> rational_roots(const pol<rational> &g,const Glist<int> &Primelist)
{
    //using namespace val;
    using namespace rationalroots;
    val::d_array<val::rational> zeros;

    if (g.iszero()) return zeros;

    pol<rational> h(g/gcd(g,g.derive()));
    pol<integer> f = primitivpart(h);

    zeros.reserve(f.degree());

    if (f.degree()==1) {
        zeros.push_back(rational(-f[0],f.leader()));
        return zeros;
    }

    if (f.getlastdeg()!=0) {
        zeros.push_back(rational(0));
        f.getdivbypower(1);
    }

    d_array<modq> modqzeros;
    pol<modq> fmodq;
    integer N = val::abs(f.getlastcoef()), D = val::abs(f.leader()),zero(0);
    int found=0;

    for (const auto &p : Primelist ) {
        if (N%integer(p) == zero || D%integer(p) ==zero) continue;
        modq::q=p;
        fmodq = toModqPolynom(f);
        if (polynomzeros(fmodq,modqzeros)) {found = 1;break;}
    }
    if (!found) return zeros;
    integer B = integer(2)*D*N,s,m,r,t;

    int l = lowestpowergreater(B,integer(modq::q),m);
    rational z, rzero;

    for (const auto &x : modqzeros) {
        if (hensel_lift(f,x,l,s)) {
            if (ratconstruction(m,s,N,D,r,t)) {
                z = rational(std::move(r), std::move(t));
                if (h(z) == rzero) zeros.push_back(rational(std::move(z)));
            }
        }
    }

    return zeros;
}


d_array<rational> rational_roots(const pol<rational> &g)
{
    return rational_roots(g,val::Primes);
}


}  // end namespace val
