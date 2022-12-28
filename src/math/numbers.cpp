
#include <numbers.h>
#include <modq.h>
#include <rand.h>
#include <rational.h>
#include <pol.h>
#include <Glist.h>
#include <thread>
#include <pol_arithmetic.h>
#include <modint.h>
#include <fstream>

namespace val
{

Glist<int> getPrimelist(const std::string &path)
{
    val::Glist<int> Prim;
    int p;
    std::ifstream file (path,std::ios::in);
    while (file) {
        file>>p>>p;
        Prim.push_back(p);
    }
    if (Prim.isempty()) {  //generate primes
		int pp=100000;
		for (int i=0 ;i<1000;++i) {
			pp=val::nextprime(pp);
			Prim.push_back(pp);
		}
    }
    return Prim;
}


// Computes greatest nat. number  m, with m<=sqrt(n/2), by bisection-method in [a,b].
val::integer sqrtofhalf(const val::integer& n)
{
 using namespace val;
 integer a,a1,c,d,b=n;

 while ((d=(b-a))>1) {
	 a1=a+d/integer(2);
	 c=a1*a1;
	 if (integer(2)*c>n) b=std::move(a1);
	 else a=std::move(a1);
 }
 if (integer(2)*(b*b)>n) return a;
 else return b;
}


int isprime(int n)
{
 int i;

 if (n<0) n*=-1;
 if (n<2) return 0;
 if (n==2) return 1;
 if (n>2 && (n%2==0)) return 0;

 for (i=2;(i*i)<=n;i++)
     if (n%i==0) return 0;
 return 1;
}


int ispseudoprime(int n,int steps)
{
 int s,t,i,j,prim;
 RandomClass RandGen;

 if (n<0) n*=-1;
 if (n<=1) return 0;
 if (n==2) return 1;
 if (n%2==0) return 0;

 modq::q=n;
 modq a,eins(1),minuseins(-1);

 t=n-1;
 for (s=0;(t%2)==0;s++) t/=2;    // ==> n-1 = 2^s * t , t odd

 for (j=0;j<steps;j++) {
     prim=0;
     a= modq(RandGen.random(2,n-1));
     a=power(a,t);
     if ((a==eins) || (a==minuseins)) {
         prim=1;
         continue;
     }
     for (i=1;i<s;i++) {
        a*=a;
        if (a==minuseins) {
            prim=1;
            break;
        }
     }
     if (prim==0) return 0;
 }
 return prim;
}



int nextprime(int n)               // depends on isprime
{
 const int grenze=2147483647;
 int i,s;

 if (n<2) return 2;
 if (n==2) return 3;
 if (n%2) s=2;
  else s=1;
 for (i=n+s;i<grenze;i+=2)
     if (ispseudoprime(i)) {
         if (isprime(i)) return(i);
     }
 return 2;
}

// Next prime > n;
val::integer nextprime(const val::integer &n1,int schritte)
{
    val::integer eins(1),zwei(2),n;
    if (n1.iszero() || n1==eins || n1 == -eins) return zwei;
    n=n1;
    if (n.signum()<0) n.changesign();
    if (n.iseven()) n++;
    while (!pseudoprim(n,schritte)) {
          n+=zwei;
    }
    return n;
}


// Random number between 0 and n-1
val::integer randinteger(const val::integer &n)
{
    int l=n.abslength();
    val::integerRandClass rnd;
    val::integer a = rnd.integerrandom(l);

    if (a.signum()<0) a.changesign();
    if (a.abslength()>=l) a%=n;
    return a;
}


val::integer multmod(const val::integer& a,const val::integer& b,const val::integer& basis)
{
 return (a*b)%basis;
}


val::integer hochmod(const val::integer &a1,const val::integer& basis,const val::integer &n1)
{
 val::integer x(1),zwei(2),eins(1),a,n;

 if (a1==eins) return x;
 a=a1;n=n1;

 while (!n.iszero()) {
       if (!n.iseven()) x=multmod(x,a,basis);

       //a*=a;
       //if (a>basis) a%=basis;
       a=multmod(a,a,basis);
       n/=zwei;
 }
 return x;
}


//Explanation in numbers.h
int pseudoprim(const val::integer &n1,int schritte)
{
 int j,prim;
 val::integer i,t,s,a,m,eins(1),zwei(2),drei(3),n;

 if (n1==zwei || n1 == -zwei) return 1;
 if (n1==drei || n1 == -drei) return 1;
 if (n1.iseven()) return 0;
 if (n1==eins || n1==-eins) return 0;
 n=n1;
 if (n.signum()<0) n.changesign();


 m=t=n-eins;
 for (s=0;t.iseven();s++) t/=zwei;    // ==> n-1 = 2^s * t , t odd

 for (j=0;j<schritte;j++) {
     prim=0;
     a= randinteger(n-drei)+zwei;
     a=hochmod(a,n,t);
     if ((a==eins) || (a==m)) {
        prim=1;
        continue;
     }
     for (i=1;i<s;i++) {
        a=multmod(a,a,n);
        if (a==m) {
            prim=1;
            break;
        }
     }
     if (prim==0) return 0;
 }
 return prim;
}

integer simultcong(const d_array<integer> &z,const d_array<integer> &m,const d_array<integer> &M)
{
 int i,n=z.length();
 integer x=1,y=0,x1=0;
 integer z1;

 d_array<integer> prod;

 if (z.isempty() || z.length()!=m.length()) return z1;

 z1=z[0];
 if (M.isempty())  {
        prod = d_array<integer>(n);
        prod[0]=m[0];
        for (i=1;i<n;i++) prod[i] = prod[i-1]*m[i];
 }
 for (i=1;i<n;i++){ 
	 if (!M.isempty()) {
        euclid(M[i-1],m[i],x,y);
        z1= x*M[i-1]*z[i] + y*m[i]*z1;
        z1%=M[i];

	 }
	 else {
        euclid(prod[i-1],m[i],x,y);
        z1= x*prod[i-1]*z[i] + y*m[i]*z1;
        z1%=prod[i];
	 }
 }

 if (z1.signum()<0) {
    if (!M.isempty()) z1+=M[n-1];
    else z1+=prod[n-1];
 }
 return z1;
}


int ratconstruction(const val::integer &m,const val::integer &a,const val::integer &k,val::integer &r,val::integer &t)
{
 integer v1=m,v2,z1,z2,q;

 r=a;t=integer(1);

 while (r>k) {
    q=v1/r;
    z1=v1 - q*r;
    z2=v2 - q*t;
    v1=std::move(r);
    v2=std::move(t);
    r=std::move(z1);
    t=std::move(z2);
 }
 z1=val::ggTspez(r,t);
 if (z1!=integer(1) && z1 != integer(-1)) {return 0;}
 if (t.signum()==-1) {t.changesign();r.changesign();}
 if (t>=k) return 0;
 return 1;
}


int ratconstruction(const val::integer &m,const val::integer &a,const val::integer &N,const val::integer &D,val::integer &r,val::integer &t)
{
 integer v1=m,v2,z1,z2,q;

 r=a;t=integer(1);

 while (r>N) {
    q=v1/r;
    z1=v1 - q*r;
    z2=v2 - q*t;
    v1=std::move(r);
    v2=std::move(t);
    r=std::move(z1);
    t=std::move(z2);
 }
 z1=val::ggTspez(r,t);
 if (z1!=integer(1) && z1 != integer(-1)) {return 0;}
 if (t.signum()==-1) {t.changesign();r.changesign();}
 if (t>D) return 0;
 return 1;
}


int hensel_lift(const pol<integer> &f,const modq &x,int n,integer &s)
{
    if (n<1) return 0;
    integer q(modq::q),a(x),b;
    s=a;

    if ((f(s)%q)!=integer(0)) return 0;

    modq div(int(f.derivation(s)%q));

    if (div==modq(0)) return 0;

    for (int i=2;i<=n;++i) {
        b = f(s)/q;
        b= b % integer(modq::q);
        a = integer(modq(int(-b))/div);
        s+= a*q;
        q*=integer(modq::q);
    }

    return 1;
}

val::pol<val::integer> simultanpolynomial(const val::d_array<val::pol<val::integer> > &Fint,const val::d_array<val::integer> &m,const val::d_array<val::integer>& M)
{
    using namespace val;
    int i,n=Fint.length(),deg;
    d_array<integer> z(n);
    integer Z;
    pol<integer> f;

    deg=-1;
    for (i=0;i<n;i++)
        if (Fint[i].degree()>deg) deg = Fint[i].degree();

    f.insert(integer(1),deg);
    deg--;
    for (;deg>=0;deg--) {
        for (i=0;i<n;i++) z[i] = (Fint[i])[deg];
        Z=simultcong(z,m,M);
        if (!Z.iszero())  f.insert(Z,deg);
    }
    return f;
}

val::pol<val::rational> rationalconstruction(const val::pol<val::integer> &fint,const val::integer& M)
{
    using namespace val;

    integer k=sqrtofhalf(M),r,t;
    pol<rational> f;
    polIterator<integer> itfint;

    if (fint.iszero()) return f;

    f.insert(rational(1),fint.degree());
    itfint=fint;
    itfint++;
    for (;itfint;itfint++) {
        if (!ratconstruction(M,itfint.actualcoef(),k,r,t)) {
            f.del();
            return f;
        }
        f.insert(rational(r,t),itfint.actualdegree());
    }

    return f;
}

// --------------------- modular gcd --------------------------------------------------------------------------

val::pol<val::modint> Tomodintpol(const val::pol<val::integer>& f,int q)
{
    val::integer iq(q);
    val::pol<val::modint> g;


    for ( const auto& t : f) {
        g.insert(val::modint(int(t.actualcoef()%iq),q),t.actualterm());
    }

    g.normalize();
    return g;
}


void ComputegcdInt(const val::pol<val::integer> &f,const val::pol<val::integer> &g,val::pol<val::integer> &gcdInt,int q)
{
    val::pol<val::modint> fmod(Tomodintpol(f,q)),gmod(Tomodintpol(g,q)),h;

    h = val::gcd(fmod,gmod);
    gcdInt = val::modPolToIntPol(h);
}


val::pol<val::rational> modular_gcd(const val::pol<val::rational> &f, const val::pol<val::rational> &g,
                                    const val::Glist<int> &Prim = val::Glist<int>(),int n,int multithread)
{
    using namespace val;
    if (f.iszero()) return g;
    if (g.iszero()) return f;

    if (n<1) n=1;

    int i=0,j,deg=-1,degvalid;
    integer zero;
    pol<rational> h,hold,onepol(rational(1)),zeropol;
    pol<integer> fint = val::primitivpart(f), gint = val::primitivpart(g), hint;
    Glist<int> Primelist;
    d_array<pol<integer> > Hint(n+1);
    d_array<integer> M(n+1),m(n+1);
    d_array<int> q(n);
    d_array<std::thread*> t(n);

    m[0]=M[0] = integer(1);

    const Glist<int> *Primlist=nullptr;

    if (Prim.isempty()) {
        Primelist = getPrimelist();
        Primlist = &Primelist;
    }
    else Primlist = &Prim;

    auto p = Primlist->begin();

    for (i=0;p!=Primlist->end();i++) {
        for (j=0;j<n;j++) {
            while (p!=Primlist->end()) {
                q[j] = p();
                if ((fint.LC()%integer(q[j])!=zero) && (gint.LC()%integer(q[j])!=zero) ) break;
                ++p;
            }
            if (p==Primlist->end()) {h.del();break;}
            m[j+1] = integer(q[j]);
            M[j+1] = M[j] * integer(q[j]);
            ++p;
		}
        if (!multithread) {
            for (j=0;j<n;++j) ComputegcdInt(fint,gint,Hint[j+1],q[j]);
        }
        else {
            for (j=0;j<n;++j) t[j] = new std::thread(ComputegcdInt,std::cref(fint),std::cref(gint),std::ref(Hint[j+1]),q[j]);
            for (j=0;j<n;++j) t[j]->join();
            for (j=0;j<n;++j) delete t[j];
        }
        if (deg==-1) deg = Hint[1].degree();
        degvalid=1;
        for (j=1;j<=n;++j) {
            if (Hint[j].degree()==0) return onepol;
            if (Hint[j].degree()<deg) {
                deg = Hint[j].degree();
                degvalid=0;
            }
            if (Hint[j].degree() > deg) degvalid=0;
        }
        if (!degvalid) {
            hold.del();
            m[0] = M[0] = integer(1);
            Hint[0].del();
            continue;
        }

        hint=simultanpolynomial(Hint,m,M);
        M[0] = std::move(M[n]);
        m[0] = M[0];
        Hint[0] = std::move(hint);
        h=rationalconstruction(Hint[0],M[0]);
        if (!h.iszero() && h==hold) {
            if ((f%h==zeropol) && (g%h==zeropol)) break;
        }
        hold=std::move(h);
    }

    return h;
}

// -----------------------------------------------------------------------------------------------------------------------------

} //end namespace val
