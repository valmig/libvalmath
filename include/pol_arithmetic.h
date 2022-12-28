#ifndef POL_ARITHMETIC_H_INCLUDED
#define POL_ARITHMETIC_H_INCLUDED

// Header file for some functions on polynomials.

#include <val_basics.h>
#include <val_utils.h>
#include <rational.h>
#include <modq.h>
#include <n_polynom.h>
#include <s_polynom.h>
#include <pol.h>

namespace val
{
class modq;
class modint;
class integer;
class rational;
template <class T> class s_polynom;
template <class T> class n_polynom;
template <class T> class pol;
class s_expo;
class n_expo;



template <class T>
s_polynom<T> To_s_polynom(const n_polynom<T> &f);

template <class T>
n_polynom<T> To_n_polynom(const s_polynom<T> &f);

template <template <typename> class poly>
val::rational content(const poly<val::rational> &g);


template <template <typename> class poly>
void primitivpart(const poly<val::rational> &g,poly<val::integer> &f,val::rational &cont);

template <template <typename> class poly>
poly<val::integer> primitivpart(const poly<val::rational> &g);

template <template <typename> class poly>
val::integer content(const poly<val::integer> &g);

template <template <typename> class poly>
void primitivpart(const poly<val::integer> &g,poly<val::integer> &f,val::integer &cont);

template <template <typename> class poly>
poly<val::integer> primitivpart(const poly<val::integer> &g);

template <template <typename> class poly>
int isintegerpol(const poly<val::rational> &g);

template <template <typename> class poly>
poly<val::rational> toRationalPolynom(const poly<val::integer> &g);


template <template <typename> class poly>
poly<val::modq> toModqPolynom(const poly<val::integer> &g);

template <template <typename> class poly>
poly<val::modint> toModintPolynom(const poly<val::integer> &g,int q);

template <typename T,template <typename> class poly>
poly<double> ToDoublePolynom(const poly<T> &g,const double& eps=1e-9);

template <template <typename> class poly>
std::string MPolToString(const poly<rational> &f);

template <typename T,template <typename> class poly>
std::string PolToString(const poly<T> &f);

template <typename T,template <typename> class poly>
poly<integer> modPolToIntPol(const poly<T> &f);


template <template <typename> class poly>
std::string PolToString(const poly<rational> &f);

template <class T>
pol<T> To_unipol(const n_polynom<T> &f,int i=0);

template <class T>
pol<T> To_unipol(n_polynom<T> &&f,int i=0);

template <class T>
pol<T> To_unipol(const s_polynom<T> &f,int i=0);

template <class T>
pol<T> To_unipol(s_polynom<T> &&f,int i=0);

template <class T>
n_polynom<T> To_n_polynom(const pol<T> &f,int i=0);

template <class T>
s_polynom<T> To_s_polynom(const pol<T> &f,int i=0);

template <class T>
int isunivariate(const val::s_polynom<T>& f);

#ifdef N_POLYNOM_H_INCLUDED

#ifdef POL_H
template <class T>
pol<T> To_unipol(const n_polynom<T> &f,int i)
{
    pol<T> g;

    for (const auto &pf : f) {
        g.insert(pf.actualcoef(),pf.actualterm()[i]);
    }
    return g;
}


template <class T>
n_polynom<T> To_n_polynom(const pol<T> &f,int i)
{
    n_polynom<T> g;
    n_expo X(0,i+1);

    for(const auto &pf : f) {
        X[i]=pf.actualdegree();
        g.insert(pf.actualcoef(),X);
        X[i]=0;
    }
    return g;
}

template <class T>
n_polynom<T> To_n_polynom(pol<T> &&f,int i)
{
    n_polynom<T> g;
    n_expo X(0,i+1);

    for(const auto &pf : f) {
        X[i]=pf.actualdegree();
        g.insert(std::move(pf.actualcoef()),X);
        X[i]=0;
    }
    return g;
}


#endif // POL_H

#ifdef S_POLYNOM_H_INCLUDED

template <class T>
s_polynom<T> To_s_polynom(const n_polynom<T> &f)
{
    val::s_polynom<T> g;

    if (f.iszero()) return g;

    int ord=n_expo::getordtype(),n=0,i;
    //auto  pf = f.begin();

    // Bestimme s_expo::n;
    for (const auto &pf : f) {      //for (;pf;pf++) {
        n=val::Max(n,pf.actualterm().dimension());
    }
    n=val::Max(n,n_polynom<T>::getstaticexpodim());
    n = val::Max(n,s_expo::getdim());
    val::s_expo::setdim(n);
    val::s_expo::setordmatrix(val::n_expo::getordmatrix());
    val::s_expo::setordtype(ord);
    val::s_expo X;

    for (const auto &pf : f) {  //for (pf=f;pf;pf++) {
        for (i=0;i<n;i++) X[i]=pf.actualterm()[i];
        g.insert(pf.actualcoef(),X);
    }
    return g;
}


template <class T>
n_polynom<T> To_n_polynom(const s_polynom<T> &f)
{
    val::n_polynom<T> g;

    if (f.iszero()) return g;

    int ord=val::s_expo::getordtype(),n=s_expo::getdim(),i;
    //auto  pf = f.begin();

    val::n_expo::setordtype(ord);
    val::n_expo::setordmatrix(s_expo::getordmatrix());
    val::n_expo X(n);

    for (const auto& pf : f) {    //for (;pf;pf++) {
        for (i=0;i<n;i++) X[i]=pf.actualterm()[i];
        g.insert(pf.actualcoef(),X);
    }
    return g;
}

#endif // S_POLYNOM_H_INCLUDED

#endif // N_POLYNOM_H_INCLUDED


#ifdef S_POLYNOM_H_INCLUDED

template <class T>
int isunivariate(const val::s_polynom<T>& f)
{
    int i=-1,j,n=val::s_expo::getdim();
    //val::s_polynomIterator<T> p;

    if (f.iszero()) return 0;
    for (const auto& p : f) {    //for (p=f;p;p++) {
        for (j=0;j<n;j++) {
            if (p.actualterm()[j]!=0) {
                if (i==-1) i=j;
                else if (i!=j) return -1;
            }
        }
    }
    return i;
}


#ifdef POL_H
template <class T>
s_polynom<T> To_s_polynom(const pol<T> &f,int i)
{
    s_polynom<T> g;

    if (s_expo::getdim()<i+1) s_expo::setdim(i+1);
    s_expo X(0);

    for(const auto &pf : f) {
        X[i]=pf.actualdegree();
        g.insert(pf.actualcoef(),X);
        X[i]=0;
    }
    return g;
}

template <class T>
s_polynom<T> To_s_polynom(pol<T> &&f,int i)
{
    s_polynom<T> g;

    if (s_expo::getdim()<i+1) s_expo::setdim(i+1);
    s_expo X(0);

    for(const auto &pf : f) {
        X[i]=pf.actualdegree();
        g.insert(std::move(pf.actualcoef()),X);
        X[i]=0;
    }
    return g;
}


template <class T>
pol<T> To_unipol(const s_polynom<T> &f,int i)
{
    pol<T> g;

    for (const auto &pf : f) {
        g.insert(pf.actualcoef(),pf.actualterm()[i]);
    }
    return g;
}

#endif // POL_H

#endif // S_POLYNOM_H_INCLUDED


#ifdef RATION_H

template <template <typename> class poly>
val::rational content(const poly<val::rational> &g)
{
    if (g.iszero()) return val::rational(1);

    val::integer z,n,eins(1),minuseins(-1);
    auto pg = g.begin();
    z = val::nominator(g.LC());
    n = val::denominator(g.LC());
    for (pg++;pg;pg++) {
        if (z!=eins && z!=minuseins) z = val::ggTspez(z,val::nominator(pg.actualcoef()));
        n = val::lcm(n,val::denominator(pg.actualcoef()));
    }
    return val::rational(std::move(z),std::move(n));
}


template <template <typename> class poly>
void primitivpart(const poly<val::rational> &g,poly<val::integer> &f,val::rational &cont)
{
    if (g.iszero()) {
        cont = val::rational(1);
        f.del();
        return;
    }
    f.del();
    cont = content(g);
    for (const auto &pg : g) f.insert(pg.actualcoef()/cont,pg.actualterm());
}


template <template <typename> class poly>
poly<val::integer> primitivpart(const poly<val::rational> &g)
{
    val::rational cont;
    poly<val::integer> f;
    primitivpart(g,f,cont);
    return f;
}

template <template <typename> class poly>
int isintegerpol(const poly<val::rational> &g)
{
    integer one(1);
    for (const auto &pg : g) if (denominator(pg.actualcoef())!=one) return 0;
    return 1;
}


template <template <typename> class poly>
poly<val::rational> toRationalPolynom(const poly<val::integer> &g)
{
    poly<val::rational> f;
    for (const auto& pg : g) f.insert(pg.actualcoef(),pg.actualterm());
    return f;
}



template <template <typename> class poly>
std::string MPolToString(const poly<rational> &f)
{
    int i,j,nenner,deg,dim;
    auto it=f.begin();
    std::string s="";
    val::rational value,one(1),minusone(-1),zero;

    if (f.iszero()) return "0";
    for (i=0;it;it++,++i) {
        nenner=0;
        dim=it.actualterm().dimension();
        deg=0;
        for (j=0;j<dim;++j) deg+=it.actualterm()[j];
        if (i==0 && deg==0) {
            s+=val::ToString(it.actualcoef());
            return s;
        }
        if (it.actualcoef()>zero) {
            if (i!=0) s+='+';
        }
        else s+='-';
        if (val::denominator(it.actualcoef())!=val::integer(1) && dim!=0) {
            nenner=1;
            s+='(';
        }
        value=val::abs(it.actualcoef());
        if (value!=one || deg==0) {
            s+=val::ToString(value);
        }
        if (nenner) s+=')';

        if (deg) {
            if (it.actualcoef()!=one && it.actualcoef()!=minusone) s+='*';
            for (j=0;j<dim;++j) {
                if (it.actualterm()[j]) break;
            }
            // First variable:
            deg= it.actualterm()[j];
            s+="x" + val::ToString(j+1);
            if (deg!=1) s+="^" + val::ToString(deg);
            for(++j;j<dim;++j) {
                deg= it.actualterm()[j];
                if (deg) {
                    s+="*x" + val::ToString(j+1);
                    if (deg!=1) s+="^" + val::ToString(deg);
                }
            }
        }
    }
    return s;
}


template <template <typename> class poly>
std::string PolToString(const poly<rational> &f)
{
    int i,nenner;
    auto it=f.begin();
    std::string s="";
    val::rational value,one(1),minusone(-1),zero;

    if (it == 0) return "0";

    for (i=0;it;it++,++i) {
        nenner=0;
        if (i==0 && it.actualdegree()==0) {
            s+=val::ToString(it.actualcoef());
            return s;
        }
        if (it.actualcoef()>zero) {
            if (i!=0) s+='+';
        }
        else s+='-';
        if (val::denominator(it.actualcoef())!=val::integer(1) && it.actualdegree()!=0) {
            nenner=1;
            s+='(';
        }
        value=val::abs(it.actualcoef());
        if (value!=one || it.actualdegree()==0) {
            s+=val::ToString(value);
        }
        if (nenner) s+=')';

        if (it.actualdegree()) {
            if (it.actualcoef()!=one && it.actualcoef()!=minusone) s+='*';
            s+='x';
            if (it.actualdegree()!=1) s+="^" + val::ToString(it.actualdegree());
        }
    }
    return s;
}


#endif // RATION_H

#ifdef INTEGER_H
template <template <typename> class poly>
val::integer content(const poly<val::integer> &g)
{
    if (g.iszero()) return val::integer(1);

    val::integer z,eins(1),minuseins(-1);
    auto pg = g.begin();

    if (g.LC() == eins || g.LC() == minuseins) return val::integer(eins);
    else z = g.LC();

    for (pg++;pg != g.end();pg++) {
        z = val::ggTspez(z,(*pg).actualcoef());
        if (z == eins || z== minuseins) return z;
    }
    return z;
}



template <template <typename> class poly>
void primitivpart(const poly<val::integer> &g,poly<val::integer> &f,val::integer &cont)
{
    cont = content(g);
    if (cont==val::integer(1) || cont==val::integer(-1)) {
        f=g;
        return;
    }
    f.del();
    for (const auto &pf : g) {
        f.insert(val::EDIV(pf.actualcoef(),cont),pf.actualterm());
    }
}


template <template <typename> class poly>
poly<val::integer> primitivpart(const poly<val::integer> &g)
{
    val::integer cont;
    poly<val::integer> f;
    primitivpart(g,f,cont);
    return f;
}


template <typename T,template <typename> class poly>
poly<integer> modPolToIntPol(const poly<T> &f)
{
    poly<integer> g;
    for (const auto& pf : f) {
        g.insert(integer(pf.actualcoef()),pf.actualterm());
    }
    return g;
}


#endif // INTEGER_H


#ifdef MODQ_H

template <template <typename> class poly>
poly<val::modq> toModqPolynom(const poly<val::integer> &g)
{
    val::integer p(modq::q);
    poly<val::modq> f;

    for (const auto& pg : g) f.insert(val::modq(pg.actualcoef()%p),pg.actualterm());
    return f;
}

#endif // MODQ_H


#ifdef MODINT_H_INCLUDED
template <template <typename> class poly>
poly<val::modint> toModintPolynom(const poly<val::integer> &f,int q)
{
    integer iq(q);
    poly<val::modint> g;

    for ( const auto& t : f) {
        g.insert(val::modint(int(t.actualcoef()%iq),q),t.actualterm());
    }
    return g;
}



#endif // MODINT_H_INCLUDED


template <typename T,template <typename> class poly>
poly<double> ToDoublePolynom(const poly<T> &g,const double& eps)
{
    poly<double> f;
    double v;

    for (const auto & pg : g) {
        v = double(pg.actualcoef());
        if (val::abs(v)>=eps) f.insert(v,pg.actualterm());
    }
    return f;
}




template <typename T,template <typename> class poly>
std::string PolToString(const poly<T> &f)
{
    int i;
    auto it=f.begin();
    std::string s="";
    T one=val::unity_element<T>(),minusone=-one,zero=val::zero_element<T>();

    if (it == 0) return "0";
    if (f.degree()==0) return val::ToString(it.actualcoef());

    for (i=0;it;it++,++i) {
        if (i!=0 && it.actualcoef()>zero) s+='+';
        if (it.actualcoef()==minusone) {
            s+='-';
            if (it.actualdegree()==0) s+='1';
        }
        else if (it.actualcoef()!=one || it.actualdegree()==0) s+=val::ToString(it.actualcoef());
        if (it.actualdegree()) {
            if (val::abs(it.actualcoef())!=one) s+='*';
            s+='x';
            if (it.actualdegree()!=1) s+="^" + val::ToString(it.actualdegree());
        }
    }
    return s;
}


template <template <typename> class poly>
std::string PolToString(const poly<modq> &f)
{
    int i;
    auto it=f.begin();
    std::string s="";
    modq one(1),zero;

    if (it == 0) return "0";
    if (f.degree()==0) return val::ToString(it.actualcoef());

    for (i=0;it;it++,++i) {
        if (i!=0 && it.actualcoef() != zero) s+='+';
        if (it.actualcoef()!=one || it.actualdegree()==0) s+=val::ToString(it.actualcoef());
        if (it.actualdegree()) {
            if (it.actualcoef()!=one) s+='*';
            s+='x';
            if (it.actualdegree()!=1) s+="^" + val::ToString(it.actualdegree());
        }
    }
    return s;
}




} // end namespace val

#endif // POL_ARITHMETIC_H_INCLUDED
