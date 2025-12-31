
#include <s_groebner.h>
#include "s_groebner/common_bb.cpp"
#include "s_groebner/bb_modq.cpp"
#include "s_groebner/bb_int.cpp"
#include "s_groebner/bbhom.cpp"
#include "s_groebner/hilbert.cpp"
#include "s_groebner/minpol.cpp"
#include "s_groebner/walk.cpp"



namespace val
{
//
void SetWriteText(strg_function &Func)
{
    common_bb::writetextoutput=&Func;
}

void SetMyMessage(strg_function &Func)
{
    common_bb::messageoutput=&Func;
}

void SetClear(void_function &Func)
{
    common_bb::clearoutput=&Func;
}

void WriteText(const std::string &s)
{
    common_bb::WriteText(s);
}


void SetPrimlistPath(const std::string& s)
{
    common_bb::primlistpath=s;
}


void SetNumberofThreads(int n)
{
    common_bb::SetComputingThreads(n);
}

//

int getFileType(const std::string &name)
{
    return common_bb::gettype(name);
}


template <>
void saveGBasis(const std::string &arg,const val::Glist<val::s_polynom<modq> >& P)
{
    common_bb::savelist(arg,P);
}


template <>
void saveGBasis(const std::string& arg,const val::Glist<val::s_polynom<integer> >& P)
{
    common_bb::savelist(arg,P);
}


template<>
int readPolynomBasis(const std::string& arg,Glist<s_polynom<integer> > & G,int onlytotdeghom)
{
    if (!G.isempty()) G.dellist();
    return common_bb::readfromfile(arg,G,onlytotdeghom);
}


template<>
int readPolynomBasis(const std::string& arg,Glist<s_polynom<modq> > & G,int onlytotdeghom)
{
    if (!G.isempty()) G.dellist();
    return common_bb::readfromfile(arg,G,onlytotdeghom);
}



int groebner(const std::string &name,Glist<s_polynom<modq> > &G)
{
    if (!G.isempty()) G.dellist();
    return bb_modq::bbgccmain(name,G);
}

int groebner(const std::string &name,Glist<s_polynom<integer> > &G)
{
    if (!G.isempty()) G.dellist();
    return bb_int::bbhommodmain(name,G);
}


int groebner(val::Glist<val::s_polynom<val::integer>> &G,int comment)
{
    return bb_int::bbhommod(G,comment);
}

int groebner(val::Glist<val::s_polynom<val::modq>> &G,int comment)
{
    G.sort();
    return bb_modq::Groebner(G,comment);
}

int primitiv_groebner(val::Glist<val::s_polynom<val::integer>> &H)
{
    if (H.isempty()) return 0;
    int nG;
    val::Glist<val::s_polynom<val::integer>> G;
    nG = common_bb::primGroebner(G,H);
    H = std::move(G);
    return nG;
}


int primitiv_groebner(val::Glist<val::s_polynom<val::modq>> &H)
{
    if (H.isempty()) return 0;
    int nG;
    val::Glist<val::s_polynom<val::modq>> G;
    nG = common_bb::primGroebner(G,H);
    H = std::move(G);
    return nG;
}


template<>
int isgroebner(const std::string &name,Glist<s_polynom<modq> > &G,int anzeige)
{
    ChronoClass Chrono;
    int m,is;
    std::string s;

    if (!G.isempty()) G.dellist();
    m=common_bb::readfromfile(name,G);

    Chrono();
    is=common_bb::isgroebner(G,m,anzeige);


    s+="\n\nTime: " + ToString(Chrono());
    if (is) s+= "\nG is groebner basis!";
    else s+="\nG is not a groebner basis!";
    s+="\nReductions: " + ToString(s_polynom<modq>::nreduction);
    common_bb::Clear();
    common_bb::WriteText(s);
    s_polynom<modq>::nreduction=0;
    return is;
}


template<>
int isgroebner(const std::string &name,Glist<s_polynom<integer> > &G,int anzeige)
{
    ChronoClass Chrono;
    int m,is;
    std::string s;

    if (!G.isempty()) G.dellist();
    m=common_bb::readfromfile(name,G);

    Chrono();
    is=common_bb::isgroebner(G,m,anzeige);

    s+="\nTime: " + ToString(Chrono());
    if (is) s+= "\nG is groebner basis!";
    else s+="\nG is not a groebner basis!";
    s+="\nMaximal integer-length: " + ToString(integer::GetMaxlength());
    s+="\nReductions: " + ToString(s_polynom<integer>::nreduction);
    common_bb::Clear();
    common_bb::WriteText(s);
    s_polynom<integer>::nreduction=0;
    return is;
}


template <>
int homgroebner(const std::string &name,Glist<s_polynom<integer> > &G)
{
    return bbhom::bbhommain(name,G);
}

template <>
int homgroebner(const std::string &name,Glist<s_polynom<modq> > &G)
{
    return bbhom::bbhommain(name,G);
}

template <>
int homgroebner(Glist<s_polynom<integer> > &G)
{
    return bbhom::bbhommain(G);
}

template <>
int homgroebner(Glist<s_polynom<modq> > &G)
{
    return bbhom::bbhommain(G);
}


template <>
int hilbertconversion(const std::string &name,Glist<s_polynom<integer> > &G,int order,const matrix<int> &M)
{
    return bbhom::hilbertconversionmain(name,G,order,M);
}

template <>
int hilbertconversion(const std::string &name,Glist<s_polynom<modq> > &G,int order,const matrix<int> &M)
{
    return bbhom::hilbertconversionmain(name,G,order,M);
}

template <>
int hilbertconversion(Glist<s_polynom<integer> > &G,int order,const matrix<int> &M)
{
    return bbhom::hilbertconversionmain(G,order,M);
}

template <>
int hilbertconversion(Glist<s_polynom<modq> > &G,int order,const matrix<int> &M)
{
    return bbhom::hilbertconversionmain(G,order,M);
}




pol<rational> Hilbertpolynomial(const std::string &name,int affin)
{
    return hilbert::Hilbertpolynomial(name,affin);
}



pol<rational> Hilbertpolynomial(const Glist<s_polynom<integer>>& G,int affin)
{
    val::vector<s_expo> In(G.length());
    int i=0,ord=s_expo::getordtype();

    for (const auto& f : G) {
        In(i) = f.LT();
        ++i;
    }
    s_expo::setordtype(-1);
    In.sort();
    s_expo::setordtype(ord);
    return hilbert::Hilbertpolynomial(In,affin);
}

pol<rational> Hilbertpolynomial(const Glist<s_polynom<modq>>& G,int affin)
{
    val::vector<s_expo> In(G.length());
    int i=0,ord=s_expo::getordtype();

    for (const auto& f : G) {
        In(i) = f.LT();
        ++i;
    }
    s_expo::setordtype(-1);
    In.sort();
    s_expo::setordtype(ord);
    return hilbert::Hilbertpolynomial(In,affin);
}



pol<rational> minimalpolynom(const Glist<s_polynom<integer> > &G,int k,int comment)
{
    return minpol::minimalpolynom(G,k,comment);
}

pol<modq> minimalpolynom(const Glist<s_polynom<modq> > &G,int k,int comment)
{
    return minpol::minimalpolynom(G,k,comment);
}



template <>
void groebnerwalk(val::Glist<val::s_polynom<modq> > &G,int neworder,int k,int l,const val::matrix<int> &M,int comment)
{
    walk::walkmain(G,neworder,k,l,M,comment);
}

template <>
void groebnerwalk(val::Glist<val::s_polynom<integer> > &G,int neworder,int k,int l,const val::matrix<int> &M,int comment)
{
    walk::walkmain(G,neworder,k,l,M,comment);
}

template <>
void groebnerwalk(const std::string& name,val::Glist<val::s_polynom<modq> > &G,int neworder,int k,int l,const val::matrix<int> &M)
{
    walk::walkmain(name,G,neworder,k,l,M);
}



template <>
void groebnerwalk(const std::string& name,val::Glist<val::s_polynom<integer> > &G,int neworder,int k,int l,const val::matrix<int> &M)
{
    walk::walkmain(name,G,neworder,k,l,M);
}


int isinG(const std::string& file,const std::string& gbfile)
{
    return common_bb::isinG(file,gbfile);
}


int Iddim(const std::string& file)
{
    return hilbert::Hilbertpolynomial(file).degree();
}


template <>
int reduceGroebner(const std::string& name,Glist<s_polynom<modq> >& G)
{
    int nG;
    std::string s;
    ChronoClass Chrono;

    s_polynom<modq>::nreduction=0;
    if (!G.isempty()) G.dellist();
    common_bb::readfromfile(name,G);
    Chrono();
    nG=common_bb::minimalGroebner(G);
    common_bb::interredBasis(G);
    s+="\nTime: " + ToString(Chrono());
    s+="\nElements in G: " + ToString(nG);
    s+="\nMonomials: " + ToString(s_polynom<modq>::getmnumber());
    s+="\nReductions: "+ ToString(s_polynom<modq>::nreduction);
    common_bb::WriteText(s);
    return nG;
}

template <>
int reduceGroebner(const std::string& name,Glist<s_polynom<integer> >& G)
{
    int nG;
    std::string s;
    ChronoClass Chrono;

    s_polynom<integer>::nreduction=0;
    if (!G.isempty()) G.dellist();
    common_bb::readfromfile(name,G);
    Chrono();
    nG=common_bb::minimalGroebner(G);
    common_bb::interredBasis(G);
    s+="\nTime: " + ToString(Chrono());
    s+="\nElements in G: " + ToString(nG);
    s+="\nMonomials: " + ToString(s_polynom<integer>::getmnumber());
    s+="\nReductions: "+ ToString(s_polynom<integer>::nreduction);
    common_bb::WriteText(s);
    return nG;
}


template <>
int reduceGroebner(Glist<s_polynom<modq> >& G, int comment)
{
    int nG;
    std::string s;
    ChronoClass Chrono;

    s_polynom<modq>::nreduction=0;
    Chrono();
    nG=common_bb::minimalGroebner(G);
    common_bb::interredBasis(G);
	if (comment) {
		s+="\nTime: " + ToString(Chrono());
		s+="\nElements in G: " + ToString(nG);
		s+="\nMonomials: " + ToString(s_polynom<modq>::getmnumber());
		s+="\nReductions: "+ ToString(s_polynom<modq>::nreduction);
		common_bb::WriteText(s);
	}
    return nG;
}

template <>
int reduceGroebner(Glist<s_polynom<integer> >& G, int comment)
{
    int nG;
    std::string s;
    ChronoClass Chrono;

    s_polynom<integer>::nreduction=0;
    Chrono();
    nG=common_bb::minimalGroebner(G);
    common_bb::interredBasis(G);
	if (comment) {
		s+="\nTime: " + ToString(Chrono());
		s+="\nElements in G: " + ToString(nG);
		s+="\nMonomials: " + ToString(s_polynom<integer>::getmnumber());
		s+="\nReductions: "+ ToString(s_polynom<integer>::nreduction);
		common_bb::WriteText(s);
	}
    return nG;
}

int homogenize(const std::string& i_name,const std::string& o_name,int comment)
{
    return common_bb::homogenize(i_name,o_name,comment);
}

int dehomogenize(const std::string& i_name,const std::string& o_name,int comment)
{
    return common_bb::dehomogenize(i_name,o_name,comment);
}

int preparenormalpos(const std::string& i_name,const std::string& o_name,int comment)
{
    return common_bb::preparenormalpos(i_name,o_name,comment);
}

} // end namespace val
