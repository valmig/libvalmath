
#include <ideal_roots.h>
#include <s_groebner.h>
#include <s_polynom.h>
#include <rational.h>
#include <pol.h>
#include <analysis.h>
#include <fstream>
#include <vector.h>
#include <complex.h>
#include <val_utils.h>
#include <d_array.h>
#include <matrix.h>
#include <pol_arithmetic.h>
#include "modint_minpol.cpp"
#include <numbers.h>
#include <rand.h>
#include <thread>


/*
#ifdef _WIN32
std::string primlispath = "C:\\includegcc\\Mathe\\primlis1.txt";
#else
std::string primlispath = "/home/miguel/includegcc/mathe/primlis1.txt";
#endif // _WIN32
*/

std::string primlispath = std::string(val::primlistpath);

namespace id_roots
{
val::Glist<int> Prime;


val::matrix<int> eliminationmatrix(const val::vector<int> &selection)
{
    int dim=selection.dimension(),i,j,zeile=0,k,l=0;
    val::matrix<int> M(0,dim,dim);

    for (i=0;i<dim;i++) if (selection[i]) l++;

    if (l==dim) {
        for (i=0;i<dim;i++) M(0,i)=1;
        for (i=1;i<dim;i++) M(i,dim-i) =-1;
        return M;
    }

    k=dim -l;

    // Create Matrix

    for (i=0;i<dim;i++) {
        if (!selection[i]) M(zeile,i)=1;
    }

    i=dim-1;
    for (j=k-1;j>0;j--) {
         zeile++;
         for (;i>=0;i--) {
             if (!selection[i]) {M(zeile,i)=-1;i--;break;}
        }
    }
    zeile++;
    for (i=0;i<dim;i++) {
        if (selection[i]) M(zeile,i)=1;
    }

    i=dim-1;
    for (j=l-1;j>0;j--) {
        zeile++;
        for (;i>=0;i--) {
            if (selection[i]) {M(zeile,i) =-1;i--;break;}
        }
    }

    return M;
}


val::vector<val::pol<double>> getdoublepolynomials(const val::Glist<val::s_polynom<val::integer>> &G,int index)
{
    using namespace val;
    int i,j,dim=G.length();
    val::vector<val::pol<double>> F(dim);
    val::s_polynom<rational> f_r;
    val::s_polynomIterator<rational> it;

    f_r = val::toRationalPolynom(G[0]);
    f_r.normalize();
    F(index) = val::ToDoublePolynom(val::To_unipol(f_r,index));

    for (i=1;i<dim;++i) {
        f_r = val::toRationalPolynom(G[i]);
        f_r.normalize();
        it = f_r.begin();
        for (j=0;j<dim;++j) {
            if ((it.actualterm())[j]) break;
        }
        for (++it;it;++it) {
            F(j).insert(double(it.actualcoef()),(it.actualterm())[index]);
        }
    }

    return F;
}


void extendordermatrix()
{
    using namespace val;
    int i,order=s_expo::getordtype(),n=s_expo::getdim();
	if (order==-1000) {
		val::matrix<int> Mnew(0,n+1,n+1);
		int j;
		Mnew(0,0)=1;
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) Mnew(i+1,j+1) = (s_expo::getordmatrix())(i,j);
        }

		val::s_expo::setordtype(-1000);
		val::s_expo::setordmatrix(Mnew);
	}
	else if (order==-1) val::s_expo::setordtype(-1);
	else if (order==0) {
        val::matrix<int> Mnew(0,n+1,n+1);
        Mnew(0,0)=1;
        for (i=0;i<n;i++) Mnew(1,i+1)=1;
        for (i=2;i<=n;i++) Mnew(i,i-2)=1;
		val::s_expo::setordtype(-1000);
		val::s_expo::setordmatrix(Mnew);
	}
	else { // order=-2
        val::matrix<int> Mnew(0,n+1,n+1);
        Mnew(0,0)=1;
        for (i=0;i<n;i++) Mnew(1,i+1)=1;
        for (i=2;i<=n;i++) Mnew(i,n+2-i) =-1;
		val::s_expo::setordtype(-1000);
		val::s_expo::setordmatrix(Mnew);
	}
    val::s_expo::setdim(n+1);
}

val::Glist<val::s_polynom<val::integer>> preparenormalposition(const val::Glist<val::s_polynom<val::integer>> &G)
{
    using namespace val;
    int i,n = s_expo::getdim();
    s_polynom<integer> f;
    Glist<val::s_polynom<integer>> H;

	val::s_expo X(0);
	val::rational coeff,zero;
	val::RandomClass Rand;

	// First polynomial : f = X0 + ciXi
	X[0] = 1;
	f.insert(integer(1),X);
	for (i=1;i<n;++i) {
        X = s_expo(0);
        X[i] = 1;
        f.insert(integer(Rand.random(-50,50)),X);
	}
	H.sinsert(std::move(f));

	for (const auto& g : G) {
        for (const auto& it : g) {
            X = s_expo(0);
            for (i=0;i<n-1;++i) X[i+1] = it.actualterm()[i];
            f.insert(it.actualcoef(),X);
        }
        H.sinsert(std::move(f));
	}

	return H;
}


} // end namespace id_roots



val::Glist<int> CreatePrimlist(const val::Glist<val::s_polynom<val::integer> > &Gint)
{
    using namespace val;
    Glist<int> Primelist;
    int lucky;
    integer zero;

    if (id_roots::Prime.isempty()) id_roots::Prime = val::getPrimelist();


    for (const auto & p : id_roots::Prime) {
        lucky=1;
        for (const auto& it : Gint) {
            if ((it.LC())%integer(p) == zero) lucky=0;
        }
        if (lucky) Primelist.push_back(p);
    }
    return Primelist;
}



void CreatemodintBasis(const val::Glist<val::s_polynom<val::integer> > &Gint,val::Glist<val::s_polynom<val::modint> > &Gmodq,int q)
{
    using namespace val;
    Gmodq.dellist();
    GlistIterator<s_polynom<integer> > itG;
    s_polynomIterator<integer> itf;
    s_polynom<modint> f;

    for (itG=Gint;itG;itG++) {
        for (itf=itG.actualvalue();itf;itf++) {
            f.insert(modint(int(itf.actualcoef()%integer(q)),q),itf.actualterm());
        }
        f.normalize();
        if (!f.iszero()) Gmodq.inserttoend(std::move(f));
    }
}


void ComputeFint(const val::Glist<val::s_polynom<val::integer> > &Gint,val::pol<val::integer>& Fint,int k,int q)
{
	using namespace val;
    Glist<s_polynom<modint> > Gmodq;
    pol<modint> fmodq;

    CreatemodintBasis(Gint,Gmodq,q);
    fmodq=minpol::minimalpolynom(Gmodq,k,q,0);
    Fint=val::modPolToIntPol(fmodq);
}
// ---------------------------------------------------------------------------------------------------------------------



// -------------------------------   squarefree part  --------------------------------------------------------------------------

val::pol<val::modq> p_root(const val::pol<val::modq> &f)
{
    using namespace val;
    if (f.degree()<=1) return f;

    int p=modq::q;
    pol<modq> g;
    polIterator<modq> pf;

    for (pf=f;pf;pf++) g.insert(pf.actualcoef(),pf.actualdegree()/p);

    return g;
}

val::pol<val::modq> squarefreepart(const val::pol<val::modq>& f0)
{
    using namespace val;
    if (f0.iszero() || f0.degree()<=1) return f0;
    pol<modq> f=f0,h,h0,h1,g;

    h=gcd(f,f.derive());
    while (h.degree()!=0) {
        if ((h.derive()).iszero()) {
            g=p_root(h);
            f*=g;
            f/=h;
        }
        else {
            h0=h;
            do {
                h1=h.derive();
                h=gcd(h,h1);
            }
            while(!h1.iszero());
            g=p_root(h);
            f*=g;
            f/=h0;
        }
        h=gcd(f,f.derive());
    }

    return f;
}


val::pol<val::rational> squarefreepart(const val::pol<val::rational>& f0,const val::Glist<int> &Prime)
{
     if (f0.iszero() || f0.degree()<=1) return f0;
     val::pol<val::rational> h1 = val::modular_gcd(f0,f0.derive(),Prime);

     return (f0/h1);
}


void computesquarefreepartmodq(const val::d_array<val::pol<val::modq>>& f,val::d_array<val::s_polynom<val::modq>>& F,int i1,int i2,int comment=0)
{
    val::pol<val::modq> sqr_f;
    for (int i=i1;i<i2;++i) {
        sqr_f = squarefreepart(f[i]);
        if (comment) val::WriteText("\n compute square-free part for i = " + val::ToString(i) + " ; degree = " + val::ToString(sqr_f.degree()));
        if (f[i].degree() != sqr_f.degree()) {
			F[i] = val::To_s_polynom(sqr_f,i);
		}
    }
}


void computeminimalpolynomial(const val::Glist<val::s_polynom<val::modq>> &G,val::d_array<val::pol<val::modq>>& f,int i1,int i2,int comment=0)
{
    for (int i=i1;i<i2;++i) {
        if (comment) val::WriteText("\n compute minimal polynomial for i = " + val::ToString(i));
        f[i] = val::minimalpolynom(G,i,0);
    }
}

// --------------------------------------------------------------------------------------------------------------------------------------------

namespace id_roots
{

int size_of_coeff(const val::pol<val::rational> &f)
{
    int sc=0,msc;
    for (const auto & m : f) {
        msc = val::Max(m().nominator().abslength(),m().denominator().abslength());
        sc = val::Max(sc,msc);
    }
    return sc;
}

void Message_to_cout(const std::string& s)
{
    std::cout<<s;
}

void (*messageoutput) (const std::string&) = &Message_to_cout;

void MyMessage(const std::string& s)
{
    messageoutput(s);
}


const val::pol<double>& RoundPol(val::pol<double>& f,const double &eps)
{
 using namespace val;
 if (f.iszero()) return f;

 int i=f.degree();

 for (;i>=0;i--) if (val::abs(f[i])>=eps) break;
 if (i==-1) {
    f=pol<double>();
    return f;
 }
 pol<double> g(f[i],i);
 i--;
 for (;i>=0;i--)
     if (val::abs(f[i])>=eps) g+=pol<double>(f[i],i);
 f=std::move(g);
 return f;
}



int isinnormalposition(const std::string &name,val::vector<val::pol<double> > &F,int &var_index)
{

	std::ifstream file(name,std::ios::in); 
	if (!file) {                                     
		MyMessage("\nCannot read file " + name + " !\n");
        return 0;
	}

	int i,j,dim,ordtype,n=0,is,uni_index=-1;
	val::matrix<int> M;
    val::s_polynom<val::rational> f;
    val::s_polynomIterator<val::rational> It_f;

	file>>dim>>ordtype;

	val::s_expo::setordtype(ordtype);
	val::s_expo::setdim(dim);
	F=val::vector<val::pol<double> >(dim);

	if (ordtype==-1000) {
		M=val::matrix<int>(dim);
		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++) file>>M(i,j);
		val::s_expo::setordmatrix(std::move(M));
	}
	else if (ordtype!=-2 && ordtype!=-1 && ordtype!=0) val::s_expo::setordtype(-2);

	val::rational coef,zero;
	val::s_expo X;
	val::vector<val::s_polynom<val::rational> > G(dim);

	do {
		file>>coef;
		if (coef!=zero) {
            file>>X;
			f.insert(std::move(coef),X);
        }
		else break;
		do {
			file>>coef;
			if (coef!=zero) {
			    file>>X;
			    f.insert(std::move(coef),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			if (n==dim) return 0;
			G(n)=std::move(f);
			n++;
		}
	}
	while (file);
	file.close();

	if (n!=dim) return 0;

	val::s_expo degX;

	// check if zero-dimensional:

	if (!val::iszerodimensional(G,degX)) return 0;


	// Search univariate polynomial:
	is=-1;
	for (i=0;i<dim;++i) {
        if ((is=val::isunivariate(G(i))) != -1) {
            uni_index=i;
            var_index=is;
            break;
        }
	}
	if (is==-1) return 0;

	// Check if all polynomials but G(uni_index) are of degree 1;
	for (i=0;i<n;i++) {
        if (i==var_index) continue;
        if (degX[i]!=1) return 0;
	}

	// Check if all polynomials, but G(uni_index), are univariate in all terms but their leading term:
	for (i=0;i<n;i++) {
        if (i==uni_index) continue;
        It_f=G(i);It_f++;
        for(;It_f;It_f++) {
            for(j=0;j<dim;j++) {
                if (j==var_index) continue;
                if (It_f.actualterm()[j]!=0) return 0;
            }
        }
	}
	// Basis is in normal position: Write polynomials as pol<double>;
	for (i=0;i<n;i++) G(i).normalize();
	for (i=0;i<n;i++) {
        It_f=G(i);
        if (i!=uni_index) {
                for (j=0;j<dim;j++)
                    if (G(i).LT()[j]!=0) break;
                It_f++;
        }
        else j=var_index;
        for (;It_f;It_f++) F[j].insert(double(It_f.actualcoef()),It_f.actualterm()[var_index]);
	}

    return 1;
}

} // end namespace id_roots




namespace val
{

// modint Minimalpolynom:    --------------------------------------------------------------------------

val::pol<val::rational> modint_minimalpolynom(const val::Glist<val::s_polynom<val::integer> > &Gint,int k,
                                              val::Glist<int>& Primlist,int n,int comment)
{
    int i=0,j,isempty=0;
    pol<rational> fold,f;
    pol<integer> g;
    d_array<pol<integer> > Fint(n+1);
    d_array<integer> M(n+1),m(n+1);
    d_array<int> q(n);
    d_array<std::thread*> t(n);

    m[0]=M[0] = integer(1);

    if (Primlist.isempty()) {
        Primlist = CreatePrimlist(Gint);
    }

    auto p = Primlist.begin();

    if (Gint.isempty()) return f;

    for (i=0;!isempty;i++) {
        for (j=0;j<n;j++) {
            if (p==Primlist.end()) {isempty=1;break;}
            q[j]=p();
            m[j+1] = integer(q[j]);
            M[j+1] = M[j] * integer(q[j]);
            ++p;
		}
        if (isempty) {f.del();break;}
        for (j=0;j<n;++j) t[j] = new std::thread(ComputeFint,std::cref(Gint),std::ref(Fint[j+1]),k,q[j]);
        for (j=0;j<n;++j) t[j]->join();
        for (j=0;j<n;++j) delete t[j];

        g=simultanpolynomial(Fint,m,M);
        M[0] = std::move(M[n]);
        m[0] = M[0];
        Fint[0] = std::move(g);
        f=rationalconstruction(Fint[0],M[0]);
        if (comment) WriteText("\n i = " + ToString(i));
        if (!f.iszero() && f==fold) break;
        if (comment) WriteText("\n f = \n" + ToString(f));
        fold=std::move(f);
    }
    return f;
}


//  ------------------------------  radical-ideal:

int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::integer>>& G,int nthreads,int comment)
{
    if (G.isempty()) return -1;

    int  i,n = val::s_expo::getdim(),isradical=1,index=0,maxdeg=0,maxsize=1000000,f_size;
    val::pol<val::rational> sqr_f;
    val::d_array<val::pol<val::rational>> f(n);
    val::d_array<val::pol<val::integer>> F(n);
    val::Glist<int> Primlist;
    val::d_array<std::thread*> thr(nthreads);

    id_roots::Prime = getPrimelist();
    Primlist = CreatePrimlist(G);

    for (i=0;i<n;++i) {
        if (comment) WriteText("\n compute minimalpolynomial for i = " + ToString(i));
        f[i] = modint_minimalpolynom(G,i,Primlist,nthreads,0);
    }

    for (i=0;i<n;++i) {
        if (comment) WriteText("\n compute square-free part for i = " + ToString(i));
        sqr_f = squarefreepart(f[i],id_roots::Prime);
        if (f[i].degree()!=sqr_f.degree()) {
            isradical = 0;
            F[i] = val::primitivpart(sqr_f);
        }
        if (comment) WriteText(" . degree = " + ToString(sqr_f.degree()));
        if (sqr_f.degree()>maxdeg) {
            maxdeg=sqr_f.degree();
            index=i;
        }
        else if (sqr_f.degree()==maxdeg) {
            if ((f_size=id_roots::size_of_coeff(sqr_f))<maxsize) {
                index=i;
                maxsize=f_size;
            }
        }
    }

    if (!isradical) {
        for (i=0;i<n;++i) {
            if (!F[i].iszero()) G.sinsert(val::To_s_polynom(F[i],i));
        }
        if (comment) WriteText("\n <G> is not radical. Compute groebner basis: ");
        val::groebner(G);
        if (comment) WriteText("\n Groebner Basis computed!");
    }
    else if (comment) WriteText("\n <G> is radical.");

    return index;
}


int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::integer>>& G,int &maxdeg,int nthreads,int comment)
{
    if (G.isempty()) return -1;

    int  i,n = val::s_expo::getdim(),isradical=1,index=0,maxsize=1000000,f_size;
    val::pol<val::rational> sqr_f;
    val::d_array<val::pol<val::rational>> f(n);
    val::d_array<val::pol<val::integer>> F(n);
    val::Glist<int> Primlist;
    val::d_array<std::thread*> thr(nthreads);

    id_roots::Prime = getPrimelist();
    Primlist = CreatePrimlist(G);

    maxdeg=0;


    for (i=0;i<n;++i) {
        if (comment) WriteText("\n compute minimalpolynomial for i = " + ToString(i));
        f[i] = modint_minimalpolynom(G,i,Primlist,nthreads,0);
    }

    for (i=0;i<n;++i) {
        if (comment) WriteText("\n compute square-free part for i = " + ToString(i));
        sqr_f = squarefreepart(f[i],id_roots::Prime);
        if (f[i].degree()!=sqr_f.degree()) {
            isradical = 0;
            F[i] = val::primitivpart(sqr_f);
        }
        if (comment) WriteText(" . degree = " + ToString(sqr_f.degree()));
        if (sqr_f.degree()>maxdeg) {
            maxdeg=sqr_f.degree();
            index=i;
        }
        else if (sqr_f.degree()==maxdeg) {
            if ((f_size=id_roots::size_of_coeff(sqr_f))<maxsize) {
                index=i;
                maxsize=f_size;
            }
        }
    }

    if (!isradical) {
        for (i=0;i<n;++i) {
            if (!F[i].iszero()) G.sinsert(val::To_s_polynom(F[i],i));
        }
        if (comment) WriteText("\n <G> is not radical. Compute groebner basis: ");
        val::groebner(G);
        if (comment) WriteText("\n Groebner Basis computed!");
    }
    else if (comment) WriteText("\n <G> is radical.");

    return index;
}


int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::modq>>& G,int nthreads,int comment)
{
    if (G.isempty()) return 1;
    int  i,j,n = val::s_expo::getdim(),isradical=1,d;
    double dd = double(n)/double(nthreads);
    val::d_array<int> i1(nthreads),i2(nthreads);
    val::d_array<val::pol<val::modq>> f(n);
    val::d_array<val::s_polynom<val::modq>> F(n);
    val::d_array<std::thread*> thr(nthreads);

    d = int(val::round(dd,0));

    for (i=0,j=0;i<nthreads;i++,j+=d) {
        i1[i] = j;
        i2[i] = j+d;
    }
    i2[nthreads-1] = n;

    for (i=0;i<nthreads;++i) {
        thr[i] = new std::thread(computeminimalpolynomial,std::cref(G),std::ref(f),i1[i],i2[i],1);
    }
    for (i=0;i<nthreads;++i) {
        thr[i]->join();
        delete thr[i];
    }

    for (i=0;i<nthreads;++i) {
        thr[i] = new std::thread(computesquarefreepartmodq,std::cref(f),std::ref(F),i1[i],i2[i],1);
    }
    for (i=0;i<nthreads;++i) {
        thr[i]->join();
        delete thr[i];
    }

    for (i=0;i<n;++i) {
        if (!F[i].iszero()) {
            isradical = 0;
            break;
        }
    }

    if (!isradical) {
        for (i=0;i<n;++i) {
            if (!F[i].iszero()) G.sinsert(std::move(F[i]));
        }
        if (comment) WriteText("\n <G> is not radical. Compute groebner basis:");
        val::groebner(G);
        if (comment) WriteText("\n Groebner Basis computed!");

    }
    else if (comment) WriteText("\n <G> is radical.");

    return isradical;
}


void SetRootsMessage(strg_function &Func)
{
    id_roots::messageoutput=&Func;
}


int idealroots(const std::string &filename,matrix<double> &MR,matrix<complex> &MC,int nthreads,int comment)
{
    if (val::getFileType(filename)!=2) {
        id_roots::MyMessage("\nWrong type of file!\n");
        return 0;
    }

    int var_index,N,NR,NC,i,j;
    val::vector<val::pol<double> > F;
    val::vector<double> Realzeros;
    val::vector<val::complex> Compzeros;

    if (!id_roots::isinnormalposition(filename,F,var_index)) {
        if (comment) WriteText("\nBasis is not in normal position!\n");
        Glist<s_polynom<integer>> G;
        val::readPolynomBasis(filename,G);
        return computeroots(G,MR,MC,nthreads,comment);
    }

    N=val::roots(F(var_index),Realzeros,Compzeros);
    NR=Realzeros.dimension();
    NC=Compzeros.dimension();

    if (NR>0) MR=matrix<double>(NR,F.dimension());
    else MR=matrix<double>();
    if (NC>0) MC=matrix<complex>(NC,F.dimension());
    else MC=matrix<complex>();

    for (i=0;i<NR;i++) {
        for (j=0;j<F.dimension();j++) {
            if (j==var_index) MR(i,j)= Realzeros(i);
            else MR(i,j)= -F(j)(Realzeros(i));
        }
    }
    for (i=0;i<NC;i++) {
        for (j=0;j<F.dimension();j++) {
            if (j==var_index) MC(i,j)= Compzeros(i);
            else MC(i,j)= -F(j)(Compzeros(i));
        }
    }
    if (comment) WriteText("\n number of roots: " + val::ToString(N));

    return N;
}


int computeroots(val::Glist<val::s_polynom<val::integer>> &G,matrix<double> &MR,matrix<complex> &MC,int nthreads,int comment)
{
    if (G.isempty() || !val::iszerodimensional(G)) {
        if (comment) WriteText("\n Basis is not zero dimensional!");
        return 0;
    }

    int i,j,N=0,NR,NC,index,n=s_expo::getdim(),maxdgree=0,order=s_expo::getordtype();
    val::vector<double> r_roots;
    val::vector<val::complex> c_roots;
    val::vector<pol<double>> F;
    val::matrix<int> OM;
    integer nofzeros;
    if (comment) WriteText("\nCheck if Basis is radical:");
    index = zero_dim_radical_ideal(G,maxdgree,nthreads,comment);
    //
    if (comment) WriteText("\nCheck if Basis has normal position index:");
    if ((nofzeros=val::Hilbertpolynomial(G).LC().nominator())==integer(maxdgree)) {
        if (comment) WriteText("\nBasis has normal position index = " + val::ToString(index));
        if (index!=n-1 || order!=-1) {
            if (comment) WriteText("\n Convert Basis:");
            if (index==n-1) {
                val::groebnerwalk(G,-1,n-1,n-1,OM,0);
            }
            else {
                val::vector<int> selection(0,n);
                selection(index) = 1;
                OM=id_roots::eliminationmatrix(selection);
                val::groebnerwalk(G,-1000,n-1,n-1,OM,0);
            }
            if (comment) WriteText("\nBasis converted!");
        }
        if (comment) WriteText("\n Compute roots: ");
        F = id_roots::getdoublepolynomials(G,index);
        N=val::roots(F[index],r_roots,c_roots);
        //std::cout<<"\n number of roots: "<<N;
        NR = r_roots.dimension();
        NC = c_roots.dimension();
        if (NR>0) MR=matrix<double>(NR,F.dimension());
        else MR=matrix<double>();
        if (NC>0) MC=matrix<complex>(NC,F.dimension());
        else MC=matrix<complex>();

        for (i=0;i<NR;++i) {
            for (j=0;j<n;++j) {
                if (j==index) MR(i,j) = r_roots(i);
                else MR(i,j) = -F(j)(r_roots(i));
            }
        }
        for (i=0;i<NC;++i) {
            for (j=0;j<n;++j) {
                if (j==index) MC(i,j) = c_roots(i);
                else MC(i,j) = -F(j)(c_roots(i));
            }
        }
        if (comment) WriteText("\n Roots computed!");

    }
    else { //
        Glist<s_polynom<integer>> H;
        Glist<int> Primlist;
        pol<rational> h;
        val::vector<int> selection(0,n+1);

        selection(0) = 1;

        if (comment) WriteText("\nBasis has no normal position index!\nSearch extension basis:");
        id_roots::extendordermatrix(); // Now s_expp::dim = n+1;
        do {
            Primlist = CreatePrimlist(H);
            H = id_roots::preparenormalposition(G);
            h = modint_minimalpolynom(H,0,Primlist,nthreads,0);
        }
        while (integer(h.degree())<nofzeros);
        if (comment) WriteText("\n Extension basis found!");

        if (comment) WriteText("\nConvert Basis: ");
        OM = id_roots::eliminationmatrix(selection);
        val::groebnerwalk(H,-1000,n,n,OM,0);
        if (comment) WriteText("\n Basis converted!");

        index = 0;
        if (comment) WriteText("\n Compute roots: ");
        F = id_roots::getdoublepolynomials(H,index);
        N=val::roots(F[index],r_roots,c_roots);
        NR = r_roots.dimension();
        NC = c_roots.dimension();
        if (NR>0) MR=matrix<double>(NR,F.dimension()-1);
        else MR=matrix<double>();
        if (NC>0) MC=matrix<complex>(NC,F.dimension()-1);
        else MC=matrix<complex>();

        for (i=0;i<NR;++i) {
            for (j=0;j<n;++j) {
                MR(i,j) = -F(j+1)(r_roots(i));
            }
        }
        for (i=0;i<NC;++i) {
            for (j=0;j<n;++j) {
                MC(i,j) = -F(j+1)(c_roots(i));
            }
        }
        if (comment) WriteText("\n Roots computed!");
        s_expo::setdim(n);
    }

    if (comment) WriteText(" Number of roots: " + val::ToString(N));

    return N;
}


} //end namespace val
