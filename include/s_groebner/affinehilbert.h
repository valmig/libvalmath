#ifndef AFFINEHILBERT_H_INCLUDED
#define AFFINEHILBERT_H_INCLUDED

#include <s_groebner/bb_int.h>
#include <s_groebner/bbhom.h>


namespace a_hilbert
{

int readfromfile(char* name,val::Glist<val::s_polynom<val::modq> > &G,int &order,val::matrix<int> &Mold,int onlytotdegcompatible=0);


// G is homogeneous and order is deg - compatible, dim (M) = dim (s_expo -1)
template <class T>
int affinehilbertconversion(val::Glist<val::s_polynom<T> >& G,int order,const val::matrix<int> &M=val::matrix<int>());

template <class T>
int affinehilbertconversionmain(char* name,val::Glist<val::s_polynom<T> >& G,int order,const val::matrix<int> &M=val::matrix<int>());



template <class T>
int affinehilbertconversion(val::Glist<val::s_polynom<T> >& G,int order,const val::matrix<int> &M)
{
    using namespace val;
    if (G.isempty()) return 0;

    int nG, porder,n=s_expo::getdim();
    matrix<int> Mnew;
    ChronoClass Chrono;
    std::string s;

    // Set target order:

    if (order==-1000) {
        int i,j;
        Mnew = matrix<int>(n);
        n--;
        for (i=0;i<=n;i++) Mnew(0,i)=1;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				Mnew(i+1,j)=M(i,j);
			}
			Mnew(i+1,n)=0;
		}
		n++;
		porder=-1000;
    }
    else if (order==-2) porder=-2;
    else if (order==-1) porder=0;
    else if (order==0) {
        int i;
        porder=-1000;
        Mnew = matrix<int>(0,n,n);
        n--;
        for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
        Mnew(0,n) =1;
        for (i=2;i<=n;i++) Mnew(i,i-2) =1;
        n++;
    }
    else porder=-2;

    //
    val::s_polynom<T>::nreduction=0;
    Chrono();
    nG=bbhom::hilbertconversion(G,porder,Mnew);
    common_bb::WriteText("\nHomogeneous Basis computed, dehomogonize and reduce!");
    nG=common_bb::dehomred(G,nG,order,M);
    s+="\n\nTime: " + val::ToString(Chrono());
    s+="\nElements in G: " + val::ToString(nG);
    s+="\nMonomials: " + val::ToString(val::s_polynom<T>::getmnumber());
    s+="\nReductions: " + val::ToString(val::s_polynom<T>::nreduction);
    common_bb::WriteText(s);
    val::s_polynom<T>::nreduction=0;

    return nG;
}


}  // end namespace


#endif // AFFINEHILBERT_H_INCLUDED
