#ifndef IDEAL_ROOTS_H_INCLUDED
#define IDEAL_ROOTS_H_INCLUDED

#include <s_groebner.h>
#include <Glist.h>
#include <s_polynom.h>


namespace val
{

template <class T> class complex_type;
typedef complex_type<double> complex;

DLL_PUBLIC void SetRootsMessage(strg_function&);

DLL_PUBLIC int idealroots(const std::string &filename,matrix<double> &RZMat,matrix<complex> &CZMat,int nthreads=3,int comment=1);

DLL_PUBLIC int computeroots(val::Glist<val::s_polynom<val::integer>> &G,matrix<double> &RZMat,matrix<complex> &CZMat,int nthreads=3,int comment=1);


template<class T>
int iszerodimensional(const val::vector<val::s_polynom<T> > &G,val::s_expo &degX);

template <class T>
int iszerodimensional(const val::Glist<val::s_polynom<T>> &G);


DLL_PUBLIC val::pol<val::rational> modint_minimalpolynom(const val::Glist<val::s_polynom<val::integer> > &Gint,int k,
                                              val::Glist<int>& Primlist,int nthreads,int comment=1);

// returns variable of univariate polynomial with maximal degree
DLL_PUBLIC int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::integer>>& G,int nthreads=3,int comment=1);

DLL_PUBLIC int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::integer>>& G,int &degree,int nthreads=3,int comment=1);

DLL_PUBLIC int zero_dim_radical_ideal(val::Glist<val::s_polynom<val::modq>>& G,int nthreads=3,int comment=1);

template <class T>
int radical_ideal(val::Glist<val::s_polynom<T>>& G,const int nthreads=3,const int comment=1);




template<class T>
int iszerodimensional(const val::vector<val::s_polynom<T> > &G,val::s_expo &degX)
{
    using namespace val;
    int k,i,j,n=s_expo::getdim(),is=1;
    //listelement<polynom<T> > *ptoG;

    if (G.isempty()) return 0;

    //degX=val::vector<int>(0,n);
    degX=val::s_expo(0);

    //for (ptoG=G.head;ptoG!=NULL;ptoG=ptoG->next) {
    for (k=0;k<G.dimension();k++) {
        for (i=0;i<n;i++) {
            if (degX[i]) continue;
            else if (G(k).LT()[i]!=0) {
                degX[i]=G(k).LT()[i];
                for (j=0;j<n;j++) {
                    if (j==i) continue;
                    if (G(k).LT()[j]!=0) {
                        degX[i]=0;
                        break;
                    }
                }
            }
        }
    }
    for (i=0;i<n;i++) {
        if (!degX[i]) {
            is=0;
            break;
        }
    }
    return is;
}




template <class T>
int iszerodimensional(const val::Glist<val::s_polynom<T>> &G)
{
    if (G.isempty()) return 0;

    int is,i,j,n=val::s_expo::getdim();

    val::s_expo X(0);

    for (const auto& f : G) {
        is = 1;
        for (i=0;i<n;++i) {
            if (f.LT()[i]) break;
        }
        if (X[i]) continue;
        j=i;
        //std::cout<<j<<" , ";
        for (i++;i<n;++i) {
            if (f.LT()[i]) {
                is=0;
                //std::cout<<" is ";
                break;
            }
        }
        if (is) X[j] = 1;
    }

    //std::cout<<"\n X = "<<X;

    for (i=0;i<n;++i) {
        if (!(X[i])) return 0;
    }

    return 1;
}

template <class T>
int radical_ideal(val::Glist<val::s_polynom<T>>& G,const int nthreads,const int comment)
{
    if (iszerodimensional(G)) {
        if (comment) WriteText("\n <G> is zerodimensional.");
        return zero_dim_radical_ideal(G,nthreads,comment);
    }
    else {
        if (comment) WriteText("\n <G> is not zerodimensional!");
    }
    //return zero_dim_radical_ideal(G,nthreads)
    return -1;
}


} // end namespace val


#endif // IDEAL_ROOTS_H_INCLUDED
