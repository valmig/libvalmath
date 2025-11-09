#ifndef BBHOM_H_INCLUDED
#define BBHOM_H_INCLUDED

#include <matrix.h>
#include <rational.h>
#include <s_groebner/common_bb.h>
#include <MyTime.h>
#include <thread>
#include <s_groebner/hilbert.h>


namespace bbhom
{
// terms have an additional variable to eventually homogenize them.
template <class T>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<T> > &P,int &wastotdegcompatible,int &washomogen);

int gauss(val::matrix<val::rational>& a,int m=-1);

void perturbmatrix(val::matrix<int> &M);

template <class T>
int totdegisinitialform(const val::s_polynom<T> &f);


template <class T>
void checktotaldegreecriteria(const val::Glist<val::s_polynom<T>> &, int &washomogen, int &wastotdegcompatble);


template <class T>
void SetKP(val::s_polynom<T> &f,common_bb::KritPairs &KP);

template <class T>
void SetKP(common_bb::KritPairs &KP,common_bb::spair &SP);

// Term order has to be totdeg-compatible.
template <class T>
int Homgroebner(val::Glist<val::s_polynom<T> > &G,int comment = 1);


// Elements in G are in increasing ordered wrt. degw and term-order.
template<class T>
int Groebnerpos(val::Glist<val::s_polynom<T> > &G,val::vector<val::integer> &w,const val::integer &dmin);

template <class T>
int bbhom(val::Glist<val::s_polynom<T> > &G);

template <class T>
int bbhom(val::Glist<val::s_polynom<T> > &G,int washomogen,int wastotdegcompatible);


template <class T>
int bbhommain(const std::string &name,val::Glist<val::s_polynom<T> >& G);


template <class T>
int bbhommain(val::Glist<val::s_polynom<T> >& G);



//

void updateG(val::Glist< val::s_polynom<val::integer> > &G,int &nG,val::Glist< val::s_polynom<val::integer> > &Gd,
             val::Glist<val::s_polynom<val::modq> > &Gp,val::Glist<val::s_polynom<val::modq> > &Gpd,val::Glist<common_bb::spair> &lspair,
             val::pol<val::integer>& hilbnumG);


void updateG(val::Glist< val::s_polynom<val::modq> > &G,int &nG,val::Glist< val::s_polynom<val::modq> > &Gd,
             val::Glist<common_bb::spair> &lspair,val::pol<val::integer>& hilbnumG);

template <class T>
void updateG(val::Glist<val::s_polynom<T> > &G,int &nG,val::Glist<val::s_polynom<T> > &Gd,val::Glist<common_bb::Pairs<T> > &lspair,
             const val::vector<val::integer> &w);

// Here M is a totdeg-compatible order-matrix.
template <class T>
int GroebnerwithHilbert(val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M,int nH=-1);


template <class T>
int hilbertconversion(val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M=val::matrix<int>());


template <class T>
int hilbertconversion(val::Glist<val::s_polynom<T> > &G,int washomogen,int wastotdegcompatible,int order,const val::matrix<int> &M=val::matrix<int>());

template <class T>
int hilbertconversionmain(const std::string &name,val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M=val::matrix<int>());

template <class T>
int hilbertconversionmain(val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M=val::matrix<int>());



//

template <class T>
int totdegisinitialform(const val::s_polynom<T> &f)
{
    int deg;
    val::s_polynomIterator<T> pf;

    if (f.iszero()) return 1;
    pf=f;
    deg=pf.actualterm().totdeg();
    pf++;
    for (;pf;pf++) if (pf.actualterm().totdeg() > deg) return 0;

    return 1;
}


template <class T>
void checktotaldegreecriteria(const val::Glist<val::s_polynom<T>> &G, int &washomogen, int &wastotdegcompatible)
{
    int deg, order = val::s_expo::getordtype(), tdegisinit = 1;
    val::s_polynomIterator<T> pf;

    washomogen = 1; wastotdegcompatible = 1;

    if (order == -1000) {
        const val::matrix<int> &M = val::s_expo::getordmatrix();
        int n = M.numberofcolumns();
        if (n == 0) wastotdegcompatible = 0;
        else wastotdegcompatible = 1;

        for (int i = 0; i < n; ++i) {
            if (M(0,i) != 1) wastotdegcompatible = 0;
        }
    }
    else if (order == -1) {
        wastotdegcompatible = 0;
    }
    else if (order != -2 && order != 0) val::s_expo::setordtype(-2);

    for (const auto &f : G) {
        pf = f;
        deg = pf.actualterm().totdeg();
        for (pf++; pf; pf++) {
            if (pf.actualterm().totdeg() != deg) {
                washomogen = 0;
                break;
            }
        }
        if (!washomogen && !wastotdegcompatible) {
            if (tdegisinit) tdegisinit=totdegisinitialform(f);
        }
    }
	if (tdegisinit) wastotdegcompatible = 1;
}
//

template <class T>
void updateG(val::Glist<val::s_polynom<T> > &G,int &nG,val::Glist<val::s_polynom<T> > &Gd,val::Glist<common_bb::Pairs<T> > &lspair,
             const val::vector<val::integer> &w)
{
 Gd.resetactual();
 while (!Gd.isempty()) {
        common_bb::update(Gd.actualvalue(),G,lspair,nG,w);
        Gd.moveheadtoend(G);
        nG++;
 }
}

//
template <class T>
int Homgroebner(val::Glist<val::s_polynom<T> > &G, int comment)
{
 int nG=0,dh,ds,d,nThreads=common_bb::ComputingThreads,i,anzs; //nd
 val::Glist<val::s_polynom<T> > H,Gd;
 val::Glist<common_bb::spair> Pair;
 val::Glist<common_bb::KritPairs> ListKPairs;
 common_bb::KritPairs KP;
 val::vector<std::thread*> Thr(nThreads);

 if (G.isempty()) return nG;
 H=std::move(G);

 H.resetactual();
 d=H.actualvalue().LT().totdeg();

 //nd=0;

 while (!H.isempty() && H.actualvalue().LT().totdeg()==d) {
     H.actualvalue().reduction(Gd);
     if (!H.actualvalue().iszero()) {
         Gd.sinsert(std::move(H.actualvalue()));
         //nd++;
     }
     H.skiphead();
 }

 common_bb::interredBasis(Gd);

 common_bb::updateG(G,nG,Gd,Pair);

 while (H.actualvalid() || !Pair.isempty()) {
     Pair.resetactual();
     if (H.actualvalid() && !H.actualvalue().iszero()) dh=H.actualvalue().LT().totdeg();
     else dh=-1;
     if (Pair.actualvalid()) ds=Pair.actualvalue().deg;
     else ds=-1;
     if (dh==-1) d=ds;
     else if (ds==-1) d=dh;
     else if (dh<=ds) d=dh;
     else d=ds;
     //-------------------------------------------------------------------------
     //common_bb::WriteText("\nDegree: " + val::ToString(d) + " s-pols: " + val::ToString(anzs));
     common_bb::KritPairs::nGd=0;
     KP.s_done=0;
     anzs=0;
     for (;H.actualvalid() && H.actualvalue().LT().totdeg()==d;H.moveactual()) {
         SetKP(H.actualvalue(),KP);
         ListKPairs.inserttoend(KP);
     }
     while(!Pair.isempty() && Pair.actualvalue().deg==d) {
         if (Pair.actualvalue().s_done) {
             Pair.skiphead();continue;
         }
         SetKP<T>(KP,Pair.actualvalue());
         ListKPairs.inserttoend(KP);
         anzs++;
         Pair.skiphead();
     }
     if (comment) common_bb::WriteText("\nDegree: " + val::ToString(d) + " s-pols: " + val::ToString(anzs));

     ListKPairs.resetactual();
     for (i=0;i<nThreads;i++)
         Thr(i) = new std::thread(common_bb::reduceKritPairs<T>,std::ref(ListKPairs),std::ref(G),std::ref(Gd),0,0);
     for (i=0;i<nThreads;i++) Thr(i)->join();

     for (i=0;i<nThreads;i++) delete Thr(i);
     ListKPairs.dellist();
     common_bb::interredBasis(Gd);

     common_bb::updateG(G,nG,Gd,Pair);
     if (comment) common_bb::WriteText("  . Elements: " + val::ToString(nG) + "    ");
 }

 return nG;
}


template <class T>
int bbhom(val::Glist<val::s_polynom<T> > &G)
{
    int perturbed=0,order=val::s_expo::getordtype(),nG=0;
    val::matrix<int> M;
    val::ChronoClass Chrono;
    std::string s;

    if (G.isempty()) return 0;

    val::s_polynom<T>::nreduction=0;
    if (order==-1) {
        val::s_expo::setordtype(0);
        perturbed=1;
    }
    else if (order==-1000) {
        perturbed=1;
        M=val::s_expo::getordmatrix();
        perturbmatrix(val::s_expo::getordmatrix());
    }
    if (perturbed) G.sort();
    Chrono();
    nG=Homgroebner(G);
    s+="\n\nTime: " + val::ToString(Chrono());
    s+="\nElements in G: " + val::ToString(nG);
    s+="\nMonomials: " + val::ToString(val::s_polynom<T>::getmnumber());
    s+="\nReductions: " + val::ToString(val::s_polynom<T>::nreduction);
    common_bb::WriteText(s);

    if (perturbed) {
        val::s_expo::getordmatrix()=std::move(M);
        val::s_expo::setordtype(order);
        G.sort();
    }

    val::s_polynom<T>::nreduction=0;
    return nG;
}


template<class T>
int Groebnerpos(val::Glist<val::s_polynom<T> > &G,val::vector<val::integer> &w,const val::integer &dmin)
{
 using namespace val;
 int nG=0;
 integer d,dh,ds;

 s_polynom<T> h;
 Glist<val::s_polynom<T>> H,Gd;
 Glist<common_bb::Pairs<T> > spair;

 if (G.isempty()) return 0;
 H=std::move(G);
 H.resetactual();

 d=H.actualvalue().getLTdegree();

 if (d<dmin) {
    while (!H.isempty() && H.actualvalue().getLTdegree() < dmin) {
        H.actualvalue().setLTzerodegree();
        G.inserttoend(std::move(H.actualvalue()));
        nG++;
        H.skiphead();
    }
 }

 while (!H.isempty() || !spair.isempty()) {
     if (H.actualvalid() && !H.actualvalue().iszero()) dh=H.actualvalue().getLTdegree();
     else dh=integer(-1);
     if (spair.actualvalid()) ds=spair.actualvalue().deg;
     else ds=integer(-1);
     if (dh==-1) d=ds;
     else if (ds==-1) d=dh;
     else if (dh<=ds) d=dh;
     else d=ds;

     while (!H.isempty() && H.actualvalue().getLTdegree()==d) {
             H.actualvalue().setLTzerodegree();
            H.actualvalue().reduction(G,0);
            H.actualvalue().reduction(Gd);
            if (!H.actualvalue().iszero()) {
                Gd.sinsert(std::move(H.actualvalue()));
            }
            H.skiphead();
     }
     while (!spair.isempty() && spair.actualvalue().deg==d) {
        if (spair.actualvalue().s_done) {
            spair.skiphead();
            continue;
        }
        h= common_bb::spol(*(spair.actualvalue().f1),*(spair.actualvalue().f2));
        h.reduction(G,0);
        h.reduction(Gd);
        if (!h.iszero()) {
            Gd.sinsert(std::move(h));
        }
        spair.skiphead();
     }

     common_bb::interredBasis(Gd);

     updateG(G,nG,Gd,spair,w);
     spair.resetactual();
 }
 return nG;
}


template <class T>
int hilbertconversion(val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M)
{
    int perturbed=0,nG=0,porder;
    val::matrix<int> M1;
    val::ChronoClass Chrono;
    std::string s;

    if (M.isempty() && order !=-1 && order !=-2 && order != 0) order=0;
    porder=order;

    if (G.isempty()) return nG;
    if (order==-1) {
        porder=0;
        perturbed=1;
    }
    else if (order==-1000) {
        perturbed=1;
        M1=M;
        perturbmatrix(M1);
    }

    val::s_polynom<T>::nreduction=0;
    Chrono();
    if (perturbed) nG=GroebnerwithHilbert(G,porder,M1);
    else nG=GroebnerwithHilbert(G,porder,M);
    s+="\n\nTime: " + val::ToString(Chrono());
    s+="\nElements in G: " + val::ToString(nG);
    s+="\nMonomials: " + val::ToString(val::s_polynom<T>::getmnumber());
    s+="\nReductions: " + val::ToString(val::s_polynom<T>::nreduction);
    common_bb::WriteText(s);

    if (perturbed) {
        val::s_expo::getordmatrix()=std::move(M);
        val::s_expo::setordtype(order);
        G.sort();
    }

    val::s_polynom<T>::nreduction=0;
    return nG;
}

template <class T>
int hilbertconversion(val::Glist<val::s_polynom<T> > &G,int washomogen,int wastotdegcompatible,int order,const val::matrix<int> &M)
{
    int perturbed=0,nG=0,porder;
    val::matrix<int> M1;
    val::ChronoClass Chrono;
    std::string s;

    if (M.isempty() && order !=-1 && order !=-2 && order != 0) order=0;
    porder=order;

    if (G.isempty()) return nG;

    val::s_polynom<T>::nreduction=0;

    if (washomogen) {
        val::s_expo::setdim(val::s_expo::getdim()-1);
        if (order==-1) {
            porder=0;
            perturbed=1;
        }
        else if (order==-1000) {
            perturbed=1;
            M1=M;
            perturbmatrix(M1);
        }

        Chrono();
        if (perturbed) nG=GroebnerwithHilbert(G,porder,M1);
        else nG=GroebnerwithHilbert(G,porder,M);
        s+="\n\nTime: " + val::ToString(Chrono());
        s+="\nElements in G: " + val::ToString(nG);
        s+="\nMonomials: " + val::ToString(val::s_polynom<T>::getmnumber());
        s+="\nReductions: " + val::ToString(val::s_polynom<T>::nreduction);
        common_bb::WriteText(s);

        if (perturbed) {
            val::s_expo::getordmatrix()=std::move(M);
            val::s_expo::setordtype(order);
            G.sort();
        }
    }
    else if (wastotdegcompatible) {
         //homogenize:
        int n=val::s_expo::getdim(),s_order=val::s_expo::getordtype();
        val::matrix<int> Mnew,Ms;
        if (s_order==-1000) {
            Ms=std::move(val::s_expo::getordmatrix());
            int i,j;
            Mnew = val::matrix<int>(n);
            n--;
            for (i=0;i<=n;i++) Mnew(0,i)=1;
            for (i=0;i<n;i++) {
                for (j=0;j<n;j++) {
                    Mnew(i+1,j)=Ms(i,j);
                }
                Mnew(i+1,n)=0;
            }
            n++;
            val::s_expo::setordtype(-1000);
            val::s_expo::setordmatrix(Mnew);
        }
        else if (s_order==-1) {val::s_expo::setordtype(0);}
        else if (s_order==0) {
            int i;
            val::s_expo::setordtype(-1000);
            Mnew = val::matrix<int>(0,n,n);
            n--;
            for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
            Mnew(0,n) =1;
            for (i=2;i<=n;i++) Mnew(i,i-2) =1;
            n++;
            val::s_expo::setordmatrix(Mnew);
        }
        else {val::s_expo::setordtype(-2);}
        for (G.resetactual();G.actualvalid();G.moveactual()) {G.actualvalue().homogenize();G.actualvalue().reord();}
        G.sort();
        // homogenized.

        // Set target-order:
        if (order==-1000) {
            int i,j;
            Mnew = val::matrix<int>(n);
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
            Mnew = val::matrix<int>(0,n,n);
            n--;
            for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
            Mnew(0,n) =1;
            for (i=2;i<=n;i++) Mnew(i,i-2) =1;
            n++;
        }
        else porder=-2;
        //
        Chrono();
        nG=GroebnerwithHilbert(G,porder,Mnew);
        common_bb::WriteText("\nHomogeneous Basis computed, dehomogenize and reduce!");
        nG=common_bb::dehomred(G,nG,order,M);
        s+="\n\nTime: " + val::ToString(Chrono());
        s+="\nElements in G: " + val::ToString(nG);
        s+="\nMonomials: " + val::ToString(val::s_polynom<T>::getmnumber());
        s+="\nReductions: " + val::ToString(val::s_polynom<T>::nreduction);
        common_bb::WriteText(s);
    }
    else {
        common_bb::MyMessage("\nBasis nor homogeneous nor totdeg-compatible!");
        return 0;
    }

    return nG;
}


template <class T>
int bbhommain(val::Glist<val::s_polynom<T> >& G)
{
    int nG,washomogen,wastotdegcompatible;

    checktotaldegreecriteria(G, washomogen, wastotdegcompatible);
    val::s_expo::setdim(val::s_expo::getdim() + 1);
    nG = bbhom(G,washomogen,wastotdegcompatible);

    return nG;
}




template <class T>
int hilbertconversionmain(val::Glist<val::s_polynom<T> > &G,int order,const val::matrix<int> &M)
{
    int nG,washomogen,wastotdegcompatible;

    checktotaldegreecriteria(G, washomogen, wastotdegcompatible);
    val::s_expo::setdim(val::s_expo::getdim() + 1);

    nG = hilbertconversion(G,washomogen,wastotdegcompatible,order,M);
    return nG;
}



} //end namespace


#endif // BBHOM_H_INCLUDED
