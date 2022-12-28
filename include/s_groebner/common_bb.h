#ifndef COMMON_BB_H_INCLUDED
#define COMMON_BB_H_INCLUDED

#include <string>
#include <s_polynom.h>
#include <rational.h>
#include <modq.h>
#include <Glist.h>
#include <atomic>
#include <mutex>
#include <vector.h>
#include <val_utils.h>
#include <thread>
//#include <MyTime.h>

namespace common_bb
{

extern int MaxThreads,ComputingThreads;
extern std::string primlistpath;

void SetComputingThreads(int n);

extern std::mutex ProtectData,ProtectOutput;

//
struct spair
{
 val::s_polynom<val::integer> *fint,*gint;
 val::s_polynom<val::modq> *fmodq,*gmodq;
 int s_done;
 int deg;
 spair(){fint=NULL;gint=NULL;fmodq=NULL;gmodq=NULL;deg=0;s_done=0;}
 spair(val::s_polynom<val::modq>* f,val::s_polynom<val::modq>* g )
      {fmodq=f;gmodq=g;fint=gint=NULL;s_done=deg=0;}
 int operator <(const spair &s) const
     {
		 if (deg<s.deg) return 1;
		 if (deg>s.deg) return 0;
		 if (fint!=NULL && gint!=NULL) {
            return (val::lcm(fint->LT(),gint->LT()) < val::lcm(s.fint->LT(),s.gint->LT()));
		 }
		 else if (fmodq!=NULL && gmodq!=NULL) {
            return (val::lcm(fmodq->LT(),gmodq->LT()) < val::lcm(s.fmodq->LT(),s.gmodq->LT()));
		 }
		 else return 0;
	 }
};


template<class T>
struct Pairs
{
    const val::s_polynom<T> *f1=NULL,*f2=NULL;
    int s_done=0;
    val::integer deg;
    int operator <(const Pairs<T> &p) const
    {
        if (deg<p.deg) return 1;
        else if (p.deg<deg) return 0;
        else return 0;
    }
    static std::atomic<int> anzs,is,anzeige;
};

template <class T> std::atomic<int> Pairs<T>::anzs(0);
template <class T> std::atomic<int> Pairs<T>::is(1);
template <class T> std::atomic<int> Pairs<T>::anzeige(0);


//

struct KritPairs
{
    val::s_polynom<val::integer> *if1,*if2,*ih;
    val::s_polynom<val::modq> *mf1,*mf2,*mh;
    int s_done;
    static std::atomic<int> nGd,p_teilt;
    static std::mutex protectGd;
};

//

DLL_PUBLIC void MyMessage(const std::string&);
DLL_PUBLIC void Message_to_cout(const std::string&);
extern DLL_PUBLIC void (*messageoutput) (const std::string&);


void DLL_PUBLIC Clear();
void DLL_PUBLIC Clear_console();
extern void DLL_PUBLIC (*clearoutput)();

void DLL_PUBLIC WriteText(const std::string&);
void DLL_PUBLIC WriteText_to_cout(const std::string&);
extern void (*writetextoutput) (const std::string&);

//
int lcmdis(const val::s_expo &t1,const val::s_expo &t2,val::s_expo &t);
//


//
int gettype(const std::string& name);

int gettype(char* name);

int getnumberofvariables(const std::string& name);


template <class T>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<T> >& P,int onlytotdeghomogen=0);

template <class T>
int readfromfile(char* name,val::Glist<val::s_polynom<T> >& P,int onlytotdeghomogen=0); // Liest Polynome aus Datei name

template <class T>
void savelist(const std::string &arg,const val::Glist<val::s_polynom<T> >& P);

template <class T>
void savelist(char* arg,const val::Glist<val::s_polynom<T> >& P)
{
    savelist(std::string(arg),P);
}


void CreatePrimelist(val::Glist<int> &Primelist,const val::Glist< val::s_polynom<val::integer> > &F,int r=-1);

int getnextprime(const val::Glist< val::s_polynom<val::integer> > &P,val::Glist<int> &Primlist,
                 const val::integer &wert,const val::Glist<val::s_polynom<val::integer> > &Gd = val::Glist<val::s_polynom<val::integer> >());

//

int Maximalinteger(const val::Glist< val::s_polynom<val::integer> >& P);

//
void restoreGp(val::Glist< val::s_polynom<val::modq> > &Gp,const val::Glist< val::s_polynom<val::integer> > &G);

//

template <class T>
void interredBasis(val::Glist< val::s_polynom<T> > &G,int degoption=0);

template <class T>
val::s_polynom<T> spol(const val::s_polynom<T>&,const val::s_polynom<T>&);

template <class T>
int dehomred(val::Glist<val::s_polynom<T> > &G,int nG,int order,const val::matrix<int> &M);

template <class T>
int minimalGroebner(val::Glist<val::s_polynom<T> > &G);

int homogenize(const std::string& i_name,const std::string& o_name,int comment=0);

int dehomogenize(const std::string& i_name,const std::string& o_name,int comment=0);

int preparenormalpos(const std::string& i_name,const std::string& o_name,int comment=0);


// Primitive Buchberger-Algorithm.
// G is already a groebner basis and its elements are normalized. Afterward H is empty.
template <class T>
int primGroebner(val::Glist<val::s_polynom<T> > &G,val::Glist<val::s_polynom<T> > &H);
//


template <class T>
int Createcritpairs(const val::Glist< val::s_polynom<T> > &G,val::matrix<char> &s_done);

template <class T>
void SimultanGroebner(val::Glist<Pairs<T> > &P,const val::Glist<val::s_polynom<T> > &G);

/*
template <class T>
void SimultanGroebner2(val::Glist<Pairs<T> > &P,const val::Glist<val::s_polynom<T> > &G);
*/

template <class T>
int isgroebner(const val::Glist<val::s_polynom<T> > &G,int mG=-1,int anzeige=1);


template <class T>
void SimultanisinG(val::Glist< val::s_polynom<T> > &H,const val::Glist< val::s_polynom<T> > &G);

template <class T>
int isinG(val::Glist< val::s_polynom<T> > &H,const val::Glist< val::s_polynom<T> > &G,int anzeige=0);

int isinG(const std::string& file,const std::string& gbfile);

int isinG(char* file,char* gbfile);

template <class T>
void reduceKritPairs(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<T> > &G,
                     val::Glist<val::s_polynom<T> > &Gd,int nd=0,int checkdiv=0);

void reduceKritPairs2(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::integer> > &G,
                      val::Glist<val::s_polynom<val::integer> > &Gd,
                      val::Glist<val::s_polynom<val::integer> > &Hd,int nd);


/*
// sequential : 1 thread
void reduceKritPairs_seq(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::integer> > &G,
                      val::Glist<val::s_polynom<val::integer> > &Gd,int nd);
*/

template <class T>
void update(val::s_polynom<T> &f,val::Glist< val::s_polynom<T> > &G,val::Glist< common_bb::spair > &lspair,int m);

template <class T>
void primupdate(val::s_polynom<T> &f,val::Glist<val::s_polynom<T> > &G,val::Glist<Pairs<T> > &PairsList);

void update(val::s_polynom<val::integer> &fint,val::s_polynom<val::modq> &fmodq,val::Glist< val::s_polynom<val::integer> > &G,
            val::Glist< val::s_polynom<val::modq> > &Gp,val::Glist<spair> &lspair,int m);

template <class T>
void update(val::s_polynom<T> &f,val::Glist<val::s_polynom<T> > &G,val::Glist<Pairs<T> > &lspair,int m,const val::vector<val::integer> &w);




template <class T>
void updateG(val::Glist< val::s_polynom<T> > &G,int &nG,val::Glist< val::s_polynom<T> > &Gd,val::Glist<common_bb::spair> &lspair);


//

template <class T>
void interredBasis(val::Glist< val::s_polynom<T> > &G,int degoption)
{

 if (G.isempty()) return;


 G.resetactual();
 G.moveactual();
 for (;G.actualvalid();G.moveactual()) {
	 G.actualvalue().reduction(G,0,1,degoption);
 }
 G.resetactual();
 return;
}


template <class T>
int dehomred(val::Glist<val::s_polynom<T> > &G,int nG,int order,const val::matrix<int> &M)
{
 int nH=0,teilt,n=val::s_expo::getdim();
 val::Glist<val::s_polynom<T> > H;
 //listelement< polynom<integer> > *pH;
 val::GlistIterator<val::s_polynom<T> > pH;

 val::s_expo::setdim(n-1);
 n--;
 val::s_expo::setordtype(order);
 val::s_expo::setordmatrix(M);

 if (nG <=1) return nG;

 G.sort(nG);

 G.moveheadtoend(H);
 nH++;

 while (!G.isempty()) {//(G.head!=NULL) {
	 G.resetactual();
	 teilt=0;
	 for (pH=H;pH;pH++) {
		 if (pH().LT() | G.actualvalue().LT()) {
			 teilt=1;
			 break;
		 }
	 }
	 if (teilt) {
		 G.skiphead();
	 }
	 else {
		 G.moveheadtoend(H);
		 nH++;
	 }
 }
 G=std::move(H);
 common_bb::interredBasis(G);
 return nH;
}


template <class T>
int minimalGroebner(val::Glist<val::s_polynom<T> > &G)
{
    using namespace val;
    int n=0,red=0;

    if (G.isempty()) return n;

    Glist<s_polynom<T> > H;
    GlistIterator<s_polynom<T> > pH;

    G.resetactual();
    if (!G.actualvalue().iszero()) {
        H.inserttoend(std::move(G.actualvalue()));
        G.skiphead();
        n++;
    }

    while (!G.isempty()) {
        red=1;
        if (!G.actualvalue().iszero()) {
            red=0;
            for (pH=H;pH;pH++)
                if (pH().LT() | G.actualvalue().LT()) {
                    red=1;
                    break;
                }
        }
        if (!red) {
           H.inserttoend(std::move(G.actualvalue()));
           n++;
        }
        G.skiphead();
    }
    G=std::move(H);
    return n;
}


// Primitiver Buchberger Alg. Dabei bilden Elemente in G bereits eine Gröbner Basis; und diese
// sind normiert. Danach in H leeer!
template <class T>
int primGroebner(val::Glist<val::s_polynom<T> > &G,val::Glist<val::s_polynom<T> > &H)
{
{
 int i;
 val::s_polynom<T> h;
 val::Glist< Pairs<T> > spair;

 if (G.isempty() && H.isempty()) return 0;

 H.resetactual();

 while (!H.isempty() || !spair.isempty()) {
	 if (!H.isempty()) {
		 h=std::move(H.actualvalue());
		 H.skiphead();
		 h.reduction(G,0);
		 if (!h.iszero()) primupdate(h,G,spair);
		 continue;
	 }
	 if (!spair.isempty()) {

		 h=spol(*(spair.actualvalue().f1),*(spair.actualvalue().f2));
		 spair.skiphead();
		 h.reduction(G,0);
		 if (!h.iszero()) primupdate(h,G,spair);
	 }
 }

 G.sort();
 i=minimalGroebner(G);
 interredBasis(G);
 return i;
}
}

//

template <class T>
int Createcritpairs(const val::Glist< val::s_polynom<T> > &G,val::matrix<int> &s_done)
{
    int i,j,m,anzs=0;
    val::s_expo t1;
    const val::s_polynom<T> *h;
    val::GlistIterator<val::s_polynom<T> > pG,p,q;

    if (G.isempty()) return 0;
    pG=G;
    if (!pG.nextvalid()) return 0;
    pG++;

    for (m=1;pG;pG++,m++) {
        for (i=0;i<m;i++) s_done[m][i]=0;
        anzs+=m;

        h=&(pG());

        for (i=0,p=G;i<m;i++,p++) {
            if (s_done[m][i]) continue;
            if (lcmdis(h->LT(),p().LT(),t1)) {
                s_done[m][i]=2;anzs--;
                continue;
            }
            for (j=0,q=G;j<m;j++,q++) {
                if (j==i || s_done[m][j]==1) continue;
                if (val::lcm(h->LT(),q().LT())|t1) {
                    s_done[m][i]=1;anzs--;
                    break;
                }
            }
        }

        for (i=0,p=G;i<m;i++,p++)
            for (j=i+1,q=p,q++;j<m;j++,q++) {
                if (s_done[j][i]) continue;
                t1 = val::lcm(p().LT(),q().LT());
                if ( (h->LT()|t1) && (val::lcm(p().LT(),h->LT())!=t1)
                    && (val::lcm(h->LT(),q().LT())!=t1) ) {
                        s_done[j][i]=1;anzs--;
                }
            }
    }

    return anzs;
}



template <class T>
void SimultanGroebner(val::Glist<Pairs<T> > &P,const val::Glist<val::s_polynom<T> > &G)
{
    const val::s_polynom<T> *f,*g;
    val::s_polynom<T> h;

    do {
        ProtectData.lock();
        if (!P.actualvalid()) {
            ProtectData.unlock();
            return;
        }
        f=P.actualvalue().f1;
        g=P.actualvalue().f2;
        P.moveactual();
        ProtectData.unlock();
        if (!Pairs<T>::is) return;
        h = spol(*f,*g);
        if (!Pairs<T>::is) return;
        h.reduction(G);
        if (!h.iszero()) {
            Pairs<T>::is=0;
            return;
        }
        Pairs<T>::anzs--;
        if (Pairs<T>::anzeige && Pairs<T>::anzs%10==0) {
            ProtectOutput.lock();
            Clear();
            WriteText("Remaining s-pols: " + val::ToString(Pairs<T>::anzs) + "  ");
            ProtectOutput.unlock();
        }
    }
    while(1);
}


/*
template <class T>
void SimultanGroebner2(val::Glist<Pairs<T> > &P,const val::Glist<val::s_polynom<T> > &G)
{
    val::s_polynom<T> h;
    val::GlistIterator<Pairs<T> > p;

    for (p=P;p;p++) {
        if (!Pairs<T>::is) return;
        h = spol(*(p().f1),*(p().f2));
        if (!Pairs<T>::is) return;
        h.reduction(G);
        if (!h.iszero()) {
            Pairs<T>::is=0;
            return;
        }
        Pairs<T>::anzs--;
        if (Pairs<T>::anzeige && Pairs<T>::anzs%10==0) {
            ProtectOutput.lock();
            Clear();
            WriteText("Noch zu pruefende S-Polynome: " + val::ToString(Pairs<T>::anzs) + "  ");
            ProtectOutput.unlock();
        }
    }
}
*/





template <class T>
int isgroebner(const val::Glist<val::s_polynom<T> > &G,int mG,int anzeige)
{
 if (mG==-1) mG=G.length();
 if (mG <= 1) return 1;

 int anzs,i,j;//k=0;
 val::matrix<int> s_done(0,mG,mG);
 int nThreads=ComputingThreads;
 val::Glist<Pairs<T> > P;

 val::GlistIterator<val::s_polynom<T> > p,q;
 Pairs<T> s;

 anzs=Createcritpairs(G,s_done);

 Pairs<T>::anzs=anzs;
 Pairs<T>::anzeige=anzeige;
 Pairs<T>::is=1;

 for (p=G,i=0;p;p++,i++) {
    q=p;q++;
    for (j=i+1;q;q++,j++) {
        if (s_done[j][i]) continue;
        s.f1 = &(p());s.f2=&(q());
        P.inserttoend(s);
    }
 }

 val::vector<std::thread* > Thr(nThreads);

 for (i=0;i<nThreads;i++) Thr(i)= new std::thread(SimultanGroebner<T>,std::ref(P),std::cref(G));
 for (i=0;i<nThreads;i++) Thr(i)->join();

 for (i=0;i<nThreads;i++) delete Thr[i];
 if (anzeige) WriteText("\nNumber of critical pairs: " + val::ToString(anzs));

 return Pairs<T>::is;
}


template <class T>
void SimultanisinG(val::Glist< val::s_polynom<T> > &H,const val::Glist< val::s_polynom<T> > &G)
{
    val::s_polynom<T> *h;

    do {
        ProtectData.lock();
        if (!H.actualvalid()) {
            ProtectData.unlock();
            return;
        }
        h=&(H.actualvalue());
        H.moveactual();
        ProtectData.unlock();
        if (!Pairs<T>::is) return;
        h->reduction(G,1);
        if (!h->iszero()) {
            Pairs<T>::is=0;
            return;
        }
        Pairs<T>::anzs--;
        if (Pairs<T>::anzeige && Pairs<T>::anzs%5==0) {
            ProtectOutput.lock();
            Clear();
            WriteText("\nRemaining s-pols: " + val::ToString(Pairs<T>::anzs));
            ProtectOutput.unlock();
        }
    } while(1);
}


// Test if <H> is subset of <G>:
template <class T>
int isinG(val::Glist< val::s_polynom<T> > &H,const val::Glist< val::s_polynom<T> > &G,int anzeige)
{
    int i;
    val::vector<std::thread*> Thr(ComputingThreads);

    Pairs<T>::is=1;

    if (anzeige) {
        Pairs<T>::anzeige=1;
        Pairs<T>::anzs=H.length();
    }
    else Pairs<T>::anzeige=0;

    H.resetactual();

    for (i=0;i<ComputingThreads;i++) Thr(i) = new std::thread(SimultanisinG<T>,std::ref(H),std::cref(G));
    for (i=0;i<ComputingThreads;i++) Thr(i)->join();
    for (i=0;i<ComputingThreads;i++) delete Thr(i);
    return Pairs<T>::is;
}


//

template <class T>
void primupdate(val::s_polynom<T> &f,val::Glist<val::s_polynom<T> > &G,val::Glist<Pairs<T> > &PairsList)
{
    using namespace val;
    if (G.isempty()) {
        G.inserttoend(std::move(f));
        return;
    }

    int ismonom=f.ismonomial();
    s_polynom<T> *g;
    Pairs<T> pair;
    GlistManipulator<s_polynom<T> > pG;
    s_expo t;

    G.inserttoend(std::move(f));
    pG.settolast(G);
    g = &(pG());

    pair.f1=g;

    for (pG=G;pG;pG++) {
        if (pG.ispointtolast(G)) break;
        //Monom-criterium
        if (ismonom && pG().ismonomial()) continue;
        // BB 1. criterium
        else if (lcmdis(g->LT(),pG().LT(),t)) continue;
        else {
            pair.f2=&(pG());
            PairsList.inserttoend(pair);
        }
    }

}



template <class T>
void update(val::s_polynom<T> &f,val::Glist<val::s_polynom<T> > &G,val::Glist<Pairs<T> > &lspair,int m,const val::vector<val::integer> &w)
{
 int i,j;
 val::s_expo t1;
 char *h_done;
 val::s_polynom<T> *hint;

 val::GlistManipulator<val::s_polynom<T> > p,q;
 val::GlistManipulator<Pairs<T> > sp;
 Pairs<T> hs;



 if (f.iszero()) return;

 if (G.isempty()) return;

 h_done = new char[m];
 for (i=0;i<m;i++) h_done[i]=0;

 hint= &f;

 for (i=0,p=G;i<m;i++,p++) {
	 if (h_done[i]) continue;
	 if (lcmdis(hint->LT(),p().LT(),t1)) {
		 h_done[i]=2;
		 continue;
	 }
	 for (j=0,q=G;j<m;j++,q++) {
		 if (j==i || h_done[j]==1) continue;
		 if (val::lcm(hint->LT(),q().LT())|t1) {
			 h_done[i]=1;
			 break;
		 }
	 }
 }


 for (sp=lspair;sp;sp++) {
	 if (sp().s_done) continue;
	 t1 = val::lcm(sp().f1->LT(),sp().f2->LT());

	 if ( (hint->LT()|t1) && (val::lcm(sp().f1->LT(),hint->LT())!=t1)
		 && (val::lcm(hint->LT(),sp().f2->LT())!=t1) ) {
		 sp().s_done=1;
	 }
 }



  // Insert in dPairs:
 for (i=0,p=G;i<m;i++,p++) {
	 if (h_done[i]) continue;
	 hs.f1=&(p());hs.f2=hint;
     hs.deg= val::s_polynom<T>::degw(val::lcm(p().LT(),hint->LT()),w);
	 lspair.sinsert(hs);
 }
 delete[] h_done;
}


template <class T>
void updateG(val::Glist< val::s_polynom<T> > &G,int &nG,val::Glist< val::s_polynom<T> > &Gd,val::Glist<common_bb::spair> &lspair)
{
 Gd.resetactual();
 while (!Gd.isempty()) {
        update(Gd.actualvalue(),G,lspair,nG);
        Gd.moveheadtoend(G);
        nG++;
 }
}





} // end namespace common_bb


#endif // COMMON_BB_H_INCLUDED
