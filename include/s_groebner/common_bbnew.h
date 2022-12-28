#ifndef COMMON_BB_H_INCLUDED
#define COMMON_BB_H_INCLUDED

#include <string>
#include <Glist.h>
#include <atomic>
#include <mutex>
#include <vector.h>
#include <val_utils.h>
#include <thread>
#include <rational.h>
#include <s_polynom.h>
#include <numbers.h>
#include <fstream>

namespace common_bb
{

extern int MaxThreads,ComputingThreads;
extern std::string primlistpath;

void SetComputingThreads(int n);

extern std::mutex ProtectData,ProtectOutput;

//  --------------------------------------------------------------------------------------------------------------
template <typename T,typename S, template<typename> class poly>
struct spair
{
 poly<T> *fint = nullptr, *gint = nullptr;
 poly<S> *fmodq = nullptr, *gmodq = nullptr;
 int s_done = 0;
 int deg = 0;
 //
 spair() = default;
 spair(poly<S>* f, poly<S>* g) : fint(nullptr), gint(nullptr), fmodq(f), gmodq(g), s_done(0), deg(0) {}
 //
 int operator <(const spair<T,S,poly> &s) const
     {
		 if (deg<s.deg) return 1;
		 if (deg>s.deg) return 0;
		 if (fint!=nullptr && gint!=nullptr) {
            return (val::lcm(fint->LT(),gint->LT()) < val::lcm(s.fint->LT(),s.gint->LT()));
		 }
		 else if (fmodq!=nullptr && gmodq!=nullptr) {
            return (val::lcm(fmodq->LT(),gmodq->LT()) < val::lcm(s.fmodq->LT(),s.gmodq->LT()));
		 }
		 else return 0;
	 }
};

// --------------------------------------------------------------------------------------------------------------

template<typename T,template <typename> class poly>
struct Pairs
{
    const poly<T> *f1=NULL,*f2=NULL;
    int s_done=0;
    val::integer deg;
    int operator <(const Pairs<T,poly> &p) const
    {
        if (deg<p.deg) return 1;
        else if (p.deg<deg) return 0;
        else return 0;
    }
    static std::atomic<int> anzs,is,anzeige;
};

template <typename T,template <typename> class poly> std::atomic<int> Pairs<T,poly>::anzs(0);
template <typename T,template <typename> class poly> std::atomic<int> Pairs<T,poly>::is(1);
template <typename T,template <typename> class poly> std::atomic<int> Pairs<T,poly>::anzeige(0);


// -------------------------------------------------------------------------------------------------------------------


template <typename T,typename S, template<typename> class poly>
struct KritPairs
{
    poly<T> *if1,*if2,*ih;
    poly<S> *mf1,*mf2,*mh;
    int s_done;
    static std::atomic<int> nGd,p_teilt;
    static std::mutex protectGd;
};

// -------------------------------------------------------------------------------------------------------------------------

DLL_PUBLIC void MyMessage(const std::string&);
DLL_PUBLIC void Message_to_cout(const std::string&);
extern DLL_PUBLIC void (*messageoutput) (const std::string&);


void DLL_PUBLIC Clear();
void DLL_PUBLIC Clear_console();
extern void DLL_PUBLIC (*clearoutput)();

void DLL_PUBLIC WriteText(const std::string&);
void DLL_PUBLIC WriteText_to_cout(const std::string&);
extern void (*writetextoutput) (const std::string&);

// ---------------------------------------------------------------------------------------------------------------------------

template <class expo>
int lcmdis(const expo &t1, const expo &t2, expo &t); // defined!

// ---------------------------------------------------------------------------------------------------------------------------


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

// -------------------------------------------------------------------------------------------------------------------------------

template <template <typename> class poly>
void CreatePrimelist(val::Glist<int> &Primelist,const val::Glist<poly<val::integer> > &F,int r=-1);  // defined


template <template <typename> class poly>
int getnextprime(const val::Glist< poly<val::integer> > &P,val::Glist<int> &Primlist,
                 const val::integer &wert,const val::Glist<poly<val::integer> > &Gd = val::Glist<poly<val::integer> >()); // defined

// -------------------------------------------------------------------------------------------------------------------------------


template <template <typename> class poly>
int Maximalinteger(const val::Glist< poly<val::integer> >& P); // defined



// Definitions -------------------------------------------------------------------------------------------------------------

template <class expo>
int lcmdis(const expo &t1, const expo &t2, expo &t)
{
    int i, is = 1;

    for (i = 0; i < t1.dimension(); ++i) {
        t[i] = val::Max(t1[i],t2[i]);
        if (t1[i] && t2[i]) is = 0;
    }
    return is;
}

// ----------------------------------------------------------------------------------------------------------------------------

template <template <typename> class poly>
void CreatePrimelist(val::Glist<int> &Primelist,const val::Glist<poly<val::integer> > &F,int r)
{
 int i,j,p;
 val::integer zero;

 if (r==-1) r=F.length();
 std::ifstream file(primlistpath,std::ios::in);//ifstream file("c:\\mathe\\primlis1.txt",ios::in);

 if (!file) {
	 //val::Error::error("Datei primlist.txt nicht gefunden oder fehlerhaft!!!");
	int pp=100000;
	val::Glist<int> Primes;
	for (i=0 ;i<1000;++i) {
		pp=val::nextprime(pp);
		Primes.push_back(pp);
	}
	for (auto p : Primes) {
        for (i = 0; i < r; ++i) {
            if (F[i].LC() %val::integer(p)==zero) break;
        }
        if (i==r) {
            Primelist.inserttoend(p);
        }

	}
	file.close();
	return;
 }
 //for (i=0;i<2000;i++)
 file>>j>>p; // Die 2 wird hier nicht beruecksichtigt
 while (file) {
	 file>>j>>p;
	 for (i=0;i<r;i++)
		 if (F[i].LC() %val::integer(p)==zero) break;
	 if (i==r) {
		 Primelist.inserttoend(p);
	 }
 }
 //cout<<"\nAnzahl der Primzahlen: "<<anz;
 file.close();
}


template <template <typename> class poly>
int getnextprime(const val::Glist< poly<val::integer> > &P,val::Glist<int> &Primlist,
                 const val::integer &wert,const val::Glist<poly<val::integer> > &Gd)
{
 int p;
 val::integer zero(0);
 //listelement< polynom<integer> > *ptoP;
 val::GlistIterator<poly<val::integer> > ptoP;

 do {
	 p=Primlist.getelement();Primlist.moveactual();
	 if ((wert%val::integer(p))==zero) continue;
	 for (ptoP=P;ptoP;ptoP++) {
		 if ((ptoP().LC())%val::integer(p)==zero) break;
	 }
	 for (ptoP=Gd;ptoP;ptoP++) {
		 if ((ptoP().LC())%val::integer(p)==zero) break;
	 }
	 if (!ptoP) break;
 }
 while(1);
 return p;

}

// -----------------------------------------------------------------------------------------------------------------------

template <template <typename> class poly>
int Maximalinteger(const val::Glist< poly<val::integer> >& P)
{
 int maxint=0;
 val::GlistIterator< poly<val::integer> > pG;


 for (pG.settohead(P);pG.actualvalid();pG.moveactual()) {
	 for (const auto& p : pG.getelement()) {
		 if (p.actualcoef().abslength()>maxint) maxint=p.actualcoef().abslength();
	 }
 }

 return maxint;
}




} // end namespace common_bb




#endif
