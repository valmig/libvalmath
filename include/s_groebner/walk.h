#ifndef WALK_H_INCLUDED
#define WALK_H_INCLUDED


#include <s_groebner/common_bb.h>
#include <s_groebner/bbhom.h>



namespace walk
{

val::matrix<int> Createordmatrix(int ord);


void perturbweight(val::vector<val::integer> &w,int deg,const val::matrix<int> &A,int k);

void setpositiv(int &wnewpos,int &wspos,int &pwtpos,const val::vector<val::integer> &ws,const val::vector<val::integer> &pwt);

int ispositv(const val::vector<val::integer> &w);



template <class T>
int totdeg(const val::Glist<val::s_polynom<T> > &G);

// Sorts elements in G wrt Mt and checks if leading terms remain equal.
template <class T>
int isinrightcone(val::Glist<val::s_polynom<T> > &G);


// Returns the number of polynomials in inwG that are monomials.
template <class T>
int Makelistsofinitials(val::Glist<val::s_polynom<T> > &G,const val::vector<val::integer> &w,val::Glist<val::s_polynom<T> > &inwG,
						val::Glist<val::s_polynom<T> > &Gnew,val::Glist<val::s_polynom<T> > &GnotMon,val::integer &dmin);


void getnextvector(val::vector<val::integer> &wold,val::vector<val::integer> &wnew,const val::rational& u);

template <class T>
void getvectornextcone(const val::Glist<val::s_polynom<T> >& G,val::rational &u,
                      val::vector<val::integer> &wold,val::vector<val::integer> &wnew,val::vector<val::integer> &wt);


template <class T>
int isGroebner(const val::Glist<val::s_polynom<T> > &G);

//

// Groebner-Lauf:
template <class T>
void GroebnerWalk(val::Glist<val::s_polynom<T> > &G,int neworder,val::matrix<int> &Mt,int k,int l,int comment=1);

template <class T>
void walkmain(val::Glist<val::s_polynom<T> > &G,int neworder,int k=2,int l=2,const val::matrix<int> &M=val::matrix<int>(),int comment=1);

template <class T>
void walkmain(const std::string& name,val::Glist<val::s_polynom<T> > &G,int neworder,int k=2,int l=2,const val::matrix<int> &M=val::matrix<int>());


//
template <class T>
void saverrorlist(const val::Glist<val::s_polynom<T> > &G,int withdegree=0)
{
    std::ofstream file("error.txt",std::ios::out | std::ios::trunc);
    val::GlistIterator<val::s_polynom<T> > pG;
    val::s_polynomIterator<T> p;
    for (pG=G;pG;pG++) {
        for (p=pG();p;p++) {
            file<<p.actualcoef()<<std::endl;
            file<<p.actualterm()<<"    ";
            if (withdegree) file<<p.actualdegree()<<std::endl;
            else file<<std::endl;
        }
        file<<"0"<<std::endl<<std::endl;
    }
}

//

template <class T>
int totdeg(const val::Glist<val::s_polynom<T> > &G)
{
 int d=0;
 val::GlistIterator<val::s_polynom<T> > pG;

 if (G.isempty()) return -1;

 for (pG=G;pG;pG++) d = val::Max(d,pG().totdeg());
 return d;
}


template <class T>
int isinrightcone(val::Glist<val::s_polynom<T> > &G)
{
 int is=1;
 val::s_expo X;

 for (G.resetactual();G.actualvalid();G.moveactual()) {
	 X = G.actualvalue().LT();
	 G.actualvalue().reord();
	 if (X != G.actualvalue().LT()) is =0;
 }
 G.resetactual();
 return is;
}


template <class T>
int Makelistsofinitials(val::Glist<val::s_polynom<T> > &G,const val::vector<val::integer> &w,val::Glist<val::s_polynom<T> > &inwG,
						val::Glist<val::s_polynom<T> > &Gnew,val::Glist<val::s_polynom<T> > &GnotMon,val::integer &dmin)
{
 int nmon=0;
 val::integer d,zero;
 val::s_polynom<T> f,h;
 val::GlistIterator<val::s_polynom<T> > pG;

 dmin=val::integer(-1);
 for (pG=G;pG;pG++) {
	 pG().makeinitforms(f,h,w);
	 inwG.inserttoend(std::move(f));
	 if (h.ismonomial()) {
         Gnew.inserttoend(std::move(h));
		 nmon++;
	 }
	 else {
		 d= val::s_polynom<T>::degw(h.LT(),w);
		 if (dmin<zero) {
			 dmin = d;
		 }
		 else if (d < dmin) dmin=d;
		 GnotMon.inserttoend(std::move(h));
	 }
 }
 return nmon;
}


template <class T>
void getvectornextcone(const val::Glist<val::s_polynom<T> >& G,val::rational &u,
                      val::vector<val::integer> &wold,val::vector<val::integer> &wnew,val::vector<val::integer> &wt)
{
 using namespace val;
 int i,n=wt.dimension();
 integer a1,a2,b;
 GlistIterator<s_polynom<T> > pG;
 s_polynomIterator<T> p,q;

 u= rational(1);
 for (i=0;i<n;i++) {
	 wold[i]=wnew[i];
	 wnew[i]=wt[i];
 }
 for (pG=G;pG;pG++) {
	 q=pG();
	 p=pG();
	 a1=s_polynom<T>::degw(p.actualterm(),wnew);
	 p++;
	 while (p) {
		 a2=s_polynom<T>::degw(p.actualterm(),wnew);
		 if (a1<a2) {
			 b=s_polynom<T>::degw(q.actualterm() / p.actualterm(),wold);
			 u=rational(b,b-(a1-a2));

			 getnextvector(wold,wnew,u);
			 a1=s_polynom<T>::degw(q.actualterm(),wnew);
		 }
		 p++;
	 }
 }
}




template <class T>
int isGroebner(const val::Glist< val::s_polynom<T> > &G)
{
 int is=1;
 val::GlistIterator<val::s_polynom<T> > pG;
 val::s_polynomIterator<T> p;

 for (pG=G;pG;pG++) {
	 p=pG();
	 p++;
	 for (;p;p++)
		 if (pG().LT() < p.actualterm()) return 0;
 }
 return is;
}


//

// 
template <class T>
void GroebnerWalk(val::Glist<val::s_polynom<T> > &G,int neworder,val::matrix<int> &Mt,int k,int l,int comment)
{
 using namespace val;
 int i,isdone=0,m=G.length(),deg=0,n=s_expo::getdim(),ncones=0,nmon,nred,oldorder=s_expo::getordtype(),mold;
 int wnewpos=0,wspos=0,pwtpos=0,anzprim=0,anzpos=0;
 val::vector<char> remain;
 double crit;
 rational u;
 GlistManipulator<s_polynom<T> > pGlist;
 Glist<s_polynom<T> > inwG,Gnew,GnotMon;
 ChronoClass Chrono;

 integer dmin;
 val::vector<integer> wold(n),wnew(n),wt(n),pwt(n);

 val::matrix<int> Ms;


 if (k>=n) k=n;
 if (k<1) k=1;
 if (l>=n) l=n;
 if (l<1) l=1;

 Chrono();
 if (s_expo::getordmatrix().isempty()) {
	 Ms=Createordmatrix(s_expo::getordtype());
 }
 else Ms = s_expo::getordmatrix();

 if (neworder != -1000) {
	 Mt=Createordmatrix(neworder);
 }


 for (i=0;i<n;i++) wt[i]= integer(Mt(0,i));

 if (k>1 || l>1) deg = totdeg(G);


 perturbweight(wnew,deg,Ms,k);
 perturbweight(pwt,deg,Mt,l);
 setpositiv(wnewpos,wspos,pwtpos,wnew,pwt);


 ncones=1;
 s_expo::setordmatrix(Mt);
 s_expo::setordtype(neworder);

 // Prüfe ob man schon im richtigen Kegel ist und berechne sonst die erste G-Basis:
 nmon=Makelistsofinitials(G,wnew,inwG,Gnew,GnotMon,dmin);
 if (!isinrightcone(GnotMon)) {
	 if (comment) common_bb::WriteText("\nInitial cone has to be treated!");
	 ncones++;
	 if (comment) common_bb::WriteText("\nNumber of cones: "+ ToString(ncones));
	 if (comment) common_bb::WriteText("\nNumber of polynomials in inwG that are monomials: " +  ToString(nmon));
	 crit = double(nmon)/double(m);
	 mold=m;
	 if (crit < 0.8 && wspos) {
		 for (pGlist=GnotMon;pGlist;pGlist++)  {
              pGlist().reord();
              pGlist().normalize();
              pGlist().setLTdegree(wnew);
		 }
		 for (pGlist=Gnew;pGlist;pGlist++) {pGlist().normalize();pGlist().setLTdegree(wnew);}
		 Gnew.sort(nmon);
		 GnotMon.sort(m-nmon);
		 Gnew.merge(GnotMon);
		 if (comment) common_bb::WriteText("\nComputing groebner basis with pos. hom. BB-algorithm!");
		 m=bbhom::Groebnerpos(Gnew,wnew,dmin);
		 if (!common_bb::isgroebner(Gnew)) {
            std::cout<<"\nIs not a groebner basis!";
            Gnew.sort();
            saverrorlist(Gnew);
            exit(-1);
		 }
		 anzpos++;

		 for (pGlist=Gnew;pGlist;pGlist++) pGlist().setLTzerodegree();
	 }
	 else {
		 if (comment) common_bb::WriteText("\nComputing groebner basis with prim. BB-Algorithmus!");
		 m=common_bb::primGroebner(Gnew,GnotMon);
		 anzprim++;
	 }
	 if (comment) common_bb::WriteText( "\nNumber of elements in Gnew = " + ToString(m));
		 // Basis Liften:
	 if (comment) {
		common_bb::WriteText("\nLifting basis");
		common_bb::WriteText("\nMaximal integer-length: " + ToString(integer::GetMaxlength()));
	 }
	 remain=val::vector<char>(0,mold);
	 s_expo::setordmatrix(Ms);s_expo::setordtype(oldorder);
	 for (Gnew.resetactual();Gnew.actualvalid();Gnew.moveactual()) {
		 Gnew.actualvalue().reord();
		 Gnew.actualvalue().liftpolynom(inwG,G,remain);
	 }
	 for (Gnew.resetactual(),G.resetactual(),i=0;Gnew.actualvalid();Gnew.moveactual()) {
        if (Gnew.actualvalue().iszero()) {
            while (G.actualvalid() && !remain[i]) {
                G.moveactual();i++;
            }
            Gnew.actualvalue()=std::move(G.actualvalue());
            G.moveactual();i++;
            if (!G.actualvalid()) break;
        }
	 }

	 G.dellist();
	 s_expo::setordmatrix(Mt); s_expo::setordtype(neworder);
	 for (Gnew.resetactual();Gnew.actualvalid();Gnew.moveactual()) {
		 Gnew.actualvalue().setdegree(wnew);
		 Gnew.actualvalue().reord(1); // Ordnet Monome bzgl. wnew und Mt
	 }
	 inwG.dellist();
	 G=std::move(Gnew);
	 G.sort();
	 // Interreduction of G:

	 if (comment) common_bb::WriteText("\nInter-reducing G!");
	 nred=s_polynom<T>::nreduction;
	 common_bb::interredBasis(G,1);
     if (comment) common_bb::WriteText("  Inter-reductions: " + ToString(s_polynom<T>::nreduction-nred));
 }
 else {
	 if (comment) common_bb::WriteText("\nInitial cone has not to be treated!");
	 s_expo::setordtype(neworder); s_expo::setordmatrix(Mt);
	 for (G.resetactual();G.actualvalid();G.moveactual()) {
		 G.actualvalue().setdegree(wnew);
		 G.actualvalue().reord(1);
	 }
	 inwG.dellist();
	 Gnew.dellist();
	 GnotMon.dellist();
	 G.sort();
 }


 while (!isdone) {
	 u=rational(0);
	 while (u!=rational(1)) {
         if (comment) common_bb::Clear();
		 getvectornextcone(G,u,wold,wnew,pwt);
		 if (u==rational(0)) break;
		 ncones++;
		 if (comment) common_bb::WriteText("\nNumber of cones: " + ToString(ncones));
		 nmon=Makelistsofinitials(G,wnew,inwG,Gnew,GnotMon,dmin);
		 if (comment) common_bb::WriteText("\nNumber of polynomials in inwG that are monomials: " + ToString(nmon));
		 crit = double(nmon)/double(m);
		 mold=m;

		 if (u==rational(1)) {
			 if ((isdone=isGroebner(G))) {
				 inwG.dellist();Gnew.dellist();GnotMon.dellist();
				 break;
			 }
			 if ((crit<0.6) && (l<n) && (!pwtpos)) {
				 l++;
				 deg++;
				 perturbweight(pwt,deg,Mt,l);
				 u=rational(0);
				 for (i=0;i<n;i++) wnew[i]=wold[i];
				 setpositiv(wnewpos,wspos,pwtpos,wnew,pwt);
				 inwG.dellist();Gnew.dellist();GnotMon.dellist();
				 continue;
			 }
		 }

		 // Criteria to choose a BB-algorithm
		 if ((crit < 0.8) && ispositv(wnew)) {
            for (pGlist=GnotMon;pGlist;pGlist++)  {
                pGlist().reord();
                pGlist().normalize();
                pGlist().setLTdegree(wnew);
            }
            for (pGlist=Gnew;pGlist;pGlist++) {pGlist().normalize();pGlist().setLTdegree(wnew);}
            Gnew.sort(nmon);
            GnotMon.sort(m-nmon);
            Gnew.merge(GnotMon);
            if (comment) common_bb::WriteText("\nComputing groebner basis with pos. hom. BB-algorithm!");
            m=bbhom::Groebnerpos(Gnew,wnew,dmin);
            anzpos++;
 		 }
		 else {
			 for (GnotMon.resetactual();GnotMon.actualvalid();GnotMon.moveactual()) GnotMon.actualvalue().reord();
			 GnotMon.resetactual();
			 if (comment) common_bb::WriteText("\nComputing groebner basis with prim. BB-algorithm!");
			 m=common_bb::primGroebner(Gnew,GnotMon);
			 anzprim++;
		 }

		 if (comment) common_bb::WriteText("\nNumber of elements in Gnew = " + ToString(m));

		 // ------------------------------------------------------------------------

		 // Basis Liften:
		 if (comment) common_bb::WriteText("\nMaximal integer-length: " + ToString(integer::GetMaxlength()));
		 remain = val::vector<char>(0,mold);
		 if (comment) common_bb::WriteText("\nLifting basis");

         for (Gnew.resetactual();Gnew.actualvalid();Gnew.moveactual()) {
            Gnew.actualvalue().setdegree(wold);
            Gnew.actualvalue().reord(1);
            Gnew.actualvalue().liftpolynom(inwG,G,remain);
         }
         for (Gnew.resetactual(),G.resetactual(),i=0;Gnew.actualvalid();Gnew.moveactual()) {
            if (Gnew.actualvalue().iszero()) {
                while (G.actualvalid() && !remain[i]) {
                    G.moveactual();i++;
                }
                Gnew.actualvalue()=std::move(G.actualvalue());
                G.moveactual();i++;
                if (!G.actualvalid()) break;
            }
        }
		 G.dellist();
        for (Gnew.resetactual();Gnew.actualvalid();Gnew.moveactual()) {
            Gnew.actualvalue().setdegree(wnew);
            Gnew.actualvalue().reord(1); // Sort monomials wrt. wnew and Mt
        }
        inwG.dellist();
        G=std::move(Gnew);
        G.sort();

        nred=s_polynom<T>::nreduction;
        if (comment) common_bb::WriteText("\nInter-reducing G!");
        common_bb::interredBasis(G,1);
	 }
	 if (l==1 || isdone) break;
	 else {
		 if (comment) common_bb::WriteText("\nGroebner walk arrived at pwt! Check if ready.");
		 isdone=isGroebner(G);
		 //saverrorlist(G);
		 //exit(-1);
	 }
	 if (!isdone) {
		 if (comment) common_bb::WriteText("\nGroebner walk not ready, walk next to wt!\n");
		 deg++;
		 perturbweight(pwt,deg,Mt,l);
		 setpositiv(wnewpos,wspos,pwtpos,wnew,pwt);
	 }
	 else {
		 if (comment) common_bb::WriteText("\nGroebner walk ready. Ordering basis!");
		 for (G.resetactual();G.actualvalid();G.moveactual()) G.actualvalue().setLTzerodegree();
		 G.sort();
	 }
 }

 for (G.resetactual();G.actualvalid();G.moveactual()) {
			 G.actualvalue().setzerodegree();
			 G.actualvalue().reord(); // Sort monomials wrt. wnew and Mt
 }
 G.sort(m);

 if (comment) {
	common_bb::WriteText("\n\nTime: " + ToString(Chrono()));
	common_bb::WriteText("\nGroebner walk successful! Number of cones: " + ToString(ncones));
	common_bb::WriteText("\nFinal perturbations-degree for target-vector: "+ ToString(l));
	common_bb::WriteText("\nprimitive BB-Alg. : "+ToString(anzprim));
	common_bb::WriteText("\nhom. pos. BB-Alg. : " + ToString(anzpos));
	common_bb::WriteText("\n\n Elements in G: "+ ToString(m));
	common_bb::WriteText("\n Monomials: " + ToString(s_polynom<T>::getmnumber()));
 }
 return;
}


template <class T>
void walkmain(const std::string& name,val::Glist<val::s_polynom<T> > &G,int neworder,int k,int l,const val::matrix<int> &M)
{
    if (!G.isempty()) G.dellist();
    common_bb::readfromfile<T>(name,G);
    walkmain(G,neworder,k,l,M);
}


} // end namespace


#endif // WALK_H_INCLUDED
