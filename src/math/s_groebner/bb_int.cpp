
#include <s_groebner/bb_int.h>
#include <pol_arithmetic.h>


namespace bb_int
{

int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::integer> > &G,int &order,val::matrix<int> &Mold,int onlytotdegcompatible) // Reads polynomials from file name
{
	int i=0,n,nG=0;
	val::s_polynom<val::rational> f;
	std::string line;
	val::matrix<int> Mnew;

	std::ifstream file(name,std::ios::in);
	if (!file) {
		common_bb::MyMessage("FILE DOES NOT EXIST");
        return 0;
	}

    do {
        getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);

	if (i!=2) {
        common_bb::MyMessage("\nWrong type of file!");
        return 0;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	file>>n>>order;
	if (order==-1000) {
		int j;
		Mold=val::matrix<int>(n);
		Mnew = val::matrix<int>(n+1);
		val::s_expo::setordtype(-1000);
		for (i=0;i<=n;i++) Mnew(0,i)=1;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				file>>Mold(i,j);
				Mnew(i+1,j)=Mold(i,j);
			}
			Mnew(i+1,n)=0;
		}
		if (onlytotdegcompatible) {
            for (i=0;i<n;i++)
                if (Mold(0,i)!=1) {
                   common_bb::MyMessage("\nOrder is not totdeg-compatible!");
                   return 0;
                }
		}
		val::s_expo::setordmatrix(Mnew);
		val::s_expo::setordtype(-1000);
	}

	else if (order==-2) {
		val::s_expo::setordtype(-2);
	}
	else if (order==-1) {
        if (onlytotdegcompatible) {common_bb::MyMessage("\nOrder is not totdeg-compatible!"); return 0;}
		val::s_expo::setordtype(0);
	}
	else if (order==0) {
        Mnew = val::matrix<int>(0,n+1,n+1);
        for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
        Mnew(0,n) =1;
        for (i=2;i<=n;i++) Mnew(i,i-2) =1;
        val::s_expo::setordtype(-1000);
        val::s_expo::setordmatrix(Mnew);
	}

	else val::s_expo::setordtype(-2);


	val::s_expo::setdim(n+1);
	val::s_expo X(0);
	val::rational coeff,zero;

	do {
		file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				file>>X[i];
			}
			f.insert(std::move(coeff),X);
		}
		else break;
		do {
			file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					file>>X[i];
				}
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			f.homogenize();
			f.reord();
			G.sinsert(val::primitivpart(f));
			f.del();
			nG++;
		}
		else break;
	}
	while (file);
	file.close();
	common_bb::WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.");
	return nG;
}


int homogenize(val::Glist<val::s_polynom<val::integer>> &G,int &order,val::matrix<int> &Mold)
{
    using namespace val;
    if (G.isempty()) return 0;

    int i,n=s_expo::getdim(),nG=0;
    Glist<s_polynom<integer>> H;
    s_polynom<integer> f;
    matrix<int> Mnew;

    order=s_expo::getordtype();

	if (order==-1000) {
		int j;
		Mold=s_expo::getordmatrix();
		Mnew = val::matrix<int>(n+1);
		val::s_expo::setordtype(-1000);
		for (i=0;i<=n;i++) Mnew(0,i)=1;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				Mnew(i+1,j)=Mold(i,j);
			}
			Mnew(i+1,n)=0;
		}
		val::s_expo::setordmatrix(Mnew);
		val::s_expo::setordtype(-1000);
	}
	else if (order==-2) {
		val::s_expo::setordtype(-2);
	}
	else if (order==-1) {
		val::s_expo::setordtype(0);
	}
	else if (order==0) {
        Mnew = val::matrix<int>(0,n+1,n+1);
        for (i=0;i<n;i++) Mnew(0,i) = Mnew(1,i) = 1;
        Mnew(0,n) =1;
        for (i=2;i<=n;i++) Mnew(i,i-2) =1;
        val::s_expo::setordtype(-1000);
        val::s_expo::setordmatrix(Mnew);
	}
	else val::s_expo::setordtype(-2);

	val::s_expo::setdim(n+1);
	val::s_expo X(0);

    for (const auto& g : G) {
        f.del();
        for (auto monom = g.begin();monom;++monom) {
            for (i=0;i<n;++i) X[i] = monom.actualterm()[i];
            f.insert(monom.actualcoef(),X);
        }
        f.homogenize();
        f.reord();
        H.push_back(std::move(f));
        nG++;
    }

    G=std::move(H);
    G.sort();
    return nG;
}


void updateG(val::Glist< val::s_polynom<val::integer> > &G,int &nG,val::Glist< val::s_polynom<val::integer> > &Gd,
             val::Glist<val::s_polynom<val::modq> > &Gp,val::Glist<val::s_polynom<val::modq> > &Gpd,val::Glist<common_bb::spair> &lspair)
{
 Gpd.resetactual();
 Gd.resetactual();
 while (!Gd.isempty()) {
        common_bb::update(Gd.actualvalue(),Gpd.actualvalue(),G,Gp,lspair,nG);
        Gd.moveheadtoend(G);
        Gpd.moveheadtoend(Gp);
        nG++;
 }
}


int HomGroebner(val::Glist<val::s_polynom<val::integer> > &H,val::Glist<val::s_polynom<val::integer> > &G,int nH)
{
 int nG=0,nd,dh,ds,d,nThreads=common_bb::ComputingThreads,i,intThreads;
 //int anzs=0,anzsd=0,anzh=0,anzhd=0;
 val::Glist<int> primlist;
 val::Glist<val::s_polynom<val::integer> > Gd;
 val::Glist<val::s_polynom<val::modq> > Gp,Gpd;
 val::s_polynom<val::integer> hint,fint,gint;
 val::s_polynom<val::modq> hmodq;
 val::Glist<common_bb::spair> Pair;
 val::GlistIterator<val::s_polynom<val::integer> > pH;
 val::Glist<common_bb::KritPairs> ListKPairs;
 common_bb::KritPairs KP;
 val::vector<std::thread*> Thr(nThreads);

 if (nH==0) return nG;

 // 1. Create suitable list of primes:
 common_bb::CreatePrimelist(primlist,H,nH);

 H.resetactual();
 // Get minimal degree:
 d=H.actualvalue().LT().totdeg();

 nd=0;

 while (!H.isempty() && H.actualvalue().LT().totdeg()==d) {
     H.actualvalue().reduction(Gd);
	 if (!H.actualvalue().iszero()) {
		 Gd.sinsert(std::move(H.actualvalue()));
		 nd++;
	 }
	 nH--;
	 H.skiphead();
 }

 common_bb::interredBasis(Gd);

 val::modq::q=common_bb::getnextprime(Gd,primlist,val::integer(1));

 for (pH=Gd;pH;pH++) {
	 hmodq.convertfrom(pH());
	 hmodq.normalize();
	 Gpd.inserttoend(std::move(hmodq));
 }
 updateG(G,nG,Gd,Gp,Gpd,Pair);

 H.resetactual();

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

	 //anzsd=0;
	 //anzhd=0;
	 nd=0;
//-------------------------------------------------------------------------
    common_bb::KritPairs::nGd=0;
    KP.s_done=0;
    for (;H.actualvalid() && H.actualvalue().LT().totdeg()==d;H.moveactual()) {
        KP.if1=KP.if2=NULL;KP.ih=&(H.actualvalue());KP.mf1=KP.mf2=NULL;KP.mh=NULL;
        ListKPairs.inserttoend(KP);
    }
    while(!Pair.isempty() && Pair.actualvalue().deg==d) {
        if (Pair.actualvalue().s_done) {
            Pair.skiphead();continue;
        }
        KP.if1=Pair.actualvalue().fint;KP.if2=Pair.actualvalue().gint;KP.ih=NULL;KP.mh=NULL;
        KP.mf1=Pair.actualvalue().fmodq;KP.mf2=Pair.actualvalue().gmodq;
        KP.s_done=Pair.actualvalue().s_done;
        ListKPairs.inserttoend(KP);
        Pair.skiphead();
    }

    ListKPairs.resetactual();

    for (i=0;i<nThreads;i++)
        Thr(i) = new std::thread(common_bb::reduceKritPairs<val::modq>,std::ref(ListKPairs),std::ref(Gp),std::ref(Gpd),0,0);
    for (i=0;i<nThreads;i++) Thr(i)->join();

    for (i=0;i<nThreads;i++) delete Thr(i);

    ListKPairs.resetactual();

    nd=common_bb::KritPairs::nGd;
    common_bb::KritPairs::p_teilt=0;
    common_bb::KritPairs::nGd=0;

    intThreads=val::Min(nThreads,nd);
    for (i=0;i<intThreads;i++)
        Thr(i) = new std::thread(common_bb::reduceKritPairs<val::integer>,std::ref(ListKPairs),std::ref(G),std::ref(Gd),0,1);
    for (i=0;i<intThreads;i++) Thr(i)->join();

    for (i=0;i<intThreads;i++) delete Thr(i);

    ListKPairs.dellist();
     //---------------------------------------------------------------------------
     if (common_bb::KritPairs::p_teilt) {
         val::modq::q=common_bb::getnextprime(G,primlist,val::integer(1),Gd);
         common_bb::restoreGp(Gp,G);
     }
	 common_bb::interredBasis(Gd);
	 common_bb::restoreGp(Gpd,Gd);
	 updateG(G,nG,Gd,Gp,Gpd,Pair);
	 //anzh+=anzhd;
	 //anzs+=anzsd;
 }
 Gp.dellist();
 primlist.dellist();
 Pair.dellist();

 return nG;
}


int bbhommod(val::Glist<val::s_polynom<val::integer> > &H,int order,val::matrix<int> &M,int nH)
{
 int nG;
 val::Glist<val::s_polynom<val::integer> > G;
 val::ChronoClass Chrono;

 if (nH==-1) nH=H.length();

 val::s_polynom<val::integer>::nreduction=0;
 val::s_polynom<val::modq>::nreduction=0;

 Chrono();
 nG=HomGroebner(H,G,nH);

 nG=common_bb::dehomred(G,nG,order,M);

 common_bb::WriteText("\nCheck if <H> = <G>...");
 if (common_bb::isinG(H,G,0)) common_bb::WriteText("\n <H> = <G> !");
 else {common_bb::MyMessage("\n<H> != <G> !");}

 common_bb::WriteText("\nCheck if G is groebner basis...");

 if (common_bb::isgroebner(G,nG,0)) common_bb::WriteText("\nG is groebner basis!");
 else common_bb::MyMessage("\nG is not groebner basis!");

 common_bb::WriteText("\n\nTime in sec.: " + val::ToString(Chrono()));
 common_bb::WriteText("\nG has " + val::ToString(nG) + " elements.");
 common_bb::WriteText("\nNumber of monomials: " + val::ToString(val::s_polynom<val::integer>::getmnumber()));
 common_bb::WriteText("\nReductions: " + val::ToString(val::s_polynom<val::integer>::nreduction));
 H=std::move(G);

 val::s_polynom<val::integer>::nreduction=0;
 val::s_polynom<val::modq>::nreduction=0;

 return nG;
}

int bbhommod(val::Glist<val::s_polynom<val::integer>> &H,int comment)
{
    using namespace val;
    if (H.isempty()) return 0;
    int nG,order;
    matrix<int> Mold;
    Glist<s_polynom<integer>> G;
    val::ChronoClass Chrono;


    nG=homogenize(H,order,Mold);
    Chrono();
    nG=HomGroebner(H,G,nG);
    nG=common_bb::dehomred(G,nG,order,Mold);
    if (comment) common_bb::WriteText("\nCheck if <H> = <G>...");
    if (!common_bb::isinG(H,G,0)) {
        if (comment) common_bb::WriteText("\n <H> != <G>");
        return 0;
    }
    else if (comment) common_bb::WriteText("\n <H> = <G> !");
    if (!common_bb::isgroebner(G,nG,0)) {
        if (comment) common_bb::WriteText("\nG is not groebner basis!");
        return 0;
    }
    else if (comment) common_bb::WriteText("\nG is groebner basis!");

    if (comment) {
        common_bb::WriteText("\n\nTime in sec.: " + val::ToString(Chrono()));
        common_bb::WriteText("\nG has " + val::ToString(nG) + " elements.");
        common_bb::WriteText("\nNumber of monomials: " + val::ToString(val::s_polynom<val::integer>::getmnumber()));
        common_bb::WriteText("\nReductions: " + val::ToString(val::s_polynom<val::integer>::nreduction));
        common_bb::WriteText("\nMaximal integer-length: " + val::ToString(val::integer::GetMaxlength()));
        common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));
    }

    H=std::move(G);

    return nG;
}


int bbhommodmain(const std::string &argv,val::Glist<val::s_polynom<val::integer> > &G)
{
 int nG,order;
 val::matrix<int> M;

 nG=readfromfile(argv,G,order,M);

 nG=bbhommod(G,order,M,nG);

 common_bb::WriteText("\nMaximal integer-length: " + val::ToString(val::integer::GetMaxlength()));
 common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));

 return nG;
}


}  // end namespace
