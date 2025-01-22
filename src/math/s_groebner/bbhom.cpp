
#include <s_groebner/bbhom.h>
#include <s_groebner/bb_int.h>

namespace bbhom
{

// Terms have an additional variable to homogenize them eventually.
template <>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::integer> > &P,int &wastotdegcompatible,int &washomogen)
{
	int i=0,n,nG=0,d,order,tdegisinit=1;
	val::s_polynom<val::integer> f;
	std::string line;
	wastotdegcompatible=1;
	washomogen=1;
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
	val::s_expo::setordtype(order);
	if (order==-1000) {
		val::matrix<int> Mold(n);
		int j;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				file>>Mold(i,j);
			}
		}
        for (i=0;i<n;i++)
             if (Mold(0,i)!=1) {
                wastotdegcompatible=0;
            }
		val::s_expo::setordmatrix(Mold);
	}
	else if (order==-1) {
        wastotdegcompatible=0;
	}
	else if (order!=0 && order!=-2) val::s_expo::setordtype(-2);

	val::s_expo::setdim(n+1);
	val::s_expo X(0);
	val::integer coeff,zero;

	do {
		file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				file>>X[i];
			}
			f.insert(std::move(coeff),X);
			d=X.totdeg();
		}
		else break;
		do {
			file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					file>>X[i];
				}
				if (X.totdeg()!=d) washomogen=0;
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
            if (!washomogen && !wastotdegcompatible) {
                if (tdegisinit) tdegisinit=totdegisinitialform(f);
            }
			f.normalize();
			P.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (file);
	file.close();
	if (tdegisinit) wastotdegcompatible=1;
	common_bb::WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.");
	return nG;
}


template <>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::modq> > &P,int &wastotdegcompatible,int &washomogen)
{
	int i=0,n,nG=0,d,order,tdegisinit=1;
	val::s_polynom<val::modq> f;
	std::string line;

	wastotdegcompatible=1;
	washomogen=1;
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

	if (i!=3) {
        common_bb::MyMessage("\nWrong type of file!");
        return 0;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	file>>val::modq::q>>n>>order;
	val::s_expo::setordtype(order);
	if (order==-1000) {
		val::matrix<int> Mold(n);
		int j;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				file>>Mold(i,j);
			}
		}

        for (i=0;i<n;i++)
             if (Mold(0,i)!=1) {
                wastotdegcompatible=0;
            }
		val::s_expo::setordmatrix(Mold);
	}
	else if (order==-1) {
        wastotdegcompatible=0;
	}
	else if (order!=0 && order!=-2) val::s_expo::setordtype(-2);

	val::s_expo::setdim(n+1);
	val::s_expo X(0);
	val::modq coeff,zero(0);

	do {
		file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				file>>X[i];
			}
			f.insert(coeff,X);
			d=X.totdeg();
		}
		else break;
		do {
			file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					file>>X[i];
				}
				if (X.totdeg()!=d) washomogen=0;
				f.insert(coeff,X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
            if (!washomogen && !wastotdegcompatible) {
                if (tdegisinit) tdegisinit=totdegisinitialform(f);
            }
			f.normalize();
			P.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (file);
	file.close();
	if (tdegisinit) wastotdegcompatible=1;
	common_bb::WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.");
	return nG;
}


int gauss(val::matrix<val::rational>& a,int m)
{
 int n,i,j,s,k=0,r=0,rows=rnumber(a);
 val::rational h,zero(0);

 if (m==-1) m=rnumber(a);
 if (m>=rnumber(a)) m=rnumber(a)-1;
 n=cnumber(a);

 if (m<=0 || n<=0) return 0;

 for (s=0;s<n;s++) {                     // Create left triangle matrix
     h=zero;
     for (i=r;i<=m;i++) {   				 // Get pivot
		 if ((h=a[i][s])!=zero) {
			 k=i;
			 break;
		 }
	 }
	 if (h!=zero) {
		 r++;
		 if (k!=(r-1))
			 a.swaprows(r-1,k);
		 h=a[r-1][s];
		 for (j=s;j<n;j++) a[r-1][j]/=h;   // divide row by pivot
		 for (i=r;i<rows;i++)
			 for (j=n-1;j>=s;j--) {
				 a[i][j]-=a[i][s]*a[r-1][j];
			 }
     }
 }
 return r;
}



void perturbmatrix(val::matrix<int> &M)
{
 int i,j,n=M.numberofcolumns(),r,m;
 val::matrix<val::rational> A(n+1,n);


 for (i=0;i<n;i++) {
     A[0][i]=val::rational(1);
 }

 for (i=1;i<=n;i++)
	 for (j=0;j<n;j++) A[i][j]=val::rational(M[i-1][j]);

 for (m=1;m<n;m++) {
	 r=gauss(A,m);
	 if (r!=m+1) break;
 }
 m--; // delete m-th row
 for (i=m;i>0;i--)
	 for (j=0;j<n;j++) M[i][j]=M[i-1][j];
 // Replace 0-th row with (1,..,1):
 for (i=0;i<n;i++) M[0][i]=1;
 return;
}


template <>
void SetKP(val::s_polynom<val::integer> &f,common_bb::KritPairs &KP)
{
    KP.if1=KP.if2=NULL;KP.ih=&f;KP.mf1=KP.mf2=NULL;KP.mh=NULL;
}


template <>
void SetKP(val::s_polynom<val::modq> &f,common_bb::KritPairs &KP)
{
    KP.if1=KP.if2=NULL;KP.mh=&f;KP.mf1=KP.mf2=NULL;KP.ih=NULL;
}


template <>
void SetKP<val::integer>(common_bb::KritPairs &KP,common_bb::spair &SP)
{
    KP.if1=SP.fint;KP.if2=SP.gint;KP.ih=NULL;KP.mh=NULL;
    KP.mf1=KP.mf2=NULL;
}


template <>
void SetKP<val::modq>(common_bb::KritPairs &KP,common_bb::spair &SP)
{
    KP.mf1=SP.fmodq;KP.mf2=SP.gmodq;KP.ih=NULL;KP.mh=NULL;
    KP.if1=KP.if2=NULL;
}

//

void updateG(val::Glist< val::s_polynom<val::integer> > &G,int &nG,val::Glist< val::s_polynom<val::integer> > &Gd,
             val::Glist<val::s_polynom<val::modq> > &Gp,val::Glist<val::s_polynom<val::modq> > &Gpd,val::Glist<common_bb::spair> &lspair,
             val::pol<val::integer>& hilbnumG)
{
 int i,order=val::s_expo::getordtype();
 val::vector<val::s_expo> I,I1;


 Gpd.resetactual();
 Gd.resetactual();
 while (!Gd.isempty()) {

        I1=val::vector<val::s_expo>(nG);
        for (i=0;i<nG;i++) hilbert::kgVdiv(G[i].LT(),Gd.actualvalue().LT(),I1[i]);
        val::s_expo::setordtype(-1);
        I1.sort();
        I=hilbert::reduce(I1);
        hilbnumG=hilbnumG-val::pol<val::integer>(val::integer(1),Gd.actualvalue().LT().totdeg()) * hilbert::Hilbertnum(I,I.dimension());
        val::s_expo::setordtype(order);

        common_bb::update(Gd.actualvalue(),Gpd.actualvalue(),G,Gp,lspair,nG);
        Gd.moveheadtoend(G);
        Gpd.moveheadtoend(Gp);
        nG++;
 }
}


void updateG(val::Glist< val::s_polynom<val::modq> > &G,int &nG,val::Glist< val::s_polynom<val::modq> > &Gd,
             val::Glist<common_bb::spair> &lspair,val::pol<val::integer>& hilbnumG)
{
 int i,order=val::s_expo::getordtype();
 val::vector<val::s_expo> I,I1;


 Gd.resetactual();
 while (!Gd.isempty()) {

        I1=val::vector<val::s_expo>(nG);
        for (i=0;i<nG;i++) hilbert::kgVdiv(G[i].LT(),Gd.actualvalue().LT(),I1[i]);
        val::s_expo::setordtype(-1);
        I1.sort();
        I=hilbert::reduce(I1);
        hilbnumG=hilbnumG-val::pol<val::integer>(val::integer(1),Gd.actualvalue().LT().totdeg()) * hilbert::Hilbertnum(I,I.dimension());
        val::s_expo::setordtype(order);

        common_bb::update(Gd.actualvalue(),G,lspair,nG);
        Gd.moveheadtoend(G);
        nG++;
 }
}


// M is totdeg-compatible order-matrix
template <>
int GroebnerwithHilbert(val::Glist<val::s_polynom<val::integer> > &G,int order,const val::matrix<int> &M,int nH)
{
 int nG=0,d,i,c,nThreads=common_bb::ComputingThreads,intThreads=0; //nd
 val::Glist<int> primlist;
 val::Glist<val::s_polynom<val::integer> > Gd,H,Hd;
 val::Glist<val::s_polynom<val::modq> > Gp,Gpd;
 val::s_polynom<val::integer> hint,fint,gint;
 val::s_polynom<val::modq> hmodq;
 val::Glist<common_bb::spair> Pair;
 val::GlistManipulator<val::s_polynom<val::integer> > pH;

 val::Glist<common_bb::KritPairs> ListKPairs;
 common_bb::KritPairs KP;
 val::vector<std::thread*> Thr(nThreads);

 val::pol<val::integer> hilbnumH,hilbnumG;
 val::s_expo t,t1;


 if (G.isempty()) return nG;

 if (nH==-1) nH=G.length();
 H=std::move(G);
 val::vector<val::s_expo> I(nH);

 for (pH=H,i=0;pH;pH++,i++) I[i]=pH().LT();
 val::s_expo::setordtype(-1);
 I.sort();
 hilbnumH=hilbert::Hilbertnum(I,nH);

 // Resort H:
 val::s_expo::setordmatrix(M);
 val::s_expo::setordtype(order);

 for (pH=H;pH;pH++) pH().reord();

 H.sort();

 common_bb::CreatePrimelist(primlist,H,nH);
 H.resetactual();
 d=H.actualvalue().LT().totdeg();

 //nd=0;
 while (!H.isempty() && H.actualvalue().LT().totdeg()==d) {
     H.actualvalue().reduction(Gd);
	 if (!H.actualvalue().iszero()) {
		 Gd.sinsert(std::move(H.actualvalue()));
		 //nd++;
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

 hilbnumG=val::pol<val::integer>(val::integer(1),0);

 updateG(G,nG,Gd,Gp,Gpd,Pair,hilbnumG);

 common_bb::WriteText("\nG initialized: "+ val::ToString(nG));
 common_bb::WriteText("\nActual degree: "+ val::ToString(d));

 H.resetactual();
 c=1;

 while (H.actualvalid() || c) {

	 common_bb::Clear();
	 common_bb::WriteText("\nnG = " + val::ToString(nG));
	 hilbert::getlowestmon(hilbnumG,hilbnumH,c,d);
	 common_bb::WriteText("\nNext degree: " + val::ToString(d));
	 common_bb::WriteText("\nElements of this degree: " + val::ToString(c));
	 if (c==0) break;

	 if (c<0) {
		 val::Error::error("\nComputation failed!!");
	 }

	 //nd=0;
	 while (H.actualvalid() && H.actualvalue().LT().totdeg() < d) {
        H.moveactual();
        nH--;
	 }

	 Pair.resetactual();
	 while (!Pair.isempty() && Pair.actualvalue().deg<d) {
        Pair.skiphead();
	 }

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
        Thr(i) = new std::thread(common_bb::reduceKritPairs<val::modq>,std::ref(ListKPairs),std::ref(Gp),std::ref(Gpd),c,0);
    for (i=0;i<nThreads;i++) Thr(i)->join();

    for (i=0;i<nThreads;i++) delete Thr(i);

    if (common_bb::KritPairs::nGd!=c) {
        val::Error::error("\n nGd != n . Computation failed!");
    }
    for (;ListKPairs.actualvalid();ListKPairs.moveactual()) ListKPairs.actualvalue().s_done=1;

    ListKPairs.resetactual();

    common_bb::KritPairs::p_teilt=0;
    common_bb::KritPairs::nGd=0;
    common_bb::Hdisbusy=false;

    intThreads=val::Min(c,nThreads);

    for (i=0;i<intThreads;i++)
        Thr(i) = new std::thread(common_bb::reduceKritPairs2,std::ref(ListKPairs),std::ref(G),std::ref(Gd),std::ref(Hd),c);
        //Thr(i) = new std::thread(common_bb::reduceKritPairs<val::integer>,std::ref(ListKPairs),std::ref(G),std::ref(Gd),c,1);
    for (i=0;i<intThreads;i++) Thr(i)->join();

    for (i=0;i<intThreads;i++) delete Thr(i);

    ListKPairs.dellist();

    if (common_bb::KritPairs::nGd!=c) {
        //exit(-1);
        val::Error::error("\n nGd != n . Computation failed!");
    }

     if (common_bb::KritPairs::p_teilt) {
         val::modq::q=common_bb::getnextprime(G,primlist,val::integer(1),Gd);
         common_bb::restoreGp(Gp,G);
     }

     common_bb::interredBasis(Gd);
	 common_bb::restoreGp(Gpd,Gd);
	 updateG(G,nG,Gd,Gp,Gpd,Pair,hilbnumG);
 }

 return nG;
}


// M is totdeg-compatible order-matrix
template <>
int GroebnerwithHilbert(val::Glist<val::s_polynom<val::modq> > &G,int order,const val::matrix<int> &M,int nH)
{
 int nG=0,d,i,c,nThreads=common_bb::ComputingThreads; //nd
 val::Glist<val::s_polynom<val::modq> > Gd,H;
 val::Glist<common_bb::spair> Pair;
 val::GlistManipulator<val::s_polynom<val::modq> > pH;

 val::Glist<common_bb::KritPairs> ListKPairs;
 common_bb::KritPairs KP;
 val::vector<std::thread*> Thr(nThreads);

 val::pol<val::integer> hilbnumH,hilbnumG;
 val::s_expo t,t1;

 if (G.isempty()) return nG;

 if (nH==-1) nH=G.length();
 H=std::move(G);
 val::vector<val::s_expo> I(nH);

 for (pH=H,i=0;pH;pH++,i++) I[i]=pH().LT();
 val::s_expo::setordtype(-1);
 I.sort();
 hilbnumH=hilbert::Hilbertnum(I,nH);

 // Resort H:
 val::s_expo::setordmatrix(M);
 val::s_expo::setordtype(order);

 for (pH=H;pH;pH++) pH().reord();

 H.sort();

 H.resetactual();
 d=H.actualvalue().LT().totdeg();

 //nd=0;
 while (!H.isempty() && H.actualvalue().LT().totdeg()==d) {
     H.actualvalue().reduction(Gd);
	 if (!H.actualvalue().iszero()) {
		 Gd.sinsert(std::move(H.actualvalue()));
		 //nd++;
	 }
	 nH--;
	 H.skiphead();
 }

 common_bb::interredBasis(Gd);

 hilbnumG=val::pol<val::integer>(val::integer(1),0);

 updateG(G,nG,Gd,Pair,hilbnumG);

 common_bb::WriteText("\nG initialized: "+ val::ToString(nG));
 common_bb::WriteText("\nActual degree: "+ val::ToString(d));

 H.resetactual();
 c=1;

 while (H.actualvalid() || c) {

	 common_bb::Clear();
	 common_bb::WriteText("\nnG = " + val::ToString(nG));
	 hilbert::getlowestmon(hilbnumG,hilbnumH,c,d);
	 common_bb::WriteText("\nNext degree: " + val::ToString(d));
	 common_bb::WriteText("\nElements of this degree: " + val::ToString(c));
	 if (c==0) break;

	 if (c<0) {
		 val::Error::error("\nComputation failed!");
	 }

	 //nd=0;
	 while (H.actualvalid() && H.actualvalue().LT().totdeg() < d) {
        H.moveactual();
        nH--;
	 }

	 Pair.resetactual();
	 while (!Pair.isempty() && Pair.actualvalue().deg<d) {
        Pair.skiphead();
	 }

    common_bb::KritPairs::nGd=0;
    KP.s_done=0;
    for (;H.actualvalid() && H.actualvalue().LT().totdeg()==d;H.moveactual()) {
        SetKP(H.actualvalue(),KP);
        ListKPairs.inserttoend(KP);
    }
    while(!Pair.isempty() && Pair.actualvalue().deg==d) {
        if (Pair.actualvalue().s_done) {
            Pair.skiphead();continue;
        }
        SetKP<val::modq>(KP,Pair.actualvalue());
        ListKPairs.inserttoend(KP);
        Pair.skiphead();
    }

    ListKPairs.resetactual();

    for (i=0;i<nThreads;i++)
        Thr(i) = new std::thread(common_bb::reduceKritPairs<val::modq>,std::ref(ListKPairs),std::ref(G),std::ref(Gd),c,0);
    for (i=0;i<nThreads;i++) Thr(i)->join();

    for (i=0;i<nThreads;i++) delete Thr(i);

    if (common_bb::KritPairs::nGd!=c) {
        val::Error::error("\n nGd != n . Computation failed!");
    }

    ListKPairs.dellist();

    common_bb::interredBasis(Gd);
	updateG(G,nG,Gd,Pair,hilbnumG);
 }

 return nG;
}


//
template <>
int bbhom(val::Glist<val::s_polynom<val::modq> > &G,int washomogen,int wastotdegcompatible)
{
    int perturbed=0,order=val::s_expo::getordtype(),nG=0;
    val::matrix<int> M;
    val::ChronoClass Chrono;
    std::string s;

    if (G.isempty()) return 0;

    val::s_polynom<val::modq>::nreduction=0;

    if (washomogen && wastotdegcompatible) {
        val::s_expo::setdim(val::s_expo::getdim()-1);
        Chrono();
        nG=Homgroebner(G);
    }
    else if (washomogen && !wastotdegcompatible) {
        val::s_expo::setdim(val::s_expo::getdim()-1);
        if (order==-1) {
            val::s_expo::setordtype(0);
            perturbed=1;
        }
        else if (order==-1000) {
            perturbed=1;
            M=val::s_expo::getordmatrix();
            perturbmatrix(val::s_expo::getordmatrix());
        }
        G.sort();
        Chrono();
        nG=Homgroebner(G);
    }
    else { // Homogenize:
        int n=val::s_expo::getdim();
        val::matrix<int> Mnew;
        if (order==-1000) {
            M=std::move(val::s_expo::getordmatrix());
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
            val::s_expo::setordtype(-1000);
            val::s_expo::setordmatrix(Mnew);
        }
        else if (order==-1) {val::s_expo::setordtype(0);}
        else if (order==0) {
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
        else {val::s_expo::setordtype(-2);order=-2;}
        for (G.resetactual();G.actualvalid();G.moveactual()) G.actualvalue().homogenize();
        G.sort();
        Chrono();
        nG=Homgroebner(G);
        common_bb::WriteText("\nReductions: " + val::ToString(val::s_polynom<val::modq>::nreduction));
        common_bb::WriteText("\nDehomogenize and reduce...");
        nG=common_bb::dehomred(G,nG,order,M);
    }

    s+="\n\nTime: " + val::ToString(Chrono());
    s+="\nElements in G: " + val::ToString(nG);
    s+="\nMonomials: " + val::ToString(val::s_polynom<val::modq>::getmnumber());
    s+="\nReductions: " + val::ToString(val::s_polynom<val::modq>::nreduction);
    common_bb::WriteText(s);

    if (perturbed) {
        val::s_expo::getordmatrix()=std::move(M);
        val::s_expo::setordtype(order);
        G.sort();
    }

    val::s_polynom<val::modq>::nreduction=0;
    return nG;
}


template <>
int bbhom(val::Glist<val::s_polynom<val::integer> > &G,int washomogen,int wastotdegcompatible)
{
    int perturbed=0,order=val::s_expo::getordtype(),nG=0;
    val::matrix<int> M;
    val::ChronoClass Chrono;
    std::string s;

    if (G.isempty()) return 0;

    val::s_polynom<val::modq>::nreduction=0;

    if (washomogen && wastotdegcompatible) {
        val::s_expo::setdim(val::s_expo::getdim()-1);
        Chrono();
        nG=Homgroebner(G);
    }
    else if (washomogen && !wastotdegcompatible) {
        val::s_expo::setdim(val::s_expo::getdim()-1);
        if (order==-1) {
            val::s_expo::setordtype(0);
            perturbed=1;
        }
        else if (order==-1000) {
            perturbed=1;
            M=val::s_expo::getordmatrix();
            perturbmatrix(val::s_expo::getordmatrix());
        }
        G.sort();
        Chrono();
        nG=Homgroebner(G);
    }
    else { // Homogenize:
        int n=val::s_expo::getdim();
        val::matrix<int> Mnew;
        if (order==-1000) {
            M=std::move(val::s_expo::getordmatrix());
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
            val::s_expo::setordtype(-1000);
            val::s_expo::setordmatrix(Mnew);
        }
        else if (order==-1) {val::s_expo::setordtype(0);}
        else if (order==0) {
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
        else {val::s_expo::setordtype(-2);order=-2;}
        for (G.resetactual();G.actualvalid();G.moveactual()) G.actualvalue().homogenize();
        G.sort();
        return bb_int::bbhommod(G,order,M);
    }


    s+="\n\nTime: " + val::ToString(Chrono());
    s+="\nElements in G: " + val::ToString(nG);
    s+="\nMonomials: " + val::ToString(val::s_polynom<val::integer>::getmnumber());
    s+="\nReductions: " + val::ToString(val::s_polynom<val::integer>::nreduction);
    common_bb::WriteText(s);

    if (perturbed) {
        val::s_expo::getordmatrix()=std::move(M);
        val::s_expo::setordtype(order);
        G.sort();
    }

    val::s_polynom<val::modq>::nreduction=0;
    return nG;
}



template <>
int bbhommain(const std::string &name,val::Glist<val::s_polynom<val::integer> >& G)
{

    int nG,wastotdegcompatible,washomogen;
    if (!G.isempty()) G.dellist();

    readfromfile(name,G,wastotdegcompatible,washomogen);

    nG=bbhom(G,washomogen,wastotdegcompatible);
    common_bb::WriteText("\nMaximal integer-length: " + val::ToString(val::integer::GetMaxlength()));
    common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));
    return nG;
}


template <>
int bbhommain(const std::string &name,val::Glist<val::s_polynom<val::modq> >& G)
{

    int nG,wastotdegcompatible,washomogen;
    if (!G.isempty()) G.dellist();

    readfromfile(name,G,wastotdegcompatible,washomogen);
    nG=bbhom(G,washomogen,wastotdegcompatible);
    return nG;
}


template <>
int hilbertconversionmain(const std::string &name,val::Glist<val::s_polynom<val::modq> > &G,int order,const val::matrix<int> &M)
{
    int nG,washomogen,wastotdegcompatible;
    if (!G.isempty()) G.dellist();

    readfromfile(name,G,wastotdegcompatible,washomogen);
    nG=hilbertconversion(G,washomogen,wastotdegcompatible,order,M);
    return nG;
}


template <>
int hilbertconversionmain(const std::string &name,val::Glist<val::s_polynom<val::integer> > &G,int order,const val::matrix<int> &M)
{
    int nG,washomogen,wastotdegcompatible;
    if (!G.isempty()) G.dellist();

    readfromfile(name,G,wastotdegcompatible,washomogen);
    nG=hilbertconversion(G,washomogen,wastotdegcompatible,order,M);
    common_bb::WriteText("\nMaximal integer-length: " + val::ToString(val::integer::GetMaxlength()));
    common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));
    return nG;
}


} // end namespace
