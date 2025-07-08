
#include <thread>
#include <s_groebner/common_bb.h>
#include <fstream>
#include <rand.h>
#include <pol_arithmetic.h>
#include <numbers.h>
#include <MyTime.h>


namespace common_bb
{
//

int MaxThreads=std::thread::hardware_concurrency(), ComputingThreads=MaxThreads;
std::string primlistpath(val::primlistpath);

std::mutex ProtectData,ProtectOutput;

//
std::atomic<int> KritPairs::nGd(0);
std::atomic<int> KritPairs::p_teilt(0);
std::atomic<bool> Hdisbusy(false);
std::mutex KritPairs::protectGd;
std::mutex protectHd;

//
void SetComputingThreads(int n)
{
    if (n<0 || n>MaxThreads) ComputingThreads=MaxThreads;
    else ComputingThreads=n;
}

//

void Message_to_cout(const std::string& s)
{
    std::cout<<s;
}

void (*messageoutput) (const std::string&) = &Message_to_cout;

void MyMessage(const std::string& s)
{
    messageoutput(s);
}


void Clear_console()
{
    std::cout<<"\r";
    return;
}

void (*clearoutput) (void) = &Clear_console;


void Clear()
{
    clearoutput();
}


void WriteText_to_cout(const std::string& s)
{
    std::cout<<s<<std::flush;
}

void (*writetextoutput) (const std::string&) = &WriteText_to_cout;

void WriteText(const std::string& s)
{
    writetextoutput(s);
}


int lcmdis(const val::s_expo &t1,const val::s_expo &t2,val::s_expo &t)
{
 int i,is=1,dim=val::s_expo::getdim();

 for (i=0;i<dim;i++) {
	 t[i]=val::Max(t1[i],t2[i]);
	 if (t1[i] && t2[i]) is=0;
 }
 return is;
}

//

int gettype (const std::string &name)
{
    int i=0;
    std::string line;

	std::ifstream file(name,std::ios::in);
	if (!file) {
		MyMessage("\nCannot read file!!!!");
        return 0;
	}
	do {
        getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);
	file.close();
	return i;
}


int gettype(char *name)
{
    return gettype(std::string(name));
}

int getnumberofvariables(const std::string& name)
{
    int m=gettype(name),n=0;

	if (m!=2 && m!=3) {
        MyMessage("\nWRONG TYPE OF FILE!");
        return 0;
	}
	std::ifstream file(name,std::ios::in);
	if (!file) {
		MyMessage("\nCannot read file!!!!");
        return 0;
	}
	file>>n;
	if (m==3) file>>n;
	return n;
}


template <>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::modq> >& P,int onlytotdeghomogen)
{
	int n=0,i=0,j,d=0,dim=0,ordtype;
	val::s_polynom<val::modq> f;
	val::matrix<int> M;
	std::string line;

	WriteText("\nReading list...");

	std::ifstream file(name,std::ios::in);
	if (!file) {
        MyMessage("FILE DOES NOT EXIST");
        return 0;
	}

    do {
        getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);

	if (i!=3) {
        MyMessage("\nWrong type of file!");
        return 0;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	file>>val::modq::q;file>>dim>>ordtype;
	val::s_expo::setordtype(ordtype);
	val::s_expo::setdim(dim);

	if (ordtype==-1000) {
		M=val::matrix<int>(dim);
		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++) file>>M(i,j);
        val::s_expo::setordmatrix(std::move(M));
	}
	else if (ordtype!=-2 && ordtype!= -1 && ordtype !=0) val::s_expo::setordtype(-2);


	val::modq coef;
	val::s_expo X;

	do {
		file>>coef;
		if (coef!=0) {
            file>>X;
            d=X.totdeg();
            f.insert(std::move(coef),X);
		}
		else break;
		do {
			file>>coef;
			if (coef!=0){
			    file>>X;
			    if (onlytotdeghomogen && X.totdeg()!=d) {
                    MyMessage("Basis is not totdeg-homogeneous!");
                    file.close();
                    P.dellist();
                    return 0;
			    }
			    f.insert(std::move(coef),X);
			}
		}
		while (coef!=0);

		if (!f.iszero()){
			f.normalize();
			P.sinsert(std::move(f));
			n++;
		}
		else break;
	}
	while (file);
	file.close();
	WriteText("\nList red!");
	WriteText("\nMonomials: " + val::ToString(val::s_polynom<val::modq>::getmnumber()));
	return n;
}

template <>
int readfromfile(char *name,val::Glist<val::s_polynom<val::modq> >& P,int onlytotdeghomogen)
{
    return readfromfile(std::string(name),P,onlytotdeghomogen);
}


template <>
int readfromfile(const std::string &name,val::Glist<val::s_polynom<val::integer> >& P,int onlytotdeghomogen)
{
	int i=0,n=0,j,d=0,dim,ordtype;
	val::rational zero;
	val::s_polynom<val::rational> f;
	val::matrix<int> M;
	std::string line;

	std::ifstream file(name,std::ios::in);
	if (!file) {
		MyMessage("FILE DOES NOT EXIST");
        return 0;
	}

    do {
        getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);

	if (i!=2) {
        MyMessage("\nWrong type of file!");
        return 0;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	file>>dim>>ordtype;

	val::s_expo::setordtype(ordtype);
	val::s_expo::setdim(dim);

	if (ordtype==-1000) {
		M=val::matrix<int>(dim);
		for (i=0;i<dim;i++)
			for (j=0;j<dim;j++) file>>M(i,j);
		val::s_expo::setordmatrix(std::move(M));
	}
	else if (ordtype!=-2 && ordtype!=-1 && ordtype!=0) val::s_expo::setordtype(-2);

	val::rational coef;
	val::s_expo X;

	do {
		file>>coef;
		if (coef!=zero) {
            file>>X;
			d=X.totdeg();
			f.insert(std::move(coef),X);
        }
		else break;
		do {
			file>>coef;
			if (coef!=zero) {
			    file>>X;
			    if (onlytotdeghomogen && X.totdeg()!=d) {
                    MyMessage("Basis is not totdeg-homogeneous!");
                    P.dellist();
                    file.close();
                    return 0;
			    }
			    f.insert(std::move(coef),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			P.sinsert(val::primitivpart(f));
			f.del();
			n++;
		}
	}
	while (file);
	file.close();
	return n;
}


template <>
int readfromfile(char *name,val::Glist<val::s_polynom<val::integer> >& P,int onlytotdeghomogen)
{
    return readfromfile(std::string(name),P,onlytotdeghomogen);
}


template <>
void savelist(const std::string &arg,const val::Glist<val::s_polynom<val::modq> >& P)
{
 int i,j,dim=val::s_expo::getdim();
 std::ofstream file(arg,std::ios::out | std::ios::trunc);
 if (!file) {
	 MyMessage("\nCANNOT WRITE IN FILE!!!");
	 return;
 }

 val::GlistIterator<val::s_polynom<val::modq> > ItP;

 file<<val::modq::q<<std::endl;
 file<<dim<<std::endl;
 file<<val::s_expo::getordtype()<<std::endl;

 if (val::s_expo::getordtype()==-1000) {
    file<<std::endl;
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) file<<val::s_expo::getordmatrix()(i,j)<<" ";
        file<<std::endl;
	 }
 }
 file<<std::endl;

 for (ItP.settohead(P);ItP.actualvalid();ItP.moveactual()) {
        file<<ItP.getelement()<<std::endl;
 }

 file<<'0'<<std::endl;
 file.close();
}


template <>
void savelist(const std::string &arg,const val::Glist<val::s_polynom<val::integer> >& P)
{
 int i,j,dim=val::s_expo::getdim();
 std::ofstream file(arg,std::ios::out | std::ios::trunc);
 if (!file) {
	 MyMessage("\nCANNOT WRITE IN FILE!!!");
	 return;
 }

 val::GlistIterator<val::s_polynom<val::integer> > ItP;

 file<<dim<<std::endl;
 file<<val::s_expo::getordtype()<<std::endl;

 if (val::s_expo::getordtype()==-1000) {
    file<<std::endl;
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) file<<val::s_expo::getordmatrix()(i,j)<<" ";
        file<<std::endl;
	 }
 }
 file<<std::endl;

 for (ItP.settohead(P);ItP.actualvalid();ItP.moveactual()) {
        file<<ItP.getelement()<<std::endl;
 }

 file<<'0'<<std::endl;
 file.close();
}




void CreatePrimelist(val::Glist<int> &Primelist,const val::Glist< val::s_polynom<val::integer> > &F,int r)
{
 int i,j,p,anz=0;
 val::integer zero;

 if (r==-1) r=F.length();
 std::ifstream file(primlistpath,std::ios::in);

 if (!file) {
	int pp=100000;
	for (int i=0 ;i<1000;++i) {
		pp=val::nextprime(pp);
		Primelist.push_back(pp);
	}
	return;
 }
 file>>j>>p;
 while (file) {
	 file>>j>>p;
	 for (i=0;i<r;i++)
		 if (F[i].LC() %val::integer(p)==zero) break;
	 if (i==r) {
		 Primelist.inserttoend(p);
		 anz++;
	 }
 }
 file.close();
}


int getnextprime(const val::Glist< val::s_polynom<val::integer> > &P,val::Glist<int> &Primlist,
                 const val::integer &wert,const val::Glist<val::s_polynom<val::integer> > &Gd)
{
 int p;
 val::integer zero(0);
 val::GlistIterator<val::s_polynom<val::integer> > ptoP;

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


//
int Maximalinteger(const val::Glist< val::s_polynom<val::integer> >& P)
{
 int maxint=0;
 val::GlistIterator<val::s_polynom<val::integer> > pG;
 val::s_polynomIterator<val::integer> p;

 for (pG.settohead(P);pG.actualvalid();pG.moveactual()) {
	 for (p=pG.getelement();p;p++) {
		 if (p.actualcoef().abslength()>maxint) maxint=p.actualcoef().abslength();
	 }
 }

 return maxint;
}

//
int isinG(const std::string& hfile,const std::string& gbfile)
{
    using namespace val;
    int i,is=0,type=gettype(gbfile),ngb,nfile,dummy,nH=0;

    if (type!=2 && type!=3) {
        MyMessage("\nWRONG TYPE OF FILE");
        return 0;
    }
    if (gettype(hfile)!=type) {
        WriteText("\nBases have different types!");
        return 0;
    }
    ngb=getnumberofvariables(gbfile);
    nfile=getnumberofvariables(hfile);

    if (nfile >ngb ) {
        WriteText("\nH has more variables than G!");
        return 0;
    }

    Glist<s_polynom<modq> > Gmodq,Hmodq;
    s_polynom<modq> fmodq;
    modq modqcoef,modqzero(0);

    Glist<s_polynom<integer> > Gint,Hint;
    s_polynom<integer> fint;
    integer intcoef,intzero(0);

    val::ChronoClass Chrono;

    s_polynom<modq>::nreduction=0;
    s_polynom<integer>::nreduction=0;

    if (type==3) readfromfile(gbfile,Gmodq);
    else readfromfile(gbfile,Gint);

    std::ifstream file(hfile,std::ios::in);

    if (type==2) file>>nfile>>dummy;
    else file>>dummy>>nfile>>dummy;

    if (dummy==-1000) {
        int j;
        for (i=0;i<nfile;i++)
            for (j=0;j<nfile;j++) file>>dummy;
    }

    s_expo X(0);

    while (file) {

		if (type==3) {
            file>>modqcoef;
            if (modqcoef==modqzero) break;
		}
		else  {
            file>>intcoef;
            if (intcoef==intzero) break;
		}
		for (i=0;i<nfile;i++) file>>X[i];
        if (type==3) fmodq.insert(modqcoef,X);
        else fint.insert(std::move(intcoef),X);

		do {
            if (type==3) {
                file>>modqcoef;
                if (modqcoef==modqzero) break;
            }
            else  {
                file>>intcoef;
                if (intcoef==intzero) break;
            }
            for (i=0;i<nfile;i++) file>>X[i];
            if (type==3) fmodq.insert(modqcoef,X);
            else fint.insert(std::move(intcoef),X);

		}
		while (1);

		if (type==3) {
            if (!fmodq.iszero()) {
                Hmodq.sinsert(std::move(fmodq));
                nH++;
            }
		}
		else {
            if (!fint.iszero()) {
                Hint.sinsert(std::move(fint));
                nH++;
            }
		}
    }
    WriteText("\nH red. Elements in H: " + ToString(nH));

    Chrono();
    if (type==3) is=isinG(Hmodq,Gmodq,1);
    else is=isinG(Hint,Gint);
    WriteText("\n\n Time: " + ToString(Chrono()));

    if (is) {
        WriteText("\n\n H is in <G>!");
    }
    else WriteText("\n\n H is not in <G>!");

    WriteText("\nReductions: ");
    if (type==3) WriteText(ToString(s_polynom<modq>::nreduction));
    else WriteText(ToString(s_polynom<integer>::nreduction));

    return is;
}


int isinG(char* hfile,char* gbfile)
{
    return isinG(std::string(hfile),std::string(gbfile));
}


//

void restoreGp(val::Glist< val::s_polynom<val::modq> > &Gp,const val::Glist< val::s_polynom<val::integer> > &G)
{
    val::GlistIterator<val::s_polynom<val::integer> > pG;
    Gp.resetactual();

    for (pG=G;pG;pG++,Gp.moveactual()) {
        Gp.actualvalue().convertfrom(pG());
        Gp.actualvalue().normalize();
    }
}

//

template <>
val::s_polynom<val::modq> spol(const val::s_polynom<val::modq>& g1,const val::s_polynom<val::modq>& g2)
{
 val::s_polynom<val::modq> s;
 val::modq eins(1);

 if (g1.iszero() || g2.iszero()) {
	 return s;
 }
 val::s_expo X=val::lcm(g1.LT(),g2.LT());
 s=g1;
 s*=(X/g1.LT());
 X/=g2.LT();
 s.minusmalmonom(g2,eins,X);
 return s;
}



template <>
val::s_polynom<val::integer> spol(const val::s_polynom<val::integer>& g1,const val::s_polynom<val::integer>& g2)
{
 val::integer gcd;
 val::s_polynom<val::integer> s;

 if (g1.iszero() || g2.iszero()) return s;

 val::s_expo X(val::lcm(g1.LT(),g2.LT()));

 gcd=val::ggTspez(g1.LC(),g2.LC());
 s=g1;

 s.malmonom(val::EDIV(g2.LC(),gcd),X/g1.LT());
 s.minusmalmonom(g2,val::EDIV(g1.LC(),gcd),X/g2.LT());

 s.normalize();
 return s;
}

//

int homogenize(const std::string& i_name,const std::string& o_name,int comment)
{
    int i=0,washomogen=1,ismodq=0,order,n,d,nG=0;

	val::s_polynom<val::rational> f;
	val::Glist<val::s_polynom<val::rational> > P;
	std::string line;

	std::ifstream i_file(i_name,std::ios::in);
	if (!i_file) {
		MyMessage("\nFILE DOES NOT EXIST\n\n");
        return 0;
	}

    do {
        getline(i_file,line);
        if (line!="") i++;
        else break;
	}
	while(i_file);

	if (i!=2 && i!=3) {
        MyMessage("\nWrong type of file!\n\n");
        return 0;
	}
	else if (i==3) ismodq=1;

	i_file.clear();
	i_file.seekg(0,std::ios::beg);

	if (ismodq) i_file>>ismodq>>n>>order;
	else i_file>>n>>order;

	val::s_expo::setordtype(order);

	if (order==-1000) {
		val::matrix<int> Mold(n);
		int j;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				i_file>>Mold(i,j);
			}
		}
		val::s_expo::setordmatrix(Mold);
	}

	else if (order!=0 && order!=-2 && order!=-1) val::s_expo::setordtype(-2);

	val::s_expo::setdim(n+1);
	val::s_expo X(0);
	val::rational coeff,zero;

	do {
        i_file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				i_file>>X[i];
			}
			f.insert(std::move(coeff),X);
			d=X.totdeg();
		}
		else break;
		do {
			i_file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					i_file>>X[i];
				}
				if (X.totdeg()!=d) washomogen=0;
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			P.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (i_file);
	i_file.close();
	if (comment) WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.\n");

    if (washomogen) {
        val::s_expo::setdim(n);
    }
    else {
         //homogenize:
        int n=val::s_expo::getdim(),s_order=val::s_expo::getordtype();
        val::matrix<int> Mnew,Ms;
        if (s_order==-1000) {
            Ms=std::move(val::s_expo::getordmatrix());
            int j;
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
        for (P.resetactual();P.actualvalid();P.moveactual()) {P.actualvalue().homogenize();P.actualvalue().reord();}
        P.sort();
    }

    std::ofstream file(o_name,std::ios::out | std::ios::trunc);

    if (ismodq) file<<ismodq<<std::endl<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;
    else file<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;

    if (val::s_expo::getordtype()==-1000) {
        int dim=val::s_expo::getdim(),j;
        file<<std::endl;
        for (i=0;i<dim;i++) {
            for (j=0;j<dim;j++) file<<val::s_expo::getordmatrix()(i,j)<<" ";
            file<<std::endl;
        }
    }
    file<<std::endl;

    for (P.resetactual();P.actualvalid();P.moveactual()) {
        file<<P.getelement()<<std::endl;
    }

    file<<'0'<<std::endl;
    file.close();

    if (comment) WriteText("\n Basis written in  " + o_name + "!\n\n");
    return 1;
}


int dehomogenize(const std::string& i_name,const std::string& o_name,int comment)
{
    int i=0,ismodq=0,order,n,d,nG=0;

	val::s_polynom<val::rational> f;
	val::Glist<val::s_polynom<val::rational> > P;
	std::string line;

	std::ifstream i_file(i_name,std::ios::in);
	if (!i_file) {
		MyMessage("\nFILE DOES NOT EXIST\n\n");
        return 0;
	}

    do {
        getline(i_file,line);
        if (line!="") i++;
        else break;
	}
	while(i_file);

	if (i!=2 && i!=3) {
        MyMessage("\nWrong type of file!\n\n");
        return 0;
	}
	else if (i==3) ismodq=1;

	i_file.clear();
	i_file.seekg(0,std::ios::beg);

	if (ismodq) i_file>>ismodq>>n>>order;
	else i_file>>n>>order;

	if (n<=1) {
        MyMessage("\nPolynomials have only one variable!\n\n");
        return 0;
	}

	if (order==-1000) {
		val::matrix<int> Mold(n),Mnew(0,n-1,n-1);
		int j;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				i_file>>Mold(i,j);
			}
		}
		for (i=1;i<n;i++) {
            for (j=0;j<n-1;j++) Mnew(i-1,j) = Mold(i,j);
		}
		val::s_expo::setordmatrix(Mnew);
	}
	else if (order==0 || order==-1) order=-1;
    else order=-2;

    val::s_expo::setordtype(order);
	val::s_expo::setdim(n-1);
	val::s_expo X;
	val::rational coeff,zero;
	n--;

	do {
        i_file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				i_file>>X[i];
			}
			i_file>>d;
			f.insert(std::move(coeff),X);
		}
		else break;
		do {
			i_file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					i_file>>X[i];
				}
				i_file>>d;
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			P.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (i_file);
	i_file.close();
	if (comment) WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.\n");

    std::ofstream file(o_name,std::ios::out | std::ios::trunc);

    if (ismodq) file<<ismodq<<std::endl<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;
    else file<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;

    if (val::s_expo::getordtype()==-1000) {
        int dim=val::s_expo::getdim(),j;
        file<<std::endl;
        for (i=0;i<dim;i++) {
            for (j=0;j<dim;j++) file<<val::s_expo::getordmatrix()(i,j)<<" ";
            file<<std::endl;
        }
    }
    file<<std::endl;

    for (P.resetactual();P.actualvalid();P.moveactual()) {
        file<<P.getelement()<<std::endl;
    }

    file<<'0'<<std::endl;
    file.close();

    if (comment) WriteText("\n Basis written in " + o_name + "!\n\n");

    return 1;
}


int preparenormalpos(const std::string& i_name,const std::string& o_name,int comment)
{
    int i=0,ismodq=0,order,n,nG=0;

	val::s_polynom<val::rational> f;
	val::Glist<val::s_polynom<val::rational> > P;
	std::string line;

	std::ifstream i_file(i_name,std::ios::in);
	if (!i_file) {
		MyMessage("\nFILE DOES NOT EXIST\n\n");
        return 0;
	}

    do {
        getline(i_file,line);
        if (line!="") i++;
        else break;
	}
	while(i_file);

	if (i!=2 && i!=3) {
        MyMessage("\nWrong type of file!\n\n");
        return 0;
	}
	else if (i==3) ismodq=1;

	i_file.clear();
	i_file.seekg(0,std::ios::beg);

	if (ismodq) i_file>>ismodq>>n>>order;
	else i_file>>n>>order;


	if (order==-1000) {
		val::matrix<int> Mold(n),Mnew(0,n+1,n+1);
		int j;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) {
				i_file>>Mold(i,j);
			}
		}
		Mnew(0,0)=1;
        for (i=0;i<n;i++) {
            for (j=0;j<n;j++) Mnew(i+1,j+1) = Mold(i,j);
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
	val::s_expo X(0);
	val::rational coeff,zero;
	val::RandomClass Rand;

	// Set first polynomial:
	X=val::s_expo(0);
	X[0]=1;
    f.insert(val::rational(1),X);
    for (i=1;i<=n;i++) {
        X=val::s_expo(0);
        X[i]=1;
        f.insert(val::rational(Rand.random(-50,50)),X);
    }

    P.sinsert(std::move(f));

	// Set rest:
	do {
        i_file>>coeff;
		if (coeff!=zero) {
			for (i=0;i<n;i++) {
				i_file>>X[i+1];
			}
			f.insert(std::move(coeff),X);
		}
		else break;
		do {
			i_file>>coeff;
			if (coeff!=zero) {
				for (i=0;i<n;i++) {
					i_file>>X[i+1];
				}
				f.insert(std::move(coeff),X);
			}
			else break;
		}
		while (1);

		if (!f.iszero()){
			P.sinsert(std::move(f));
			nG++;
		}
		else break;
	}
	while (i_file);
	i_file.close();
	if (comment) WriteText("\nRed list G from file.\nG has " + val::ToString(nG) + " elements.\n");

    std::ofstream file(o_name,std::ios::out | std::ios::trunc);

    if (ismodq) file<<ismodq<<std::endl<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;
    else file<<val::s_expo::getdim()<<std::endl<<val::s_expo::getordtype()<<std::endl;

    if (val::s_expo::getordtype()==-1000) {
        int dim=val::s_expo::getdim(),j;
        file<<std::endl;
        for (i=0;i<dim;i++) {
            for (j=0;j<dim;j++) file<<val::s_expo::getordmatrix()(i,j)<<" ";
            file<<std::endl;
        }
    }
    file<<std::endl;

    for (P.resetactual();P.actualvalid();P.moveactual()) {
        file<<P.getelement()<<std::endl;
    }

    file<<'0'<<std::endl;
    file.close();

    if (comment) WriteText("\n Basis written in " + o_name + "!\n\n");

    return 1;
}

//

template <>
void reduceKritPairs(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::modq> > &G,
                     val::Glist<val::s_polynom<val::modq> > &Gd,int nd,int)
{
    val::s_polynom<val::modq> *f1,*f2,*mh,f;
    val::s_polynom<val::integer> *h;
    int *s_done;
    do {
        ProtectData.lock();
        while (P.actualvalid() && P.actualvalue().s_done) P.moveactual();
        if (!P.actualvalid()) {
            ProtectData.unlock();
            return;
        }
        f1=P.actualvalue().mf1;f2=P.actualvalue().mf2;h=P.actualvalue().ih;mh=P.actualvalue().mh;
        s_done = &(P.actualvalue().s_done);
        P.moveactual();
        ProtectData.unlock();
        if (nd && KritPairs::nGd==nd) {*s_done=1;return;}
        if (h!=NULL) {
            f.convertfrom(*h);
            f.normalize();
        }
        else if (mh!=NULL) f=std::move(*mh);
        else {
            f=spol(*f1,*f2);
        }
        if (nd && KritPairs::nGd==nd)  {
            *s_done=1;
            return;
        }
        f.reduction(G,0);
        if (!f.iszero()) {
            KritPairs::protectGd.lock();
            if (nd && KritPairs::nGd==nd) {KritPairs::protectGd.unlock();*s_done=1;return;}
            f.reduction(Gd,1);
            if (!f.iszero()) {
                Gd.sinsert(std::move(f));
                KritPairs::nGd++;
            }
            else *s_done=1;
            KritPairs::protectGd.unlock();
        }
        else *s_done=1;
    }
    while (1);
}



template <>
void reduceKritPairs(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::integer> > &G,
                     val::Glist<val::s_polynom<val::integer> > &Gd,int nd,int checkdiv)
{
    val::s_polynom<val::integer> *f1,*f2,*h,f;
    val::integer q(val::modq::q), zero;

    do {
        ProtectData.lock();
        while (P.actualvalid() && P.actualvalue().s_done) P.moveactual();
        if (!P.actualvalid()) {
            ProtectData.unlock();
            return;
        }
        f1=P.actualvalue().if1;f2=P.actualvalue().if2;h=P.actualvalue().ih;
        P.moveactual();
        ProtectData.unlock();
        if (nd && KritPairs::nGd==nd) {
            return;
        }
        if (h!=NULL) {
            f=std::move(*h);
        }
        else {
            f=spol(*f1,*f2);
        }
        if (nd && KritPairs::nGd==nd) {
            return;
        }
        f.reduction(G,0);
        if (!f.iszero()) {
            KritPairs::protectGd.lock();
            if (nd && KritPairs::nGd==nd) {
                KritPairs::protectGd.unlock();
                return;
            }
            f.reduction(Gd,1);
            if (!f.iszero()) {
                if (checkdiv && f.LC()%q==zero) KritPairs::p_teilt=1;
                Gd.sinsert(std::move(f));
                KritPairs::nGd++;
            }
            KritPairs::protectGd.unlock();
        }
    }
    while (1);
}


void reduceKritPairs2(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::integer> > &G,
                      val::Glist<val::s_polynom<val::integer> >&Gd,
                      val::Glist<val::s_polynom<val::integer> > &Hd,int nd)
{
    val::s_polynom<val::integer> *f1,*f2,*h,f;
    val::integer q(val::modq::q), zero;

    do {
        ProtectData.lock();
        while (P.actualvalid() && P.actualvalue().s_done) P.moveactual();
        if (!P.actualvalid()) {
            ProtectData.unlock();
            return;
        }
        f1=P.actualvalue().if1;f2=P.actualvalue().if2;h=P.actualvalue().ih;
        P.moveactual();
        ProtectData.unlock();
        if (KritPairs::nGd==nd) {
            return;
        }
        if (h!=NULL) {
            f=std::move(*h);
        }
        else {
            f=spol(*f1,*f2);
        }
        if (KritPairs::nGd==nd) {
            //WriteText(" HK_r ");
            return;
        }

        f.reduction(G,0);
        protectHd.lock();
        Hd.inserttoend(std::move(f));
        KritPairs::nGd++;
        protectHd.unlock();

        do {
            if (KritPairs::protectGd.try_lock() ) {
                Hdisbusy=true;
                protectHd.lock();
                while (!Hd.isempty()) {
                    f=std::move(Hd.actualvalue());
                    Hd.skiphead();
                    protectHd.unlock();
                    f.reduction(Gd,0);
                    if (f.LC()%q==zero) KritPairs::p_teilt=1;
                    Gd.sinsert(std::move(f)); WriteText(".");
                    protectHd.lock();
                }
                protectHd.unlock();
                Hdisbusy=false;
                KritPairs::protectGd.unlock();
                break;
            }
            else if (Hdisbusy) break;
        }
        while (1);

        if (KritPairs::nGd==nd) return;

    }
    while (1);
}

/*
// sequential : 1 thread
void reduceKritPairs_seq(val::Glist<KritPairs> &P,val::Glist<val::s_polynom<val::integer> > &G,
                      val::Glist<val::s_polynom<val::integer> > &Gd,int nd)
{
    val::s_polynom<val::integer> f;
    val::integer q(val::modq::q), zero;

    do {
        while (P.actualvalid() && P.actualvalue().s_done) P.moveactual();
        if (!P.actualvalid()) {
            return;
        }

        if (P.actualvalue().ih!=NULL) {
            f=std::move(*(P.actualvalue().ih));
        }
        else {
            f=spol(*(P.actualvalue().if1),*(P.actualvalue().if2));
        }
        P.moveactual();

        f.reduction(G,0);
        f.reduction(Gd,1);
        if (f.LC()%q==zero) KritPairs::p_teilt=1;
        Gd.sinsert(std::move(f));
        KritPairs::nGd++;
        if (KritPairs::nGd==nd) return;
    }
    while (1);
}
*/


template<>
void update(val::s_polynom<val::integer> &f,val::Glist< val::s_polynom<val::integer> > &G,val::Glist< common_bb::spair > &lspair,int m)
{
 int i,j;
 val::s_expo t1;
 char *h_done;
 val::s_polynom<val::integer> *hint;
 val::GlistManipulator<val::s_polynom<val::integer> > p,q;
 val::GlistManipulator<spair> sp;
 spair hs;


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
	 t1 = val::lcm(sp().fint->LT(),sp().gint->LT());

	 if ( (hint->LT()|t1) && (val::lcm(sp().fint->LT(),hint->LT())!=t1)
		 && (val::lcm(hint->LT(),sp().gint->LT())!=t1) ) {
		 sp().s_done=1;
	 }
 }

  // Insertion in dPairs:
 for (i=0,p=G;i<m;i++,p++) {
	 if (h_done[i]) continue;
	 hs.fint=&(p());hs.gint=hint;
     hs.deg= (val::lcm(p().LT(),hint->LT())).totdeg();
	 lspair.sinsert(hs);
 }
 delete[] h_done;
}


template<>
void update(val::s_polynom<val::modq> &f,val::Glist< val::s_polynom<val::modq> > &G,val::Glist< common_bb::spair > &lspair,int m)
{
 int i,j;
 val::s_expo t1;
 char *h_done;
 val::s_polynom<val::modq> *hmodq;

 val::GlistManipulator<val::s_polynom<val::modq> > p,q;
 val::GlistManipulator<spair> sp;
 spair hs;


 if (f.iszero()) return;

 if (G.isempty()) return;

 h_done = new char[m];
 for (i=0;i<m;i++) h_done[i]=0;

 hmodq= &f;

 for (i=0,p=G;i<m;i++,p++) {
	 if (h_done[i]) continue;
	 if (lcmdis(hmodq->LT(),p().LT(),t1)) {
		 h_done[i]=2;
		 continue;
	 }
	 for (j=0,q=G;j<m;j++,q++) {
		 if (j==i || h_done[j]==1) continue;
		 if (val::lcm(hmodq->LT(),q().LT())|t1) {
			 h_done[i]=1;
			 break;
		 }
	 }
 }


 for (sp=lspair;sp;sp++) {
	 if (sp().s_done) continue;
	 t1 = val::lcm(sp().fmodq->LT(),sp().gmodq->LT());

	 if ( (hmodq->LT()|t1) && (val::lcm(sp().fmodq->LT(),hmodq->LT())!=t1)
		 && (val::lcm(hmodq->LT(),sp().gmodq->LT())!=t1) ) {
		 sp().s_done=1;
	 }
 }


  // Insertion in dPairs:
 for (i=0,p=G;i<m;i++,p++) {
	 if (h_done[i]) continue;
	 hs.fmodq=&(p());hs.gmodq=hmodq;
     hs.deg= (val::lcm(p().LT(),hmodq->LT())).totdeg();
	 lspair.sinsert(hs);
 }
 delete[] h_done;
}



void update(val::s_polynom<val::integer> &fint,val::s_polynom<val::modq> &fmodq,val::Glist< val::s_polynom<val::integer> > &G,
            val::Glist< val::s_polynom<val::modq> > &Gp,val::Glist<spair> &lspair,int m)
{
 int i,j;
 val::s_expo t1;
 char *h_done;
 val::s_polynom<val::integer> *hint;
 val::s_polynom<val::modq> *hmodq;

 val::GlistManipulator<val::s_polynom<val::integer> > p,q;
 val::GlistManipulator<val::s_polynom<val::modq> > pGp;
 val::GlistManipulator<spair> sp;
 spair hs;


 if (fint.iszero()) return;

 if (G.isempty()) return;

 h_done = new char[m];
 for (i=0;i<m;i++) h_done[i]=0;

 hint= &fint;

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
			 break;  // Verlasse j-Schleife
		 }
	 }
 }


 for (sp=lspair;sp;sp++) {
	 if (sp().s_done) continue;
	 t1 = val::lcm(sp().fint->LT(),sp().gint->LT());

	 if ( (hint->LT()|t1) && (val::lcm(sp().fint->LT(),hint->LT())!=t1)
		 && (val::lcm(hint->LT(),sp().gint->LT())!=t1) ) {
		 sp().s_done=1;
	 }
 }


 hmodq = &fmodq;
  // Einfuegen in dPairs:
 for (i=0,p=G,pGp=Gp;i<m;i++,p++,pGp++) {
	 if (h_done[i]) continue;
	 hs.fint=&(p());hs.gint=hint;hs.fmodq=&(pGp());hs.gmodq=hmodq;
     hs.deg= (val::lcm(p().LT(),hint->LT())).totdeg();
	 lspair.sinsert(hs);
 }
 delete[] h_done;
}


//

} // end namespace common_bb
