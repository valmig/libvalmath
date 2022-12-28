
#include <s_groebner/hilbert.h>
#include <val_basics.h>
#include <s_groebner/common_bb.h>


namespace hilbert
{

int isdisjunct(const val::s_expo& a,const val::s_expo& b)
{
 int i,n=val::s_expo::getdim();

 for (i=0;i<n;i++) if (val::Min(a[i],b[i])!=0) return 0;
 return 1;
}


int notdisjunct(const val::s_expo& a,const val::s_expo& b)
{
 int i,n=val::s_expo::getdim();

 for (i=0;i<n;i++)
	 if (a[i] && b[i]) return 1;
 return 0;
}


// c = lcm(a,b)/b, Vor.: c=0
void kgVdiv(const val::s_expo &a,const val::s_expo &b,val::s_expo &c)
{
 int i,n=val::s_expo::getdim();

 for (i=0;i<n;i++)
	 if (a[i]>b[i]) c[i]=a[i]-b[i];
	 else c[i]=0;
 return;
}


// Prem:: vector is in increasing  order (ordtype=-1) sorted
int reduce(const val::vector<val::s_expo> &feld,val::vector<char> &iszero)
{
 int i,j,anz=0,m=iszero.dimension();

 for (i=0;i<m;i++) {
	 if (iszero[i]) continue;
	 else anz++;
	 for (j=i+1;j<m;j++) {
		 if (iszero[j]) continue;
		 if (feld[i]|feld[j]) {iszero[j]=1;}
	 }
 }
 return anz;
}


// Prem:: vector is in increasing  order (ordtype=-1) sorted
val::vector<val::s_expo> reduce(const val::vector<val::s_expo> &feld)
{
 int i,j,anz=0,m=feld.dimension();
 val::vector<char> iszero(0,m);

 for (i=0;i<m;i++) {
	 if (iszero[i]) continue;
	 else anz++;
	 for (j=i+1;j<m;j++) {
		 if (iszero[j]) continue;
		 if (feld[i]|feld[j]) {iszero[j]=1;}
	 }
 }

 val::vector<val::s_expo> I2(anz);
 for (i=j=0;i<m;i++)
	 if (iszero[i]) continue;
	 else {
		 I2[j]=feld[i];
		 j++;
	 }
 return I2;
}


// Prem:: vector is in increasing  order (ordtype=-1) sorted
// Computes Hilbert-numerator (unreduced) of K[X]/I
val::pol<val::integer> Hilbertnum(const val::vector<val::s_expo> &I,int anzI,int MCIcrit)
{
 int deg,i,j,crit=1;
 val::vector<char> iszero;
 val::vector<val::s_expo> Idivm,I2;
 val::pol<val::integer> f,f1;

 // Basis-cases:
 if (anzI==0) return (val::pol<val::integer>(val::integer(1),0));
 if (anzI==1) return (val::pol<val::integer>(val::integer(1),0)-val::pol<val::integer>(val::integer(1),I[0].totdeg()));

 if (anzI<=val::s_expo::getdim() && MCIcrit) {
	 int MCI=1;
	 for (i=0;i<anzI && MCI;i++){
		 for (j=i+1;j<anzI;j++) {
			 if (notdisjunct(I[i],I[j])) {
				 if ((i<anzI-1) && (j<anzI-1)) crit=0;
				 MCI=0;
				 break;
			 }
		 }
	 }
	 if (MCI) {
		 f=val::pol<val::integer>(val::integer(1),0) - val::pol<val::integer>(val::integer(1),I[0].totdeg());
		 for (i=1;i<anzI;i++) 
			 f-=f.multbypower(I[i].totdeg());
		 return f;
	 }
 }



 Idivm= val::vector<val::s_expo>(anzI-1);
 for (i=0;i<anzI-1;i++) kgVdiv(I[i],I[anzI-1],Idivm[i]);
 deg=I[anzI-1].totdeg();
 Idivm.sort();

 I2=reduce(Idivm);
 f=Hilbertnum(I,anzI-1,crit);
 f1=Hilbertnum(I2,I2.dimension());
 f1.getmultbypower(deg);
 f-=f1;
 return f;
}


// Computes (c,d) with: f(z)-q(z) = c*z^d + a1*z^(d+1)+...
void getlowestmon(const val::pol<val::integer> &f,const val::pol<val::integer> &g,int &c,int &d)
{
 val::pol<val::integer> h;

 h=f-g;
 if (h.iszero()) {
    c=d=0;
    return;
 }
 val::polIterator<val::integer> ph;
 ph=h;
 ph.moveactualtolast();
 c=int(ph.actualcoef());
 d=ph.actualdegree();
}


val::pol<val::rational> binomial(int z,int d)
{
 int i;

 val::pol<val::rational> f(val::rational(1)); // => f = 1 .

 for (i=0;i<d;i++)
	 f*=val::rational(1,i+1)*(val::pol<val::rational>(val::rational(1),1)+val::pol<val::rational>(val::rational(z-i,1),0));

 return f;
}


//Prem.: I is reduced and sorted ( -1 )
val::pol<val::rational> Hilbertpolynomial(const val::vector<val::s_expo> &I, int affin)
{
 using namespace val;
 int i,d,e=0,N,order=s_expo::getordtype();
 pol<integer> f,g;
 pol<rational> h;

 if (I.dimension()==0) return h;

 s_expo::setordtype(-1);
 f=Hilbertnum(I,I.dimension());
 s_expo::setordtype(order);

 g=pol<integer>(integer(1),0)-pol<integer>(integer(1),1);

 while (f.eval(integer(1))==integer(0)) {
	 f/=g;
	 e++;
 }
 d=s_expo::getdim()-e;
 if (affin) d++;
 N=deg(f);

 if (d==0) return h;

 for (i=0;i<=N;i++)
	 h+=rational(integer(f[i]),integer(1))*binomial(d-i-1,d-1);

 return h;
}


val::pol<val::rational> Hilbertpolynomial(const std::string &filename,int affin)
{
    using namespace val;

    pol<rational> h;
    std::ifstream file(filename,std::ios::in);
    if (!file) {
        common_bb::MyMessage("\nCANNOT READ FILE!!!!");
        return h;
    }

    int i=0,n,q,ord,nG=0;
    std::string line;
    integer coef;

    do {
        std::getline(file,line);
        if (line!="") i++;
        else break;
	}
	while(file);

	if (i!=2 && i!=3) {
        common_bb::MyMessage("\nWrong type of file!");
        return h;
	}

	file.clear();
	file.seekg(0,std::ios::beg);

	if (i==2) file>>n>>ord;
	else file>>q>>n>>ord;

	s_expo::setdim(n);
	s_expo::setordtype(-1);
	s_expo X;
	Glist<s_expo> G;

	if (ord==-1000) {
        int j;
        for (i=0;i<n;i++)
            for (j=0;j<n;j++) file>>q;
	}
	while (file) {
		file>>coef;
		if (coef!=0) {
            file>>X;
            G.inserttoend(X);
            nG++;
		}
		else break;
		do {
			file>>coef;
			if (coef!=0){
			    file>>X;
			}
		}
		while (coef!=0);
	}

	if (G.isempty()) return h;

	val::vector<s_expo> I(nG);
	for (i=0;i<nG;i++) I[i] = std::move(G[i]);

	I.sort();
	h=Hilbertpolynomial(I,affin);
	return h;
}


} // end namespace
