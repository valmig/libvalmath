

#include <s_groebner/walk.h>


namespace walk
{



void transformweight(const val::vector<val::rational> &wr,val::vector<val::integer> &wi)
{
 using namespace val;
 int i,n=wi.dimension();
 integer kgN;

 // Get lcm of denominators of v[i]
 kgN=denominator(wr(0));
 for (i=1;i<n;i++) kgN=lcm(kgN,denominator(wr[i]));

 // Evtl. one could check if numbers are to high for int type.
 for (i=0;i<n;i++) wi[i]=nominator(wr[i]) * val::EDIV(kgN,denominator(wr(i)));
}


void normalize(val::vector<val::integer> &w)
{
 int i,n=w.dimension();
 val::integer ggt,eins(1),zero;

 ggt=w[0];
 for (i=1;i<n;i++)
	 if (ggt==eins) return;
	 else ggt=val::ggTspez(ggt,w[i]);
 if (ggt!=zero)
	 for (i=0;i<n;i++) w[i].EDIVBY(ggt);
 return;
}


val::matrix<int> Createordmatrix(int ord)
{
    int i,n=val::s_expo::getdim();
    val::matrix<int> M(0,n,n);

    if (ord==-1) for (i=0;i<n;i++) M(i,i) =1;
    else if (ord==0) {
        for (i=0;i<n;i++) M(0,i) = 1;
        for (i=1;i<n;i++) M(i,i-1) =1;
    }
    else {  // ord =-2
        for (i=0;i<n;i++) M(0,i) = 1;
        for (i=1;i<n;i++) M(i,n-i) =-1;
    }
    return M;
}


void perturbweight(val::vector<val::integer> &w,int deg,const val::matrix<int> &A,int k)
{
 using namespace val;
 int i,j,m=0,n=w.dimension(),max;
 rational eps,faktor(1);

 if (k==1) {
	 for (i=0;i<n;i++) w[i] = integer(A(0,i));
 }
 else {
	 val::vector<rational> wh(n);
	 for (i=0;i<n;i++) wh[i]=rational(A(0,i));
	 if (k>n) k=n;
	 // Presuming  Maxcoef(A)=1 !
	 for (i=1;i<k;i++) {
		 max=1;
		 for (j=0;j<n;j++)
			 if (max < hilfinteger::abs(A(i,j))) max = hilfinteger::abs(A(i,j));
		 m+=max;
	 }
	 eps=rational(1,deg*m+1);
	 for (i=1;i<k;i++) {
		 faktor*=eps;
		 for (j=0;j<n;j++) wh[j]+=faktor*rational(A(i,j));
	 }
	 transformweight(wh,w);
 }
 normalize(w);
}

void setpositiv(int &wnewpos,int &wspos,int &pwtpos,const val::vector<val::integer> &ws,const val::vector<val::integer> &pwt)
{
 int i,n=ws.dimension();
 val::integer Null;

 wnewpos=wspos=pwtpos=1;
 for (i=0;i<n;i++) {
	 if (ws[i]==Null) wspos=0;
	 if (pwt[i]==Null) pwtpos=0;
	 if (pwt[i]==Null && ws[i]==Null) wnewpos=0;
 }
}


int ispositv(const val::vector<val::integer> &w)
{
    int n=w.dimension();
    for (int i=0;i<n;i++) {
        if (w[i].isNull()) return 0;
    }
    return 1;
}


void getnextvector(val::vector<val::integer> &wold,val::vector<val::integer> &wnew,const val::rational& u)
{
 using namespace val;
 int i,n=wold.dimension();
 integer gcd,eins(1);

 for (i=0;i<n;i++) wnew[i]=(denominator(u)-nominator(u)) * wold[i]+nominator(u)*wnew[i];
 gcd=wnew[0];
 for (i=1;i<n;i++)
	 if (gcd==eins) break;
	 else gcd=ggTspez(gcd,wnew[i]);

 if (gcd==eins) return;
 for (i=0;i<n;i++) wnew[i].EDIVBY(gcd);
 return;
}


template <>
void walkmain(val::Glist<val::s_polynom<val::modq> > &G,int neworder,int k,int l,const val::matrix<int> &M,int comment)
{
    val::matrix<int> Mt=M;
    GroebnerWalk(G,neworder,Mt,k,l,comment);
}


template <>
void walkmain(val::Glist<val::s_polynom<val::integer> > &G,int neworder,int k,int l,const val::matrix<int> &M,int comment)
{
    val::matrix<int> Mt=M;
    GroebnerWalk(G,neworder,Mt,k,l,comment);
    if (comment) common_bb::WriteText("\nMaximal integer-length in G: " + val::ToString(common_bb::Maximalinteger(G)));
}


} // end namespace
