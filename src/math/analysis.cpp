
#ifndef ANALYSIS_CPP
#define ANALYSIS_CPP

#include <analysis.h>
#include <vector.h>
#include <val_basics.h>
#include <complex.h>


namespace val
{


const double PI = 3.1415926535897932384626433832795;

template <>
std::atomic<int> term<double>::mnumber(0);



// constants:

namespace analysis
{
// error bound
double fehler(1e-8);
// ln2: ------------------------------------------------------------------------
static const double MLN2 = 0.6931471805599453094172321214581766; // ln 2
// ln 2 split into two parts
static const double L2U = 0.69314718055966295651160180568695068359375;
static const double L2L = 0.28235290563031577122588448175013436025525412068e-12;
//-------------------------------------------------------------------------------

static const double PI4_A= 0.7853981554508209228515625; // PI/4 split into three parts
static const double PI4_B= 0.794662735614792836713604629039764404296875e-8;
static const double PI4_C =0.306161699786838294306516483068750264552437361480769e-16;
// 4/PI
static const double M_4_PI = 1.273239544735162542821171882678754627704620361328125;  // 4/PI


int isclosezero(const pol<double> &f,const double &eps = 1e-9);

pol<double> gcd(const pol<double> &f,const pol<double> &g, const double &eps = 1e-9);

const pol<double>& RoundPol(pol<double>& f,const double &eps=1e-9);

const pol<double>& divbylinearfaktors(pol<double> &f,const val::vector<double> &Root);

// Generates a Sturm-chain to polynomial f. Requires: gcd(f,f') = 1;
// Returns length of chain - 1.
int Sturm(const pol<double>& f,val::vector<pol<double> > &p);


// Returns number of sign-changes of the Sturm-chain p at x.
int Wechsel(const val::vector<pol<double> > &p,const double &x);


// Localize all intervals, where f has only one root.
void lokalisiere(const val::vector<pol<double> > &p,val::vector<double[]> &a,const double &a1,const double &a2,int& gefunden,const double &eps=1e-9);

// -------------- Fktn zur Berechnung komplexer Nstn -------------------------------------------------------

// Computes y, so that the polynomial x^2 - sx -t has a root x = 1/2 * s +- y.
// Returns if the root is real or not.
int getrootdeg2(const double &s,const double &t,double &y);


// division of f by x^2 -sx -t.
void divbyquadratfaktor(pol<double> &f,const double& s, const double& t, double eps=1.-9);


// double Horner - method.
void dhorner(const val::vector<double> &a,val::vector<double> &b,const double& s,const double& t);

int bairstow(const pol<double> &f, double &s, double &t,int N=100,double E=1.e-9);

// Check if  f = aX^2n + b , where  ab>0.
int iscomplexcyclomatic(const pol<double>& f);




// ----------------------------------------------------------




// ---------------------------------------------------------------------------------------------------

int isclosezero(const pol<double> &f,const double &eps)
{
	
	for (const auto & m : f) {
		if (abs(m.actualcoef()) >= eps) return 0;
	}	
	return 1;
}

pol<double> gcd(const pol<double> &f,const pol<double> &g, const double &eps)
{
	pol<double> akt,rest,vorrest;

	if (deg(f)<deg(g)) {akt=g;rest=f;}
	else {akt=f;rest=g;}
	while (!isclosezero(rest))
	{
		vorrest=std::move(rest);
		rest= akt % vorrest;
		akt= std::move(vorrest);
		//std::cout<<"\n rem: \n"<<rest<<std::endl;
	}
	akt.normalize();
	return akt;
}


const pol<double>& RoundPol(pol<double>& f,const double &eps)
{
 if (f.iszero()) return f;

 int i=f.degree();

 for (;i>=0;i--) if (abs(f[i])>=eps) break;
 if (i==-1) {
    f=pol<double>();
    return f;
 }
 pol<double> g(f[i],i);
 i--;
 //g[i--]=f[i];
 for (;i>=0;i--)
     if (abs(f[i])>=eps) g+=pol<double>(f[i],i);
 f=std::move(g);
 return f;
}



const pol<double>& divbylinearfaktors(pol<double> &f,val::vector<double> &Root)
{
    if (f.iszero()) return f;
    int i,l=Root.dimension();
    pol<double> g;

    for (i=0;i<l;i++){
        g=pol<double>(1.0,1) - pol<double>(Root[i],0);
        f/=g;
    }
    return RoundPol(f);
}



int Sturm(const pol<double>& f,val::vector<pol<double> > &p)
{
 int i=1;

 p=val::vector<pol<double> >();
 if (f.iszero()) return 0;
 p=val::vector<pol<double> >(deg(f)+1);

 p[0]=f;
 p[1]=f.derive();

 while ( deg(p[i])>0 ) {
       i++;
       p[i]=-(p[i-2]%p[i-1]);
 }
 return i;
}


// Number of changes of sign of the Sturm-Chain p at x.
int Wechsel(const val::vector<pol<double> > &p,const double &x)
{
 if (p.isempty()) return 0;

 int i,n=0,l=p.dimension()-1;
 int vz1,vz2;

 vz1= (p[0].eval(x)>0)? 1:-1;
 for (i=1;i<=l;i++)
     if (abs(p[i].eval(x))>0) {
        vz2= (p[i].eval(x)>0)? 1:-1;
        if (vz2!=vz1) n++;
        vz1=vz2;
     }

 return n;
}


void lokalisiere(const val::vector<pol<double> > &p,val::vector<double[2]> &a,const double &a1,const double &a2,int& gefunden,const double &eps)
{
 int i,n;
 double b,c,d;

 if (p.isempty() || a.isempty()) return;
 if (gefunden == a.dimension()) return;

 n=Wechsel(p,a1)-Wechsel(p,a2);


 if (n==1) {
    a[gefunden][0]=a1;
    a[gefunden][1]=a2;
    gefunden++;
    return;
 }

 if (a2-a1<eps) {
    c=(a2-a1)/(n+1);
    b=a1;
    for (i=1;i<=n && gefunden<a.dimension();i++) {
        b+=c;
        a[gefunden][0]=a[gefunden][1]=b;
        gefunden++;
    }
    return;
 }

 b=(a1+a2)*0.5;
 d=(a2-a1)*0.5;

 while (p[0].eval(b)==0) {
        d*=0.5;
        b=a1+d;
 }

 if (Wechsel(p,a1)>Wechsel(p,b)) analysis::lokalisiere(p,a,a1,b,gefunden,eps);
 if (Wechsel(p,b)>Wechsel(p,a2)) analysis::lokalisiere(p,a,b,a2,gefunden,eps);
}


//  ---------------------------------------------------------------------------

int getrootdeg2(const double &s,const double &t,double &y)
{
 int reell=0;

 y=0.25*s*s+t;
 if (y<0) y*=-1;
 else reell=1;

 y=val::sqrt(y);
 return reell;
}


// divides f by polynomial x^2 -sx - t.
void divbyquadratfaktor(pol<double> &f,const double& s, const double& t, double eps)
{
    pol<double> g(1,2);

    g+=pol<double>(-s,1) + pol<double>(-t,0);

    f/=g;

    RoundPol(f,eps);
}



// Double Horner.
void dhorner(const val::vector<double> &a,val::vector<double> &b,const double& s,const double& t)
{
   if (a.dimension()<5) return;

   int i,grad=a.dimension()-1;
   if (b.dimension()<grad+1) b=val::vector<double>(0.0,grad+1);

   b[grad]= a[grad];
   b[grad-1]= s*a[grad]+a[grad-1];
   for (i= grad-2;i>0;i--) b[i]= a[i]+s*b[i+1]+t*b[i+2];
   b[0]= a[0]+t*b[2];
   b[grad-1]= b[grad-1]+s*b[grad];
   for (i=grad-2;i>2;i--) b[i]=b[i]+s*b[i+1]+t*b[i+2];
   b[2]= b[2]+t*b[4];
   return;
}


int bairstow(const pol<double> &f, double &s, double &t,int N,double E)
{
 int fertig=0,grad,m;
 double eps,delta,r;
 //term<double> *ppol;
 polIterator<double> ppol;

 grad=deg(f);
 if (grad<4) return 0;

 vector<double> a(0.0,grad+1),b(0.0,grad+1);


 for (ppol=f;ppol!=0;ppol++)                    // write f as vector 
     a[ppol.actualdegree()]=ppol.actualcoef();

 eps=delta=0.0;

 for (m=0;m<N;m++,s+=eps,t+=delta) {
	 dhorner(a,b,s,t);
	 eps= (b[0]*b[3]-b[1]*b[2])/(b[2]*(b[2]+s*b[3])-t*b[3]*b[3]);
	 delta= (t*b[1]*b[3]-b[0]*(b[2]+s*b[3]))/(b[2]*(b[2]+s*b[3])-t*b[3]*b[3]);
	 r= val::sqrt(eps*eps+delta*delta);
	 if (r<E) {
		 fertig=1;
		 break;
	 }
 }

 if (fertig) return m;
 else{
	 return 0;
 }
}


int iscomplexcyclomatic(const pol<double>& f)
{
    if (f.iszero()) return 0;

    int grad=deg(f);
    polIterator<double> pf;

    pf=f;

    if (grad<2 || (grad%2!=0)) return 0;
    pf++;
    if (!pf.actualvalid()) return 0;
    if (pf.actualdegree()!=0) return 0;
    if (leader(f)*f.getlastcoef()<0.0) return 0;
    return 1;
}



} // end namespace analysis

// ---------------------------------------------------------------------------------------------------

// Decomposition of a double value  wrt IEEE754 in sign, exponent and mantissa
void sign_exp_mantissa(const double &a,int& sign,int& exp,__int64& mantissa)
{
    int i;
    __int64 b,h(1),e(1);

    mantissa=0;exp=0;sign=0;

    b =  *((__int64*) &a);

    for (i=0;i<52;i++,h=h<<1) {
        if ((h&b) !=0) {mantissa|=h;}
    }

    for (;i<63;i++,h<<=1,e<<=1) {
        if ((h&b) !=0) exp|=e;
    }

    if ((h&b)!=0) sign=1;
    else sign =0;
}


int isNaN(const double& a)
{
    int exp,sign;
    __int64 mantissa;

    sign_exp_mantissa(a,sign,exp,mantissa);
    if (exp == 2047 && mantissa!=0) return 1;
    else return 0;
}


void normalize(const double& a,double& m,int& exp) // a = (2^exp) * m , mit 0.5 <= |m| < 1.
{
    int sign;
    __int64 mantissa,hexp;
    m=0.0;
    double *p;

    if (a==0.0) {
        exp=0;
        return;
    }

    sign_exp_mantissa(a,sign,exp,mantissa);

    if (exp == 2047) {
        m=a;
        return;
    }
    if (exp==0) {  // underflow
        return;
    }
    hexp= 1022;
    exp-=1022;
    hexp<<=52;
    mantissa|=hexp;

    p= (double*) (&mantissa);
    m=*p;
    if (sign==1) m=-m;
}



// --------------------------------------------------------------------------------------------------------------------------




double integral(const pol<double> &p,const double& a,const double& b)
{
    return p.integrate(a,b);
}

pol<double> realRoots(const pol<double> &g,val::vector<double> &Roots,const double& eps,const double &x1,const double &x2)
{
 pol<double> f;
 Roots=val::vector<double>();

 if (g.degree()<=0) return f;

 double A,b,a1,a2;
 val::vector<pol<double> > p;
 int N=0,grad,i,gefunden=0,z;

 f=g;
 

 f.normalize();
 

 if (x1!=0 && x1==x2) {
    A=f.eval(x1);
    if (abs(A)<eps) {
       Roots=val::vector<double>(1);
       Roots[0]=x1;
    }
    return f;
 }
 
 //std::cout<<f<<std::endl<<f.derive()<<std::endl<<f%f.derive()<<std::endl;

 f/=analysis::gcd(f,f.derive());

 grad=f.degree();
  

 if (grad==1) {
    A=-f[0]/f[1];
    if (x1!=x2 && (A<x1 || A>x2)) return f;
    else {
        Roots = val::vector<double>(A,1);
        return f;
    }
 }
 

 z = f.getlastdeg();
 f.getdivbypower(z);

 if (z>0 && x1<=0.0 && x2>=0.0) z=1;
 else z=0;
 
 analysis::Sturm(f,p);

 if (x1==0 && x2==0) {
    grad=deg(f);
    a2=a1=abs(f[0]);
    for (i=1;i<grad;i++) {
	a1+=abs(f[i]);
	a2=Max(a2,abs(f[i])+1);
    }
    a1=Max(a1,1.0);
    A=Min(a1,a2)+1; // f(A) * f(-A) != 0
    a1=-A;
    a2=A;
 }
 else {
    a1=x1;
    a2=x2;
    if (a1>a2) val::swap(a1,a2);
    if (f(a1)==0.0) a1-=eps;
    if (f(a2)==0.0) a2+=eps;
 }

 N=analysis::Wechsel(p,a1)-analysis::Wechsel(p,a2);
 
 
 Roots = val::vector<double>(0.0,N+z);
 val::vector<double[2]> a(N);

 analysis::lokalisiere(p,a,a1,a2,gefunden,eps);

 val::vector<int> computed(0,N);
 

 for (i=0;i<N;i++) {
      if (computed[i]) continue;
      a1=a[i][0];
      a2=a[i][1];
      if (f.eval(a1)==0) {
          Roots[i]=a1;
          continue;
      }
      if (f.eval(a2)==0) {
          if (analysis::Wechsel(p,a1)-analysis::Wechsel(p,a2+eps*0.5)==1) {
                Roots[i]=a2;
                if (i != N-1) a[i+1][0]=a2+eps*0.5;
                continue;
          }
          else {
              Roots[i+1]=a2; computed[i+1]=1;
              a2-=eps*0.5;
          }
      }
      while (a2-a1>=eps) {
            b=(a1+a2)/2;
            if (f.eval(b)==0) break;
	        if (f.eval(a1)*f.eval(b)<0) a2=b;
	        else a1=b;
      }
      b=(a1+a2)*0.5;
      //if ( (suche) || ((b>=x1) && (b<=x2)) )
      Roots[i]=b;
 }


 if (z) {
    Roots[N] = 0.0;
    f.getmultbypower(1);
 }

 return f;
}


int roots(const pol<double> &g,val::vector<double> &Real,val::vector<val::complex> &Comp,const double &eps)
{
    Real=val::vector<double>();
    Comp=val::vector<val::complex>();
    if (g.degree()<=0) return 0;

    pol<double> f=realRoots(g,Real,eps,0,0);
    int N=Real.dimension(),grad,i,ngefunden=0,fertig;
    double s,t,y;

    analysis::divbylinearfaktors(f,Real);
    grad=f.degree();
    N+=grad;
    if (grad<2) return N;

    grad/=2;
    Comp=val::vector<val::complex> (0.0,grad);

    for (i=grad;i>1;i--) {

        if (analysis::iscomplexcyclomatic(f)) {
            int i,n=deg(f),m=n/2;
            double b=f.getlastcoef()/leader(f),roh;
            b=exp(b,(double(1)/double(n))); // b = b^(1/n);
            for (i=0;i<m;i++) {
                roh=double(2*i + 1)*val::PI/double(n);
                Comp[ngefunden] = val::complex(b*val::cos(roh),b*val::sin(roh));
                ngefunden++;
            }
            return N;
        }

        fertig=0;
        s=1,t=1;
        fertig = analysis::bairstow(f,s,t,100,eps);
        if (fertig) {
            if (!analysis::getrootdeg2(s,t,y)) {
                Comp[ngefunden] = val::complex(0.5*s,y);
                ngefunden++;
            }
        }
        else {
            int d=deg(f);
            s=f[d-1]/f[d];
            if (s==0) s=1e-3;
            t=f[d-2]/f[d];
            if (t==0) t=1e-3;
            fertig= analysis::bairstow(f,s,t,100,eps);
            if (fertig) {
                if (!analysis::getrootdeg2(s,t,y)) {
                    Comp[ngefunden] = val::complex(0.5*s,y);
                    ngefunden++;
                }
            }
        }

        if (!fertig) {
            if (ngefunden!=Comp.dimension()) {
                if (ngefunden>0) {
                    val::vector<val::complex> hilf(ngefunden);
                    for (i=0;i<ngefunden;i++) hilf[i]=std::move(Comp[i]);
                    Comp=std::move(hilf);
                }
                else Comp=val::vector<val::complex>();
            }
            return N;
        }
        analysis::divbyquadratfaktor(f,s,t,eps);
    }
    if (deg(f)==2) {
        s=-f[1],t=-f[0];
        if (!analysis::getrootdeg2(s,t,y)) {
            Comp[ngefunden] = val::complex(0.5*s,y);
            ngefunden++;
        }
    }
    return N;
}

int extreme(const pol<double>& g,val::vector<double> &Maxima,val::vector<double> &Minima,const double &eps,const double &x1,const double& x2)
{
    Minima=Maxima=val::vector<double>();
    if (g.degree()<=1) return 0;

    pol<double> f=g.derive();
    vector<double> Roots;
    int i,k,N,NHP=0,NTP=0;

    realRoots(f,Roots,eps,x1,x2);
    N=Roots.dimension();

    if (N==0) return 0;

    vector<int> isextreme(0,N);

    for (i=0;i<N;i++) {
        if (f(Roots[i]-eps) < 0) {
            if (f(Roots[i]+eps) > 0) {isextreme[i]=-1;NTP++;}
        }
        else {
            if (f(Roots[i]+eps) < 0) {isextreme[i]=1;NHP++;}
        }
    }
    if (NHP) Maxima=val::vector<double>(NHP);
    for (i=k=0;i<N;i++) {
        if (isextreme[i]==1) {
            Maxima[k]=Roots[i];
            k++;
        }
    }
    if (NTP) Minima=val::vector<double>(NTP);
    for (i=k=0;i<N;i++) {
        if (isextreme[i]==-1) {
            Minima[k]=Roots[i];
            k++;
        }
    }
    return NHP+NTP;
}


//   -----------------------------------------------------------------------------------


double power(const double &a,int z)
{
    int n=abs(z);
    double y,x=1.0;

    if (a==0.0) {
        if (z==0.0) return 1.0;
        else if (z>0) return 0.0;
        else return val::Inf;
    }

    if (z==0) return 1.0;
    if (z>0) y=a;
    else y=1.0/a;

    while (n!=0) {
       if (n%2!=0) {
           x*=y;
       }
       y*=y;
       n=n/2;
    }
    return x;
}


double exp(const double &d)  // Expontial-Funktion: siehe Numerische Math: isc10simd.pdf
{
    if (d==val::Inf || d>=1e9) return val::Inf;
    if (d==-val::Inf || d<-1e9) return 0.0;
    if (isNaN(d)) return val::NaN;

    int q = int(d/analysis::MLN2),i,N=4;
    double s = d - q*analysis::L2U - q*analysis::L2L,t;

    s*=power(2.0,-N); // Reduction of argument

    // Evaluating Taylor series: first 8 summands of exp(x)-1
    t = ( ( ( s /40320 + 1.0/5040)*s + 1.0/720)* s + 1.0/120)* s ;
    t = ( ( ( ( t + 1.0/24)* s + 1.0/6)* s + 1.0/2)* s + 1)*s ;

    for ( i =0; i<N; i++)  t *= (2+t );

    return (t+1)*power(2.0,q);
}



double log(const double& d)   // Logarithmus naturalis: siehe Numerische Math: isc10simd.pdf
{
    int e,i;
    double m,x,y=0.0;

    if (d==0.0) return -val::Inf;  //-inf
    if (d<0.0) return val::NaN; //(-1.0/0.0 + 1.0/0.0);  //nan
    if (isNaN(d)) return val::NaN;
    if (d==val::Inf) return val::Inf;

    normalize(d,m,e);

    if (m<0.7071) {m*=2.0;e--;}

    x=(m-1)/(m+1);

    for (i=19;i>=1;i-=2) {
        y=y*x*x + 2.0/i;
    }

    return e*analysis::MLN2 + x*y;  // e*log(2) + x*y
}


double exp(const double &a,const double &x)
{
    if (a<0) return val::NaN;//(0.0/0.0);  //nan
    else if (a==0) {
        if (x==0) return 1.0;
        else if (x>0) return 0.0;
        else return val::NaN;//(0.0/0.0);
    }
    else return exp(log(a)*x);
}

double log(const double &a,const double &x)
{
    if (a<=0) return val::NaN;//(0.0/0.0);   // nan
    else if (x<0) return val::NaN; //(0.0/0.0);  //nan
    else if (x==0) return (-val::Inf);//(-1.0/0.0);  //-inf
    else return (log(x)/log(a));
}



// Trigonometric functions:

void xsincos0 (const double &t,double &scs,double &scc ) // 0  s < =4
{

    int i,N=3;
    double s=t;
    // Argument reduction
    s *= power(2.0,-N);
    s *= s ; // Evaluating Taylor series
    s = ( ( ( ( s /1814400 - 1.0/20160)* s + 1.0/360)* s - 1.0/12)* s + 1)* s ;
    for ( i =0; i<N; i++) s = (4-s ) * s ; // Applying double angle formula
    s /= 2.0 ;
    scs = sqrt ((2-s )* s ) ;
    scc = 1 - s ;
}


void xsincos(const double &d,double &scs,double &scc)
{
    double s = abs (d) ;
    int q = int( s * analysis::M_4_PI ) , r = q + ( q & 1 ) ;

    s -= r * analysis::PI4_A ;
    s -= r * analysis::PI4_B ;
    s -= r * analysis::PI4_C ;
    xsincos0 (s,scs,scc);
    if ( ( ( q + 1) & 2) != 0) { s = scc ; scc = scs ; scs = s ; }
    if ( ( ( q & 4) != 0) != (d < 0 ) ) scs = -scs;
    if ( ( ( q + 2) & 4) != 0) scc = -scc ;
}

double sin(const double &d)
{
    if (isNaN(d) || d==val::Inf || d==-val::Inf) return val::NaN;

    double scs,scc;

    xsincos(d,scs,scc);
    return scs;
}

double cos(const double &d)
{
    if (isNaN(d) || d==val::Inf || d==-val::Inf) return val::NaN;

    double scs,scc;

    xsincos(d,scs,scc);
    return scc;
}


double tan(const double &d)
{
    if (isNaN(d) || d==val::Inf || d==-val::Inf) return val::NaN;

    double scs,scc;

    xsincos(d,scs,scc);
    return scs/scc;
}



double arctan(const double &d)
{
    if (isNaN(d)) return val::NaN;
    if (d==val::Inf) return (val::PI/2);
    if (d==-val::Inf) return (-val::PI/2);
    if (d==0.0) return 0.0;
    double x,y,z = abs(d);


    if (z < 1.0) x = 1.0/z;
    else x = z ;
    int i ; // Applying cotangent half angle formula

    for (i=0 ;i<2; i++){
            x += sqrt(1+x*x);
    }
    x=1.0/x ;

    y=0.0; // Evaluating Taylor series
    for (i=10;i>=0;i--)  {
        y = y*x*x + power(-1,i)/double(2*i+1);
    }
    y *= x * power(2.0,2);
    if (z>1.0) y = PI/2 -y;
    return (d<0) ? -y : y ;
}

double arcsin(const double &x)
{
    if (x>1.0 || x<-1.0) return val::NaN;
    if (x==1.0) return PI/2;
    else if (x==-1.0) return -PI/2;
    return arctan(x/sqrt((1-x)*(1+x)));
}

double arccos(const double &x)
{
    if (x<-1.0 || x>1.0) return val::NaN;
    if (x==1.0) return 0.0;
    if (x==-1.0) return val::PI;
    if (x==0.0) return PI/2;
    else if (x<0) {
        double y=-x;
        return PI - arctan(sqrt((1-y)*(1+y))/y);
    }
    return arctan(sqrt((1-x)*(1+x))/x);
}

// hyperbolic functions:

double sinh(const double &x)
{
	return (exp(x) - exp(-x))*0.5;
}

double cosh(const double &x)
{
	return (exp(x) + exp(-x))*0.5;
}

double tanh(const double &x)
{
	double px = exp(x), nx = exp(-x);
	return (px -nx)/(px + nx);
}

double arsinh(const double &x)
{
	return log(x + sqrt(x*x +1));
}

double arcosh(const double &x)
{
	//if (x <= double(1)) return val::NaN;
	return log(x + sqrt(x*x -1));
}

double artanh(const double &x)
{
	return log((1+x)/(1-x))*0.5;
}


// complex functions

double atan2(const double& y, const double &x)
{
    if (x > 0) return val::arctan(y/x);
    else if (x < 0) {
        if (y >= 0) return val::arctan(y/x) + val::PI;
        else return  val::arctan(y/x) - val::PI;
    }
    else if (y > 0) return 0.5*val::PI;
    else if (y < 0) return -0.5*val::PI;
    else return val::NaN;
}

double arg(const val::complex &z)
{
return val::atan2(z.imaginary(), z.real());
}

val::complex exp(const val::complex &z)
{
    double r = val::exp(z.real());
    return val::complex(r * val::cos(z.imaginary()), r * val::sin(z.imaginary()));
}

complex sin(const complex &z)
{
    return complex(sin(z.real()) * cosh(z.imaginary()), cos(z.real()) * sinh(z.imaginary()));
}

complex cos(const complex &z)
{
    return complex(cos(z.real()) * cosh(z.imaginary()), -sin(z.real()) * sinh(z.imaginary()));
}

complex sinh(const complex &z)
{
    return complex(sinh(z.real()) * cos(z.imaginary()), cosh(z.real()) * sin(z.imaginary()));
}

complex cosh(const complex &z)
{
    return complex(cosh(z.real()) * cos(z.imaginary()), sinh(z.real()) * sin(z.imaginary()));
}

// ------------------------------------------------------------------------------------------------------------------


} // end namespace val

#endif // ANALYSIS_CPP
