
#include <val_basics.h>


namespace val
{

char primlistpath[] = "";




int gcd(int a,int b)
{
 if (a==0) return b;
 if (b==0) return a;

 int vorrest=a;

 if (a<b) {
	 a=b;
	 b=vorrest;
	 vorrest=a;
 }
 while (b != 0)
 {
       vorrest=b;
       b= a % b;
       a= vorrest;
 }

 return vorrest;
}

double getInf()
{
    double x=0.0;
    return (1.0/x);
}

double getNaN()
{
    double x=0.0;
    return (0.0/x);
}

const double Inf=getInf();
const double NaN=getNaN();

// function for computation of sqrt
namespace wurzel {

  inline double polynom(const double& x)
  {
     return double((10.0/11.0)*x + 161/880);
  }

  double potenz(const double& a,int z)   //int n
  {
    int n=val::abs(z);
    double y,x=1.0;

    if (a==0.0) {
        if (z==0.0) return 1.0;
        else if (z>0) return 0.0;
        else return val::Inf;//(1.0/0.0);
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

  void zerleg(const double& a,long& m,double& x)
  {
	  m=0;
	  do {
		  x= a/potenz(10,m);
		  if (x>=1.0) m+=2;
          if (x<0.01) m-=2;
	  }
      while (x<0.01 || x>=1);
  }

  double fehler(1e-8);

}// end namespace wurzel

double sqrt(const double& a)
{
  long int m;
  double alpha,beta,x;

  if (a<0.0) {
     return val::NaN; //nan
  }
  if (a==val::Inf) return val::Inf;
  if (a==0.0 || a==1.0) return a;
  wurzel::zerleg(a,m,x);
  beta= wurzel::polynom(x);
  alpha= x/beta;
  do {
     beta=(beta+alpha)*0.5;
     alpha= x/beta;
  }
  while ((beta-alpha)>=wurzel::fehler);
  return (wurzel::potenz(10.0,m/2)*beta);
}


double round(const double& x,int k)
{
    double y=val::abs(x),faktor=1,d=0,z=1,y1;
    unsigned limit=~0;
    int i;

    if (x==NaN || x==-NaN) return x;
    if (x==Inf || x==-Inf) return x;
    if (y > 1.0e30) return x;

    //limit = ~limit;
    //std::cout<<" "<<limit<<"  ";

    for(i=0;i<k;++i) z*=10.0;

    y*=z;

    //std::cout<<" "<<y<<"  ";

    while (faktor<y) faktor*=10;
    if (faktor>y) faktor/=10;
    //std::cout<<" "<<faktor<<"  ";
    if (y<double(limit)) {
        //std::cout<<" "<<limit<<"  ";
        d=unsigned(y);
        y-=d;
        if (y>=0.5) ++d;
        if (x<0) return -d/z;
        else return d/z;
    }
    y1=y;
    while (faktor>=1) {
        d+= double(unsigned(y1/faktor)) *faktor;
        y1-=double(unsigned(y1/faktor)) *faktor;
        faktor/=10;
    }
    //std::cout<<" "<<d<<"  ";
    y-=d;
    if (y>=0.5) ++d;
    if (x<0) return -d/z;
    else return d/z;

}



} // end namespace val
