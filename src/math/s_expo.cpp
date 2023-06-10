

#include <s_expo.h>


namespace val
{



int s_expo::dim=1;

int s_expo::ordtype=-2;

matrix<int> s_expo::ordmatrix;
const matrix<int> &s_expo::c_ordmatrix=s_expo::ordmatrix;

//
s_expo::s_expo(int c)
{
    coeff=new int[dim];
    for (int i=0;i<dim;i++) coeff[i]=c;
}

s_expo::s_expo(const s_expo& Y)
{
    coeff=new int[dim];
    for (int i=0;i<dim;i++) coeff[i]=Y.coeff[i];
}

s_expo::s_expo(s_expo&& Y)
{
    coeff=Y.coeff;
    Y.coeff=NULL;
}

//

const s_expo& s_expo::operator =(const s_expo &Y)
{
    for (int i=0;i<dim;i++) coeff[i]=Y.coeff[i];
    return *this;
}


const s_expo& s_expo::operator =(s_expo &&Y)
{
    if (coeff==Y.coeff) return *this;
    delete[] coeff;
    coeff=Y.coeff;
    Y.coeff=NULL;
    return *this;
}

//

int s_expo::operator==(const s_expo &Y) const
{
    for (int i=0;i<dim;i++) if (coeff[i]!=Y.coeff[i]) return 0;
    return 1;
}


int s_expo::operator <(const s_expo& x) const
{
 int i,j,m1,m2;

 m1=1;m2=0;
 if (!ordtype)  { // totdeg - lex
    m1=m2=0;
    for (i=0;i<dim;i++) {
		m1+=coeff[i];
		m2+=x.coeff[i];
    }

    if (m1<m2) return 1;
    if (m1>m2) return 0;
 }

 if (ordtype==-1 || m1==m2) {          // lex
    for (i=0;i<dim;i++)
		if (coeff[i] != x.coeff[i])
			return (coeff[i] < x.coeff[i]);
    return 0;
 }

 if (ordtype==-2) {    // DegRevLex
    m1=m2=0;
    for (i=0;i<dim;i++) {
		m1+=coeff[i];
		m2+=x.coeff[i];
    }
    if (m1<m2) return 1;
    if (m1>m2) return 0;

	for (i=dim-1;i>0;i--) {
		 if (coeff[i]>x.coeff[i]) return 1;
		 else if (coeff[i]<x.coeff[i]) return 0;
	 }
	 return 0;
 }

 if (ordtype!=-1 && !ordmatrix.isempty()) {
	 int y,z;

	 for (i=0;i<dim;i++) {
		 y=z=0;
		 for (j=0;j<dim;j++) {
			 y+=coeff[j]*s_expo::c_ordmatrix(i,j);
			 z+=x.coeff[j]*s_expo::c_ordmatrix(i,j);
		 }
		 if (y!=z) return (y<z);
	 }
	 return 0;
 }

 return 0;
}


int expocompare (const s_expo& a,const s_expo& b)
{
 int i,j,wa,wb;


 if (s_expo::ordtype==-1) {
	 for (i=0;i<s_expo::dim;i++) {
		 if (a.coeff[i]<b.coeff[i]) return -1;
		 else if(a.coeff[i]>b.coeff[i]) return 1;
	 }
	 return 0;
 }
 else if (s_expo::ordtype==0) {
	 i=0;
	 if (i==0) {
		 wa=wb=0;
		 for (j=0;j<s_expo::dim;j++) {
			 wa+=a.coeff[j];
			 wb+=b.coeff[j];
		 }
		 if (wa<wb) return -1;
		 else if (wa>wb) return 1;
		 i++;
	 }
	 for (;i<s_expo::dim;i++) {
		 if (a.coeff[i-1]<b.coeff[i-1]) return -1;
		 else if(a.coeff[i-1]>b.coeff[i-1]) return 1;
	 }
	 return 0;
 }
 else if (s_expo::ordtype==-2) {
      wa=wb=0;
      for (j=0;j<s_expo::dim;j++) {
			 wa+=a.coeff[j];
			 wb+=b.coeff[j];
      }
      if (wa<wb) return -1;
      else if (wa>wb) return 1;
	  for (i=s_expo::dim-1;i>0;i--) {
		 if (a.coeff[i]>b.coeff[i]) return -1;
		 else if (a.coeff[i]<b.coeff[i]) return 1;
	  }
	  return 0;
 }
 else {
	 for (j=0;j<s_expo::dim;j++) {
		 wa=wb=0;
		 for (i=0;i<s_expo::dim;i++) {
			 wa+=s_expo::c_ordmatrix(j,i)*a.coeff[i];
			 wb+=s_expo::c_ordmatrix(j,i)*b.coeff[i];
		 }
		 if (wa<wb) return -1;
		 else if (wa>wb) return 1;
	 }
	 return 0;
 }
}

//

int& s_expo::operator[] (int i)
{
    if (i<0 || i>=dim) Error::error("ERROR s_expo: index out of range!!");
    return coeff[i];
}

int s_expo::operator[] (int i) const
{
    if (i<0 || i>=dim) return 0;
    else return coeff[i];
}

//

s_expo s_expo::operator *(const s_expo& X) const
{
    s_expo Z;
    for (int i=0;i<dim;i++) Z.coeff[i] = coeff[i] + X.coeff[i];
    return Z;
}

const s_expo& s_expo::operator *=(const s_expo& X)
{
    for (int i=0;i<dim;i++) coeff[i] += X.coeff[i];
    return *this;
}

int s_expo::operator |(const s_expo &X) const
{
    for (int i=0;i<dim;i++) if (coeff[i]>X.coeff[i]) return 0;
    return 1;
}

s_expo s_expo::operator /(const s_expo &X) const
{
    s_expo Z;
    for (int i=0;i<dim;i++) Z.coeff[i] = coeff[i] - X.coeff[i];
    return Z;
}

const s_expo& s_expo::operator /=(const s_expo &X)
{
    for (int i=0;i<dim;i++) coeff[i]-=X.coeff[i];
    return *this;
}

s_expo gcd(const s_expo& X,const s_expo& Y)
{
    s_expo Z;
    for (int i=0;i<s_expo::dim;i++) Z.coeff[i] = val::Min(X.coeff[i],Y.coeff[i]);
    return Z;
}


s_expo lcm(const s_expo& X,const s_expo& Y)
{
    s_expo Z;
    for (int i=0;i<s_expo::dim;i++) Z.coeff[i] = val::Max(X.coeff[i],Y.coeff[i]);
    return Z;
}

//

int s_expo::totdeg() const
{
    int d=0;
    for (int i=0;i<dim;i++) d+=coeff[i];
    return d;
}

//


std::istream& operator >>(std::istream& is,s_expo& x)
{
   for (int i=0;i<s_expo::dim;i++) {
	   x.coeff[i]=0;
	   is>>x.coeff[i];
   }
   return is;
}


std::ostream& operator <<(std::ostream& os,const s_expo& x)
{
  for (int i=0;i<s_expo::dim-1;i++)  os<<x.coeff[i]<<' ';
  os<<x.coeff[s_expo::dim-1];
  return os;
}



} // end namespace
