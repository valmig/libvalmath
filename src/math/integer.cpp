#ifndef INTEGER_CPP
#define INTEGER_CPP

#include <integer.h>
#include <error.h>
#include <val_basics.h>

#define __int64 long long


namespace val
{


namespace hilfinteger
{
 int numberofbits=8*sizeof(unsigned int);
 unsigned int highon=1<<(numberofbits-1),      // highon : highest bit is 1, the rest 0
              lowon=1;                         // lowon : lowest bit is 1, the rest 0
}


void hilfinteger::multunsigned(unsigned a,unsigned b,unsigned &high,unsigned &low)
{
 __int64 c;

 c=(__int64) (a) * (__int64)(b);
 high=unsigned(c>>32);
 low=unsigned(c);
}


// Condition : a1<d a=a1*2^32+a2
unsigned hilfinteger::divrest64(unsigned a1,unsigned a2,unsigned d1)
{
 int neg=0;
 __int64 a,q,d,r;

 a=((__int64) (a1)<<32) + (__int64)(a2);
 d=(__int64)(d1);
 if (a<0) {
     neg=1;
     a=(d<<32)-a;
 }
 q=a/d;
 if (neg) {
     r=a%d;
     q=((__int64)(1)<<32)-q;
     if (r) {
         q--;
     }
 }
 return unsigned(q);
}


// Prem: a=a2a1a0, b=b1b0 to Basis 2^numberofbits, a2<b1
unsigned hilfinteger::trialdiv(unsigned *a,unsigned *b)
{
 int groesser=0,i,anz=0;
 unsigned q;
 unsigned trialproduct[3];

 q=divrest64(a[2],a[1],b[1]);
 do {
     //Multiplication from a with q:
     multbyunsigned(b,2,q,trialproduct);
     for (i=2;i>=0;i--)
         if (trialproduct[i]!=a[i]) break;

     if (i==-1) groesser=0;
     else if (trialproduct[i]>a[i]) groesser=1;
     else groesser=0;

     if (groesser) {
         anz++;
         q--;
         if (anz>2) {
             Error::warning("\nWARNING: hilfinteger::trialdiv: more than 1 correction!");
         }
     }
 }
 while (groesser);
 return q;
}


unsigned hilfinteger::modularinvers(unsigned a0)
{

 // Binary computation by Hensel 
    unsigned alpha,zweih_i=1,s=1;
    int i;

    for (i=1;i<numberofbits;i++) {
        zweih_i<<=1;  // zweih_i = 2^i
        alpha = a0 * s;
        if (alpha & zweih_i) s|= zweih_i;  // ih i-th Bit ==1 , so s += 2^i
    }

    return s;
}


// Condition: a divides b!
integer hilfinteger::exactdivision(integer b,integer a)
{
 if (a.dat==NULL)  {Error::error("\nERROR: hilfinteger::exactdivision: Division by zero!");}
 if (b.dat==NULL) return integer();

 int s,k,lq,la,lb; 
 unsigned a1,d=1;
 integer q,y,x,c;

 for (s=0;!a.dat[s];s++) ;
 if (s) {
     a.truncby(s);
     b.truncby(s);
 }
 for (k=0;!(a.dat[0]&d);k++) d=d<<1;
 if (k) {
     a.shiftrightby(k);
     b.shiftrightby(k);
 }
 // a is odd now!
 a1=hilfinteger::modularinvers(a.dat[0]);
 //std::cout<<"\na1= "<<a1;

 // Get length of q:
 la=hilfinteger::abs(a.laenge);
 lb=hilfinteger::abs(b.laenge);
 lq=lb-la;
 y.dat=b.dat;
 y.dat+=lq;
 y.laenge=la;
 q.laenge=lq;
 if (!y.abslower(a)) q.laenge++;
 q.dat=new unsigned[q.laenge];
 lq++;
 //
 y.dat=b.dat;
 y.laenge=lq;
 x.dat=a.dat;
 c.dat=new unsigned[hilfinteger::min(la,lq)+1];

 for (k=0;k<q.laenge;k++) {
     q.dat[k]=a1*y.dat[0];
     x.laenge=hilfinteger::min(la,lq-k);
     c.laenge=hilfinteger::multbyunsigned(x.dat,x.laenge,q.dat[k],c.dat);
     y.subto(c);
     y.dat++;
     y.laenge--;
 }
 y.dat=NULL;
 x.dat=NULL;
 q.laenge*=signum(a.laenge)*signum(b.laenge);
 return q;
}


unsigned hilfinteger::ggT(unsigned a,unsigned b)
{
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


// Division with rem. a=q*b+r. Condition: q=r=0! and abs(b.laenge)=1
integer hilfinteger::restsingle(const integer &a,const integer& b)
{
 int la=hilfinteger::abs(a.laenge);

 if (!la) return integer();
 if (la==1) {
     return integer(a.dat[0]%b.dat[0]);
 }

 unsigned wert=0,q;
 integer r;

 if (a.dat[la-1]<b.dat[0]) {
     la--;
     wert=a.dat[la];
 }
 for (int i=la-1;i>=0;i--) {
     if (wert) {
         q=hilfinteger::divrest64(wert,a.dat[i],b.dat[0]);
     }
     else q=a.dat[i]/b.dat[0];
     wert=a.dat[i]-q*b.dat[0];
 }
 if (wert) {
     r.dat=new unsigned[1];
     r.dat[0]=wert;
     r.laenge=1;
 }
 return r;
}


// =============== Element-Functions of class integer ========================
// =========================================================================

int integer::mem_optimized=1;
int integer::Maxlength=0;
integer::Output_Type integer::Output_Style=INTEGER;
int integer::Output_Digits=16;



// ================= Constructors and destructor =========================


integer:: integer(int n)
{

 if (n==0) {
     dat=NULL;
     laenge=0;
 }
 else {
     dat=new unsigned int[1];
     dat[0]=hilfinteger::abs(n);
     laenge=hilfinteger::signum(n);
     if (!Maxlength) Maxlength=1;
 }
 return;
}


integer:: integer(unsigned n)
{
 if (!n) {
     dat=NULL;
     laenge=0;
 }
 else {
     dat=new unsigned int[1];
     dat[0]=n;
     laenge=1;

     if (!Maxlength) Maxlength=1;
 }
 return;
}



integer:: integer(const integer& x)
{

 if (x.dat==NULL) {
     dat=NULL;
     laenge=0;
     return;
 }

 int i,l=hilfinteger::abs(x.laenge);

 laenge=x.laenge;
 dat= new unsigned int[l];
 for (i=0;i<l;i++) dat[i]=x.dat[i];

 if (l>Maxlength) Maxlength=l;
}


//std::move-constructor
integer::integer(integer&& x)
{
    dat=x.dat;
    laenge=x.laenge;
    x.dat=NULL;x.laenge=0;
}

integer::integer(const integer& x,int n)
{
 if (x.dat==NULL) {
     dat=NULL;
     laenge=0;
     return;
 }

 int i,l=hilfinteger::abs(x.laenge);

 if (n<0) n=0;

 laenge=l+n;
 dat = new unsigned int[laenge];
 for (i=0;i<n;i++) dat[i]=0;
 for (;i<laenge;i++) dat[i]=x.dat[i-n];

 if (laenge>Maxlength) Maxlength=laenge;
 if (hilfinteger::signum(x.laenge)==-1) laenge=-laenge;
}


// =========================================================================


// ======================== Cast-operators ================================


integer::operator int() const
{
 if (dat==NULL) return 0;
 if (laenge<0) return -int(dat[0]);
 else return dat[0];
}


integer::operator unsigned() const
{
 if (dat==NULL) return 0;
 //if (laenge<0) return -int(dat[0]);
 else return dat[0];
}

integer::operator double() const
{
    double b(0.0);
    int l = hilfinteger::abs(laenge),exponent = 32*(l)-1,k,i,j;
    unsigned bit = 1<<31;
    __int64 *bint= (__int64 *) &b,Exp(0);


    if (l==0) return b;
    if (l==1) {
        b= double(dat[0]);
        if (laenge<0) b=-b;
        return b;
    }
    if (l>=33) {
        if (laenge > 0) return val::Inf;//(1.0/0.0); //inf
        else return -val::Inf;//(-1.0/0.0); //-inf
    }

    // Find first left bit != 0
    for (k=l-1;k>=0;k--,exponent-=32)
        if (dat[k]) break;
    for (i=31;i>=0;i--,exponent--,bit>>=1) {
        if (dat[k] & bit) break;
    }

    for (j=52;j>=0;) {
         j--;i--;
         bit>>=1;
         if (i==-1) {
            i=31;
            k--;
            bit=1<<31;
         }
         if (k==-1) break;
         if (dat[k] & bit) (*bint)=(*bint) | 1;
         if (j>=0) (*bint) <<=1;
    }
     (*bint)>>=1;
    j++;

    if (exponent >1024) {
        if (laenge > 0) return val::Inf;//(1.0/0.0); //inf
        else return -val::Inf;//(-1.0/0.0); //-inf
    }

    (*bint)<<=j;
    if (laenge<0) Exp++;
    Exp<<=11;
    Exp=Exp | ((__int64) 1023 +  (__int64) exponent);
    Exp<<=52;
    (*bint) = Exp | (*bint);
    return b;
}


// =================== Assignment operators ===============================

integer& integer::operator =(const integer& x)
{
 if (dat==x.dat) return *this;
 int i,l=hilfinteger::abs(x.laenge);

 if (x.dat==NULL) {
     delete[] dat;
     dat=NULL;
     laenge=0;
     return *this;
 }

 if (hilfinteger::abs(laenge)!=l) {
     delete[] dat;
     //laenge=x.laenge;
     dat= new unsigned int[l];
 }
 laenge=x.laenge;
 for (i=0;i<l;i++) dat[i]=x.dat[i];
 return *this;
}


integer& integer::operator=(integer&& x)
{
    if (dat==x.dat) return *this;
    if (dat!=NULL) delete[] dat;
    dat=x.dat;
    laenge=x.laenge;
    x.dat=NULL;x.laenge=0;
    return *this;
}


integer& integer::operator =(int n)
{

 if (n==0) {
     delete[] dat;
     dat=NULL;
     laenge=0;
     return *this;
 }
 else if (laenge!=1 && laenge!=-1){
     delete[] dat;
     dat=new unsigned int[1];
 }
 dat[0]=hilfinteger::abs(n);
 laenge=hilfinteger::signum(n);

 if (!Maxlength) Maxlength=1;
 return *this;
}

// ==========================================================================


// ====================== sonst. Memberfktn =================================


int integer::iseven() const
{
 if (dat==NULL) return 1;
 if (dat[0]&1) return 0;
 else return 1;
}


int integer::isodd() const
{
 if (dat==NULL) return 0;
 if (dat[0]&1) return 1;
 else return 0;
}

// ===========================================================================


// ===================== comparisons ==============================

int integer::operator ==(const integer& x) const
{
 if (laenge!=x.laenge) return 0;

 int i,l=hilfinteger::abs(laenge);

 for (i=0;i<l;i++)
     if (dat[i]!=x.dat[i]) return 0;
 return 1;
}


int integer::operator !=(const integer& x) const
{
 return !( (*this)==x );
}


int integer::operator ==(int n) const
{
 if (laenge>1 || laenge<-1) return 0;
 if (dat==NULL) return (!n);
 else return (dat[0]==unsigned(hilfinteger::abs(n)) && laenge==hilfinteger::signum(n));
}



int integer::operator !=(int n) const
{
 return !((*this)==n);
}


int integer::operator <(const integer& x) const
{
 int i,sign=signum(),signx=x.signum(),l=hilfinteger::abs(laenge);


 if (sign!=signx) return 1+signx;

 // now sign==signx
 if (l<hilfinteger::abs(x.laenge)) return 1+sign;
 if (l>hilfinteger::abs(x.laenge)) return 1-sign;

 // now sign==signx a. laenge==x.laenge
 if (laenge==0) return 0;

 for (i=l-1;i>=0;i--)
     if (dat[i]!=x.dat[i]) break;
 if (i==-1) return 0;
 return (dat[i]<x.dat[i])? 1+sign : 1-sign;
}


int integer::operator <=(const integer& x) const
{
 int i,sign = signum(),signx=x.signum(),l=hilfinteger::abs(laenge);

 if (sign!=signx) return 1+signx;

 // now sign==signx
 if (l<hilfinteger::abs(x.laenge)) return 1+sign;
 if (l>hilfinteger::abs(x.laenge)) return 1-sign;

 // now sign==signx a. abs(laenge)==abs(x.laenge)
 if (laenge==0) return 1;

 for (i=l-1;i>=0;i--)
     if (dat[i]!=x.dat[i]) break;
 if (i==-1) return 1;
 return (dat[i]<x.dat[i])? 1+sign : 1-sign;
}


int integer::operator >(const integer& x) const
{
 return !(*this<=x);
}


int integer::operator >=(const integer& x) const
{
 return !(*this<x);
}


int integer:: operator <(int n) const
{
 integer x(n);
 return (*this<x);
}


int integer:: operator <=(int n) const
{
 integer x(n);
 return (*this<=x);
}


int integer:: operator >(int n) const
{
 integer x(n);
 return (*this>x);
}


int integer:: operator >=(int n) const
{
 integer x(n);
 return (*this>=x);
}


// returns: abs(*this)<y 
int integer::abslower(const integer &y) const
{
 int l=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);

 if (l<ly) return 1;
 else if (ly<l) return 0;

 for (int i=l-1;i>=0;i--) {
     if (dat[i]<y.dat[i]) return 1;
     else if (dat[i]>y.dat[i]) return 0;
 }

 return 0;
}




// ==========================================================================



// ========================== Bit-operators ================================

// multiplication with 2^n
integer integer::operator <<(int n) const
{
 if (!n) return *this;
 if (dat==NULL) return integer();
 if (n<0) return operator>>(-n);

 int k=n/hilfinteger::numberofbits,r=n%hilfinteger::numberofbits;

 if (!r) return integer(*this,k);

 int i,s=hilfinteger::numberofbits-r,l=hilfinteger::abs(laenge);//lx;
 unsigned int wert;
 integer x;

 i=l-1;
 wert=dat[i]>>s;
 l+=k;
 if (wert) {
     l++;
     x.dat=new unsigned int[l];
     x.dat[l-1]=wert;
 }
 else x.dat= new unsigned int[l];

 for (;i>0;i--) {
     wert=dat[i-1]>>s;
     x.dat[i+k] = (dat[i]<<r) | wert;
 }
 x.dat[k]=dat[i]<<r;
 for (i=k-1;i>=0;i--) x.dat[i]=0;

 if (l>Maxlength) Maxlength=l;
 if (laenge<0) x.laenge=-l;
 else x.laenge=l;

 return x;
}


// integer division by 2^n
integer integer::operator >>(int n) const
{

 if (!n) return *this;
 if (dat==NULL) return integer();
 if (n<0) return operator<<(-n);

 int k = n/hilfinteger::numberofbits,l=hilfinteger::abs(laenge);

 if (k>=l) return integer();

 int r=n%hilfinteger::numberofbits,lx=l-k;
 unsigned int wert,last;
 integer x;

 if (!r) {
     x.dat=new unsigned int[lx];
     for (int i=0;i<lx;i++) x.dat[i]=dat[i+k];
     if (laenge<0) x.laenge=-l;
     else x.laenge=l;
     return x;
 }

 int s=hilfinteger::numberofbits-r;

 last=dat[l-1]>>r;

 if (last) x.laenge=lx;
 else x.laenge=lx-1;

 if (!x.laenge) return x;
 x.dat=new unsigned int[x.laenge];

 for (int i=0;i<lx-1;i++) {
     wert=dat[i+k+1]<<s;
     x.dat[i]=wert | (dat[i+k]>>r);
 }
 if (last) x.dat[lx-1]=last;

 if (laenge<0) x.laenge=-x.laenge;

 return x;
}


// =======================================================================================

// ====================== private - Addition und Substraktion ===========================


// Prem.: abs(laenge)>=abs(y.laenge), z.dat=NULL. Result: z=signum(*this)*(abs(*this)+abs(y))
void integer::add(const integer& y,integer &z) const
{
 // Triviale Faelle:
 if (dat==NULL) return;
 if (y.dat==NULL) {
     z=*this;
     return;
 }

 int i,l=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge),L;
 unsigned int *feld,*pz,*px,*py,wert=0;

 L=l;
 if (l==ly) L++;

 feld=new unsigned int[L];
 pz = feld; px=dat; py=y.dat;

 for (i=0;i<ly;++i,px++,py++,pz++) {
     *pz= *px + wert;
     if (!(*pz) && wert) {
         *pz = *py;
         wert=1;
     }
     else if (( *pz += *py )< *py ) wert=1;
     else wert=0;
 }
 for (;i<l;++i,px++,pz++) {
     *pz = *px + wert;
     if (!(*pz) && wert) wert=1;
     else wert=0;
 }

 if (wert) {
     if (l<L) {
        *pz=wert;
        z.dat=feld;
        z.laenge=L;
     }
     else {
        L++;
        z.dat= new unsigned[L];
        z.laenge=L;
        for (i=0;i<l;++i) z.dat[i]=feld[i];
        z.dat[l] =wert;
        delete[] feld;
     }
 }
 else {
    z.dat=feld;
    z.laenge=l;
 }

 feld=NULL;
 if (z.laenge>Maxlength) Maxlength=z.laenge;

 if (laenge<0) z.laenge*=-1;

 return;
}



// Prem.: abs(laenge)>=abs(y.laenge)
integer& integer::addto(const integer &y)
{
    if (dat==NULL || y.dat==NULL) return *this;

    int i,l=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
    unsigned wert=0,*px=dat,*py=y.dat;

    for (i=0;i<ly;++i,++px,++py) {
        *px += wert;
        if (!(*px) && wert) {
            *px += *py ;
            wert=1;
        }
        else if ((*px += *py )< *py) wert=1;
        else wert=0;
    }
    for (;i<l && wert;++i) {
        *px +=wert;
        if (!(*px) && wert) wert=1;
        else wert=0;
    }
    if (wert) {
        l++;
        unsigned *feld=new unsigned[l];
        for (i=0;i<l-1;++i) feld[i]=dat[i];
        feld[l-1]=wert;
        delete[] dat;
        dat = feld;
        feld=NULL;
        if (laenge<0) laenge = -l;
        else laenge =l;
        if (l>Maxlength) Maxlength=l;
    }
    return *this;
}



// Pre.: abs(*this)>=abs(y), z.dat=NULL. Result: z=signum(*this)*(abs(*this)-abs(y))
void integer::sub(const integer& y,integer &z) const
{
 // Trivial cases:
 if (dat==NULL) return;  // => y.dat==NULL
 if (y.dat==NULL) {
     z=*this;
     return;
 }

 int i,l=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 unsigned int *feld,wert=0,*pz,*px,*py;

 feld=new unsigned int[l];

 px=dat; py=y.dat; pz=feld;

 for (i=0;i<ly;++i,++px,++py,++pz) {
     wert+= *py;
     if ((*py) && !wert) {
         wert=1;
         *pz = *px;
     }
     else {
         *pz = *px -wert;
         if (*pz > *px) wert=1;
         else wert=0;
     }
 }
 for (;i<l;++i,++px,++pz) {
     *pz = *px -wert;
     if (*pz > *px ) wert=1;
     else wert=0;
 }

 if (feld[l-1]) {
     z.dat=feld;
     feld=NULL;
 }
 else {
     for (i=l-2;i>=0;--i)
         if (feld[i]) break;
     if (i==-1) {
         delete[] feld;
         return ;
     }
     l=i+1;
     if (mem_optimized) {
        z.dat=new unsigned int[l];
        for (i=0;i<l;++i) z.dat[i]=feld[i];
        delete[] feld;
     }
     else {
        z.dat=feld;
        feld=NULL;
     }
 }
 if (laenge<0) z.laenge=-l;
 else z.laenge=l;

 return ;
}



//Prem.: *this >=0 *this>=y. Afterwards *this=*this - |y| , keeping zeros at the end.
void integer::subto(const integer& y)
{

 // Triviale Faelle:
 if (dat==NULL) return;  // => y.dat==NULL
 if (y.dat==NULL) {
     return;
 }

 int i,l=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 unsigned int feld=0,wert=0;


 for (i=0;i<ly;++i) {
     wert+=y.dat[i];
     if (y.dat[i] && !wert) {
         wert=1;
         feld=dat[i];
     }
     else {
         feld=dat[i]-wert;
         if (feld>dat[i]) wert=1;
         else wert=0;
     }
     dat[i]=feld;
 }
 for (;i<l;++i) {
     feld=dat[i]-wert;
     if (feld>dat[i]) wert=1;
     else {dat[i]=feld;break;}//wert=0;
     dat[i]=feld;
 }

 if (dat[l-1]) {
     return;
 }
 else {
     for (i=l-2;i>=0;--i)
         if (dat[i]) break;
     l=i+1;
 }
 laenge=l;

 return ;
}


// ======================================================================================================


// ========================== Additive operators =====================================


integer integer::operator +(const integer &y) const
{
 if (y.dat==NULL) return *this;
 if (dat==NULL) return y;

 int signx,signy,lx=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 integer z;

 signx = (laenge<0)? -1:1;
 signy = (y.laenge<0)? -1:1;

 if (signx==signy) {
     if (lx<ly) y.add(*this,z);
     else add(y,z);
 }
 else {
     if (abslower(y)) {  //(this<y) {
         y.sub(*this,z);
     }
     else {
         sub(y,z);
     }
 }
 return z;
}


integer integer::operator -(const integer &y) const
{
 if (y.dat==NULL) return *this;
 if (dat==NULL) {
     integer z(y);
     z.laenge=-y.laenge;
     return z;
 }

 int signx,signy,lx=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 integer z;

 signx = (laenge<0)? -1:1;
 signy = (y.laenge<0)? -1:1;

 if (signx!=signy) {
     if (lx<ly) {
         y.add(*this,z);
         z.laenge=-z.laenge;
     }
     else add(y,z);
 }
 else {
     if (abslower(y)) {
         y.sub(*this,z);
         z.laenge=-z.laenge;
     }
     else sub(y,z);
 }

 return z;
}


integer integer::operator -() const
{
 integer z(*this);

 z.laenge=-z.laenge;
 return z;
}


integer& integer::operator +=(const integer& y)
{
 if (y.dat==NULL) return *this;
 if (dat==NULL) {
     *this=y;
     return *this;
 }

 int signx,signy,lx=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 unsigned *ph;
 integer z;


 signx = (laenge<0)? -1:1;
 signy = (y.laenge<0)? -1:1;

 if (signx==signy) {
     if (lx<ly) y.add(*this,z);
     else return addto(y);
 }
 else {
     if (abslower(y)) {                     //this<y
         y.sub(*this,z);
     }
     else {

         subto(y);
         if (laenge==0) {
            delete[] dat;
            dat=NULL;
            return *this;
         }
         if (signx<0 && laenge>0) laenge=-laenge;
         return *this;
     }
 }
 ph=dat;
 dat=z.dat;
 laenge=z.laenge;
 z.dat=ph;

 return *this;
}


integer& integer::operator -=(const integer& y)
{
 if (y.dat==NULL) return *this;
 if (dat==NULL) {
     *this=y;
     laenge=-y.laenge;
     return *this;
 }

 int signx,signy,lx=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge);
 unsigned *ph;
 integer z;

 signx = (laenge<0)? -1:1;
 signy = (y.laenge<0)? -1:1;

 if (signx!=signy) {
     if (lx<ly) {
         y.add(*this,z);
         z.laenge=-z.laenge;
     }
     else return addto(y);
 }
 else {
     if (abslower(y)) {
         y.sub(*this,z);
         z.laenge=-z.laenge;
     }
     else sub(y,z);
 }
 ph=dat;
 dat=z.dat;
 laenge=z.laenge;
 z.dat=ph;
 return *this;
}



integer& integer::operator++() // praefix ++
{
 *this=*this+integer(1);
 return *this;
}


integer& integer::operator++(int) // postfix ++
{
 *this = *this + integer(1);
 return *this;
}


integer& integer::operator --() //praefix --
{
 *this=*this-integer(1);
 return *this;
}


integer& integer::operator --(int)  // postfix --
{
 *this=*this - integer(1);
 return *this;
}

// ==========================================================================



// ============= multiprecission - Mulltiplikation =====================================



// Result: z = abs(*this) * abs (x) . Prem.: c has right dimension. Returns length.
// length of c is at least l+1.
int hilfinteger::multbyunsigned(unsigned *dat,int l,unsigned x,unsigned *c)
{
 unsigned high,low,wert=0;
 int i;

 for (i=0;i<l;++i,++c) {
     hilfinteger::multunsigned(dat[i],x,high,low);
     *c =wert+low;
     if (*c < wert) wert=1;
     else wert=0;
     wert+=high;
 }
 if (wert) {
     *c = wert;
     l++;
 }
 else *c = 0;

 return l;
}



// ================================================================================


// ============================== multiplicative functions ======================================



// Prem.: z=0
void integer::classicmult(const integer &y,integer &z) const
{
 if (dat==NULL || y.dat==NULL) return;

 int i,j,lx=hilfinteger::abs(laenge),ly=hilfinteger::abs(y.laenge),signx,signy,la,lb,l;//lh;
 unsigned wert,*c,*h=NULL,*a,*b,*ph;

 l=lx+ly;

 if (lx<ly) {
     la=ly;lb=lx;
     a=y.dat;b=dat;
 }
 else {
     la=lx;lb=ly;
     a=dat;b=y.dat;
 }
 c=new unsigned int[l];
 c[l-1]=0;

 if (lb>1) {
     h=new unsigned int[la+1];
 }
 // 1.:
 hilfinteger::multbyunsigned(a,la,b[0],c);

 for (i=1;i<lb;++i) {
     hilfinteger::multbyunsigned(a,la,b[i],h);
     // add now:
     wert=0;
     ph=h;
     for (j=0;j<la;++j,++ph) {
         *ph += wert;
         c[j+i] += *ph;
         if (wert && !(*ph)) wert=1;
         else if (c[j+i] < *ph) wert=1;
         else wert=0;
     }
     c[i+la] = wert + *ph;
 }
 if (lb>1) delete[] h;

 if (!c[l-1]) l--;
 z.laenge=l;
 z.dat=c;
 c=NULL;
 //

 if (l>Maxlength) Maxlength=l;
 signx = (laenge<0)? -1:1;
 signy = (y.laenge<0)? -1:1;
 z.laenge*=signx*signy;

 return;
}


integer shift(const integer &x,int n)
{
 if (x.dat==NULL) return integer(); // =0

 int lx=hilfinteger::abs(x.laenge),l,i;
 integer z;

 l=lx+n;
 z.dat=new unsigned int[l];
 if (x.laenge<0) z.laenge=-l;
 else z.laenge=l;
 for (i=0;i<n;i++) z.dat[i]=0;
 for (i=0;i<lx;i++) z.dat[i+n]=x.dat[i];

 if (l>integer::Maxlength) integer::Maxlength=l;
 return z;
}


integer trunc(const integer &x,int n)
{
 int lx=hilfinteger::abs(x.laenge),l,i;
 integer z;

 if (lx<=n) return z;

 l=lx-n;
 z.dat=new unsigned int[l];
 for (i=0;i<l;i++) z.dat[i]=x.dat[i+n];
 if (x.laenge<0) z.laenge=-l;
 else z.laenge=l;
 return z;
}


integer rest(const integer& x,int n) // Rule: x = shift(trunc(x,n),n) + rest(x,n)
{
 if (n<=0) return integer(); //= 0

 int lx=hilfinteger::abs(x.laenge),l,i;
 integer z;

 l=hilfinteger::min(lx,n);
 for (i=l-1;i>=0;i--) if (x.dat[i]!=0) break;
 if (i==-1) return z;
 l=i+1;
 z.dat=new unsigned int[l];
 for (;i>=0;i--) z.dat[i]=x.dat[i];
 if (x.laenge<0) z.laenge=-l;
 else z.laenge=l;
 return z;
}


/*
#ifdef STDTHREAD
void hilfinteger::karatsuba(const integer &u, const integer &v,integer &result)
{
 int s,n,lx=hilfinteger::abs(u.length()),lv=hilfinteger::abs(v.length());

 n=hilfinteger::max(lx,lv);
 if (n<32) {
     u.classicmult(v,result);
     return ;
 }
 else { //karatsuba:
     s=n/2;
     integer w(trunc(u,s)), //w=trunc(*this,s);
             x(rest(u,s)),  //x=rest(*this,s);
             y(trunc(v,s)),     //y=trunc(v,s);
             z(rest(v,s)),      //z=rest(v,s);
             r,p,q;
             //r((w+x)*(y+z)), //r=(w+x)*(y+z);
             //p(w*y),         //p=w*y;
             //q(x*z);

     karatsuba(w+x,y+z,r);
     karatsuba(w,y,p);
     karatsuba(x,z,q);
     result = (shift(p,2*s)+shift(r-p-q,s)+q);
 }
}
#endif
*/


integer integer::operator *(const integer &v) const
{
 int n,lx=hilfinteger::abs(laenge),lv=hilfinteger::abs(v.laenge);

 n=hilfinteger::min(lx,lv);
 if (n<32) {
     integer z;
     classicmult(v,z);
     return z;
 }
 else { //karatsuba:
     int s=n/2;
     integer w(trunc(*this,s)), //w=trunc(*this,s);
             x(rest(*this,s)),  //x=rest(*this,s);
             y(trunc(v,s)),     //y=trunc(v,s);
             z(rest(v,s)),      //z=rest(v,s);
             r((w+x)*(y+z)), //r=(w+x)*(y+z);
             p(w*y),         //p=w*y;
             q(x*z);

     return(shift(p,2*s)+shift(r-p-q,s)+q);
 }
}


integer& integer::operator *=(const integer& a)
{
 int n,l=hilfinteger::abs(laenge),la=hilfinteger::abs(a.laenge);

 n=hilfinteger::max(l,la);
 if (n<32) {
     unsigned *ph;
     integer z;

     classicmult(a,z);
     ph=dat;
     dat=z.dat;
     laenge=z.laenge;
     z.dat=ph;
     return *this;
 }
 else{
     *this = *this * a;
     return *this;
 }
}


// Division with remainder. a=q*b+r. Prem.: q=r=0! and abs(b.laenge)=1
void integer::divrestsingle(const integer &a,const integer& b,integer &q,integer &r)
{
 int la=hilfinteger::abs(a.laenge);

 if (!la) return;
 if (la==1) {
     q=integer(a.dat[0]/b.dat[0]);
     r=integer(a.dat[0]%b.dat[0]);
     return;
 }

 int lq;
 unsigned wert=0;


 if (a.dat[la-1]<b.dat[0]) {
     lq=la-1;
     la--;
     wert=a.dat[la];
 }
 else lq=la;
 q.dat= new unsigned[lq];
 q.laenge=lq;
 for (int i=la-1;i>=0;i--) {
     if (wert) {
         q.dat[i]=hilfinteger::divrest64(wert,a.dat[i],b.dat[0]);
     }
     else q.dat[i]=a.dat[i]/b.dat[0];
     wert=a.dat[i]-q.dat[i]*b.dat[0];
 }
 if (wert) {
     r.dat=new unsigned[1];
     r.dat[0]=wert;
     r.laenge=1;
 }
}


// Division  a=q*b+wert. Prem.: q=0! and abs(b.laenge)=1
void integer::divsingle(const integer& a,const integer &b,integer& q)
{
 int la=hilfinteger::abs(a.laenge);

 if (!la) return;
 if (la==1) {
     q=integer(a.dat[0]/b.dat[0]);
     return;
 }

 int lq;
 unsigned wert=0;

 if (a.dat[la-1]<b.dat[0]) {
     lq=la-1;
     la--;
     wert=a.dat[la];
 }
 else lq=la;
 q.dat= new unsigned[lq];
 q.laenge=lq;
 for (int i=la-1;i>=0;i--) {
     if (wert) {
         q.dat[i]=hilfinteger::divrest64(wert,a.dat[i],b.dat[0]);
     }
     else q.dat[i]=a.dat[i]/b.dat[0];
     wert=a.dat[i]-q.dat[i]*b.dat[0];
 }
}


void integer::restsingle(const integer& a,const integer& b,integer & r)
{
 int la = hilfinteger::abs(a.laenge);

 if (!la) return;
 if (la==1) {
     r=integer(a.dat[0]%b.dat[0]);
     return;
 }

 unsigned wert=0,q;

 if (a.dat[la-1]<b.dat[0]) {
     la--;
     wert=a.dat[la];
 }
 for (int i=la-1;i>=0;i--) {
     if (wert) {
         q=hilfinteger::divrest64(wert,a.dat[i],b.dat[0]);
     }
     else q=a.dat[i]/b.dat[0];
     wert=a.dat[i]-q*b.dat[0];
 }
 if (wert) {
     r.dat=new unsigned[1];
     r.dat[0]=wert;
     r.laenge=1;
 }
}



void integer::divrest(const integer &a,const integer &b,integer &q,integer &r)
{

 // triviale Faelle:
 if (b.dat==NULL) {
     Error::error("\nERROR: val::divrest: Division by zero - integer!");
 }
 if (a.dat==NULL || a.abslower(b)) {
     q=0;
     r=a;
     return;
 }

 int lb=hilfinteger::abs(b.laenge),signa,signb;

 signa = (a.laenge<0)? -1:1;
 signb = (b.laenge<0)? -1:1;
 if (lb==1) {
     divrestsingle(a,b,q,r);
     q.laenge*=signa*signb;
     r.laenge*=signa;
     return;
 }

 int i,la=hilfinteger::abs(a.laenge),lq,lr,normal;
 unsigned *pr,*px,*py,q1,f=1;
 integer x,y,z,qz;



 //-------------------------------------------------------
 // Get length of q:

 lq=la-lb;
 y.dat=a.dat;
 y.dat+=lq;
 y.laenge=lb;
 if (!y.abslower(b)) lq++;
 q.dat=new unsigned[lq];
 q.laenge=lq;
 //------------------------------------------------------
 // Normalize if necessary:
 if (b.dat[lb-1]>=hilfinteger::highon) { // highon=2^(numberofbits-1)
     normal=0;
     y.dat=b.dat;
     y.laenge=lb;
     x.dat=a.dat;
     x.laenge=la;
 }
 else {
     normal=1;
     y.laenge=lb;
     y.dat=new unsigned[lb+1];
     f=unsigned((__int64)((__int64) (1)<<hilfinteger::numberofbits)/(__int64)(b.dat[lb-1]+1));
     x.dat=new unsigned[la+1];
     hilfinteger::multbyunsigned(b.dat,lb,f,y.dat);
     x.laenge=la=hilfinteger::multbyunsigned(a.dat,la,f,x.dat);
 }

 //----------------------------------------------------------



 pr=new unsigned[lr=lb+1];
 // set py for trialdivision:
 py=y.dat;py+=lb-2;
 // ---------------
 z.dat=new unsigned[lb+1];
 z.laenge=lb;
 for (i=la-lb;i<la;i++) z.dat[i-(la-lb)]=x.dat[i];

 if (z>=y) {
     q.dat[lq-1]=1;
     lq--;
     z.subto(y);
 }
 la-=lb;


 for (;lq>0;lq--,la--) {

     for (i=z.laenge;i>0;i--) z.dat[i]=z.dat[i-1];
     z.laenge++;

     z.dat[0]=x.dat[la-1];
     if (z<y) {
         q.dat[lq-1]=0;
         //std::cout<<std::endl<<0;
         continue;
     }
     if (z.laenge==lb) {
         q.dat[lq-1]=1;
         z.subto(y);
         continue;
     }
     // trial division:
     // Apply trial division:
     // put  px in right position;
     px=z.dat;
     px+=lb-2; // == z.laenge-3;
     if (px[2]==py[1]) {
         // 1. case:
         if (px[1]<py[0]) {
             q1=(py[0]-px[1])/py[1];
             if (q1==0) q1=~0;
             else q1=0-q1;
         }
         else q1=~0;
     }
     else {
         q1=hilfinteger::trialdiv(px,py);
     }
     qz.laenge=hilfinteger::multbyunsigned(y.dat,lb,q1,pr);
     qz.dat=pr;
     if (z.abslower(qz)) {
         q1--;
         qz.laenge=hilfinteger::multbyunsigned(y.dat,lb,q1,pr);
         qz.dat=pr;
     }
     q.dat[lq-1]=q1;
     z.subto(qz);
     qz.dat=NULL;
 }

 for (i=z.laenge-1;i>=0;i--)
     if (z.dat[i]) break;
 z.laenge=i+1;

 if (normal) {
     qz.laenge=0;
     divrestsingle(z,integer(f),r,qz);
 }
 else {
     y.dat=NULL;
     x.dat=NULL;
     if (z.laenge) {
         r.laenge=z.laenge;
         r.dat=new unsigned[r.laenge];
         for (i=0;i<r.laenge;i++) r.dat[i]=z.dat[i];
     }
 }
 delete[] pr;
 q.laenge*=signa*signb;
 r.laenge*=signa;

 return;
}


void divrem(const integer& a,const integer& b,integer& q,integer& r)
{
    q=r=integer();
    integer::divrest(a,b,q,r);
}


integer integer::operator /(const integer &b) const
{
 integer q;

 if (hilfinteger::abs(b.laenge)==1) {
     divsingle(*this,b,q);
     q.laenge*=hilfinteger::signum(laenge)*hilfinteger::signum(b.laenge);
     return q;
 }
 else {
     integer r;
     divrest(*this,b,q,r);
     return q;
 }
}



integer& integer::operator /=(const integer &b)
{
 integer q;

 if (hilfinteger::abs(b.laenge)==1) {
     divsingle(*this,b,q);
     q.laenge*=hilfinteger::signum(laenge)*hilfinteger::signum(b.laenge);
 }
 else {
    integer r;
    divrest(*this,b,q,r);
 }

 unsigned *ph=dat;

 dat=q.dat;
 laenge=q.laenge;
 q.dat=ph;

 return *this;
}



integer integer::operator %(const integer &b) const
{
 integer r;

 if (hilfinteger::abs(b.laenge)==1) {
     restsingle(*this,b,r);
     r.laenge*=hilfinteger::signum(laenge);
     return r;
 }
 else {
     integer q;
     divrest(*this,b,q,r);
     return r;
 }
}


integer& integer::operator %=(const integer &b)
{
 integer r;

 if (hilfinteger::abs(b.laenge)==1) {
     restsingle(*this,b,r);
     r.laenge*=hilfinteger::signum(laenge);
 }
 else  {
    integer q;
    divrest(*this,b,q,r);
 }

 unsigned *ph=dat;

 dat=r.dat;
 laenge=r.laenge;
 r.dat=ph;

 return *this;
}



void integer::truncby(int k)
{
 if (k==0) return;
 int l = hilfinteger::abs(laenge);

 if (k>=l) {
     laenge=0;
     return;
 }

 l-=k;
 for (int i=0;i<l;i++) dat[i]=dat[i+k];
 if (laenge<0) laenge=-l;
 else laenge=l;
 return;
}

// Prem.: At most numberofbits-1 bits are shifted
void integer::shiftrightby(int r)
{
 if (!r) return;
 int l= hilfinteger::abs(laenge);
 if (!l) return;

 int i,s=hilfinteger::numberofbits-r,lx=l;
 unsigned last,wert;

 last=dat[l-1]>>r;

 if (!last) l=lx-1;

 if (!l) return;

 for (i=0;i<lx-1;i++) {
     wert=dat[i+1]<<s;
     dat[i]=wert | (dat[i]>>r);
 }
 if (last) dat[lx-1]=last;

 if (laenge<0) laenge=-l;
 else laenge=l;

}


integer integer::gcdbin(const integer &a1,const integer &b1)
{
 integer a(a1),b(b1);

 int s,k,d=1,i,la=hilfinteger::abs(a.laenge),lb=hilfinteger::abs(b.laenge);
 integer h;

 if (la>lb) {
     h=a%b;
     delete[] a.dat;
     a.dat=b.dat;
     a.laenge=b.laenge;
     b.dat=h.dat;
     b.laenge=h.laenge;
     h.dat=NULL;
 }
 else if (lb>la) {
     h=b%a;
     delete[] b.dat;
     b.dat=a.dat;
     b.laenge=a.laenge;
     a.dat=h.dat;
     a.laenge=h.laenge;
     h.dat=NULL;
 }

 if (a.laenge==0) return b;
 if (b.laenge==0) return a;

 for (s=0;!a.dat[s] && !b.dat[s];s++) ;
 if (s) {
     a.truncby(s);
     b.truncby(s);
 }

 for (k=0;!(a.dat[0]&d) && !(b.dat[0]&d);k++) d=d<<1;
 if (k) {
     a.shiftrightby(k);
     b.shiftrightby(k);
 }

 // At least one of both is odd.
 while (a.laenge) {
     for (i=0;!a.dat[i];i++) ;
     if (i) a.truncby(i);
     for (i=0,d=1;!(a.dat[0]&d);i++) d=d<<1;
     if (i)  a.shiftrightby(i);
     for (i=0;!b.dat[i];i++) ;
     if (i)  b.truncby(i);
     for (i=0,d=1;!(b.dat[0]&d);i++) d=d<<1;
     if (i) b.shiftrightby(i);
     // now both are odd
     if (a.abslower(b)) {
         b.subto(a);
         b.shiftrightby(1);
     }
     else {
         a.subto(b);
         a.shiftrightby(1);
     }
 }
 integer c(shift(b,s));
 c=c<<k;
 return c;
}

//Prem.: a divides b
integer EDIV(const integer& b,const integer& a)
{
 if (hilfinteger::abs(a.laenge)==1) {
     integer q;
     if (a.dat[0]==1) {
        q=b;
        if (a.laenge==-1) q.changesign();
        return q;
     }
     integer::divsingle(b,a,q);
     q.laenge*=hilfinteger::signum(a.laenge)*hilfinteger::signum(b.laenge);
     return q;
 }
 else return hilfinteger::exactdivision(b,a);
}


// Prem.: a divides *this
integer& integer::EDIVBY(const integer& a)
{
 unsigned *ph;

 if (hilfinteger::abs(a.laenge)==1) {
     integer q;
     if (a.dat[0]==1) {
        if (a.laenge==1) return *this;
        changesign();
        return *this;
     }

     divsingle(*this,a,q);
     q.laenge*=hilfinteger::signum(a.laenge)*hilfinteger::signum(laenge);
     ph=dat;
     dat=q.dat;
     laenge=q.laenge;
     q.dat=ph;
     return *this;
 }
 else {
     integer q;
     q=hilfinteger::exactdivision(*this,a);
     ph=dat;
     dat=q.dat;
     laenge=q.laenge;
     q.dat=ph;
     return *this;
 }
}


integer gcdeuk(integer a,integer b)
{
 int i;
 integer h;

 if (a.dat==NULL) return b;
 if (b.dat==NULL) return a;
 h.dat=new unsigned[h.laenge=hilfinteger::abs(b.laenge)];

 while (b.dat!=NULL) {
     h.laenge=hilfinteger::abs(b.laenge);
     for (i=0;i<h.laenge;i++) h.dat[i]=b.dat[i];
     b=a%b;
     a=h;
 }
 return a;
}


integer ggTspez(const integer &a,const integer &b)
{
    if (a.dat==NULL) return b;
    if (b.dat==NULL) return a;
    if (hilfinteger::abs(b.laenge)>1 && hilfinteger::abs(a.laenge)>1) return integer::gcdbin(a,b);
    if (hilfinteger::abs(b.laenge)==1 && hilfinteger::abs(a.laenge)==1) {
        if (a.dat[0]==1 || b.dat[0]==1) return integer(1);
        return integer(hilfinteger::ggT(a.dat[0],b.dat[0]));
    }

    if (hilfinteger::abs(b.laenge)>1) {
        integer r(hilfinteger::restsingle(b,a));
        if (r.dat==NULL) return a;
        return integer(hilfinteger::ggT(a.dat[0],r.dat[0]));
    }
    else {
        integer r(hilfinteger::restsingle(a,b));
        if (r.dat==NULL) return b;
        return integer(hilfinteger::ggT(b.dat[0],r.dat[0]));
    }
}


integer gcd(const integer &a,const integer &b)
{
    return ggTspez(a,b);
}


integer lcm(const integer& a,const integer& b)
{
    if (a.isNull() || b.isNull()) return integer(0);
    else if (a.abslength()>b.abslength()) return EDIV(a,ggTspez(a,b))*b;
    else return EDIV(b,ggTspez(a,b))*a;
}


integer abs(const integer& a)
{
    integer x(a);
    if (x.laenge<0) x.laenge*=-1;
    return x;
}


// non-recursive function  eudlid. Returns gcd(a0,b0)!
// x0*a0 + y0*b0 = gcd(a0,b0)
integer euclid(integer a0,integer a1,integer& x0,integer& y0)
{
 integer h,q,x1=integer(0),y1=integer(1);

 x0=integer(1);y0= integer(0);

 while (!a1.iszero()) {
       q=a0/a1;
       h=std::move(a1);
       a1=a0%h;
       a0=std::move(h);
       h= std::move(x1);
       x1=x0 - q*h;
       x0=std::move(h);
       h= std::move(y1);
       y1= y0 - q*h;
       y0= std::move(h);
 }

 if (a0.signum()<0) {
    a0.changesign();
    x0.changesign();
    y0.changesign();
 }

 return a0;
}


integer integer::char_to_integer(char* buf,int l)
{
  int i,sign=1,exp=0,kommastellen=0,signexp=1;
  integer ten(10),x;

  if (buf==NULL) return x;

  if (l==-1) for (l=0;buf[l]!='\0';l++) ;

  if (l==0) return x;


  // delete all 0 and - at the beginning:
  for (i=0;i<l;) {
      if (buf[i]=='-') {
          sign*=-1;
          i++;
      }
      else if (buf[i]=='0') {
          i++;
      }
      else break;
  }
  if (i==l) {
     return x;
  }
  // get digits:
  for (;i<l;i++) {
      if (buf[i]=='.' || buf[i]==',') {
        if (kommastellen) continue;
        kommastellen=1;
        continue;
      }
      if (buf[i]=='e') {i++;break;}
      if (kommastellen) kommastellen++;
      x*=ten;
      x+=integer(int(buf[i]-48));
  }

  if (kommastellen) kommastellen--;

  // Get exponent:
  for (;i<l;i++) {
      if (buf[i]=='+' || buf[i]=='0') continue;
      else if (buf[i]=='-') {signexp*=-1;continue;}
      else break;
  }
  for (;i<l;i++) {
      exp*=10;
      exp+= int(buf[i]-48);
  }

  if (signexp==-1) exp*=-1;

  exp-=kommastellen;

  if (exp>=0) {
    for (i=0;i<exp;i++) x*=ten;
  }
  else {
    exp*=-1;
    for (i=0;i<exp && x.dat!=NULL ;i++) x/=ten;
  }

  if (sign==-1) x.laenge*=-1;
  return x;
}


std::istream& operator >>(std::istream& is,integer& x)
{
  char *buf,*hbuf;
  int i,l,k=1,wert=0;
  integer ten(10);

  buf = new char[1000];

  l=0;
  while ((wert!=10) && (wert!=32) && is) {
    wert=is.get();
    if ((wert==10) || (wert==32)){
        if (l!=0) buf[l]='\0';
        else wert=0;
    }
    else buf[l++]=wert;
    if (l>=(k*1000)) {
        k++;
        hbuf= new char[k*1000];
        for (i=0;i<((k-1)*1000);i++) hbuf[i]=buf[i];
        delete[] buf;
        buf=hbuf;
        hbuf=NULL;
    }

    if (!is) l--;
  }
  x=integer::char_to_integer(buf,l);

  delete[] buf;
  return is;
}


std::ostream& operator <<(std::ostream& os,integer x)
{
 integer q,r,faktor(10),ten(10),eins(1);

 if (x.dat==NULL) {
     os<<'0';
     return os;
 }

 if (x.laenge<0) os<<'-';

 if (integer::Output_Style==integer::INTEGER) {
        if (x.laenge==1 || x.laenge==-1) {
            os<<x.dat[0];
            return os;
        }
        while (faktor.abslower(x)) faktor*=ten;

        while (x.dat!=NULL) {
            faktor/=ten;
            integer::divrest(x,faktor,q,r);
            if (q.dat==NULL) os<<'0';
            else os<<q.dat[0];
            x=std::move(r);
            q=0;
        }
        while (eins.abslower(faktor)) {
            os<<'0';
            faktor/=ten;
        }
 }
 else {
     int i=1,exp=0,numberofzeros=0,iskomma=0;

     while (faktor.abslower(x)) {faktor*=ten;exp++;}

     while (x.dat!=NULL && i<=integer::Output_Digits) {
            faktor/=ten;exp--;
            integer::divrest(x,faktor,q,r);
            if (q.dat==NULL) numberofzeros++;
            else if (q.dat[0]==10) {
                os<<'1';i++;
            }
            else {
                 if (!iskomma && i>1) {os<<'.';iskomma=1;}
                 for(;numberofzeros>0;numberofzeros--) os<<'0';
                 os<<q.dat[0];
            }
            i++;
            x=std::move(r);
            q=0;
        }
        while (eins.abslower(faktor) && i<=integer::Output_Digits) {
              faktor/=ten;
              i++;
              exp--;
        }
        if ((exp +i-1)>0) os<<"e+0"<<exp+(i-1);
 }

 return os;
}

} //namespace val

#endif  //INTEGER_CPP
