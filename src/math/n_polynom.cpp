

#include <n_polynom.h>
#include <val_utils.h>
#include <rational.h>
#include <modq.h>
#include <Glist.h>


namespace val
{

//  -----------------------------------------------------------------------------------
n_expo::n_expo(int value,int d)
{
    if (d<=0) {
        coeff=nullptr;
        dim=0;
        return;
    }
    coeff= new int[dim=d];
    for (int i =0;i<d;++i) coeff[i] = value;
}

n_expo::n_expo(const n_expo &X)
{
    dim = X.dim;
    if (X.coeff==nullptr) {
        coeff=nullptr;
        return;
    }
    coeff = new int[dim];
    for (int i=0;i<dim;++i) coeff[i] = X.coeff[i];
}

n_expo::n_expo(n_expo &&X)
{
    coeff = X.coeff;
    dim = X.dim;
    X.coeff = nullptr;
    X.dim = 0;
}

//
const n_expo& n_expo::operator=(const n_expo& X)
{
    if (coeff==X.coeff) return *this;
    if (dim != X.dim) {
        delete[] coeff;
        coeff = new int[dim=X.dim];
    }
    for (int i=0;i<dim;++i) coeff[i] = X.coeff[i];
    return *this;
}


const n_expo& n_expo::operator=(n_expo&& X)
{
    if (coeff==X.coeff) return *this;
    if (coeff!=nullptr) delete[] coeff;
    coeff = X.coeff;
    dim=X.dim;
    X.coeff=nullptr;
    X.dim= 0;
    return *this;
}

int n_expo::iszero() const
{
    if (coeff==nullptr) return 1;
    for (int i=0;i<dim;++i)
        if (coeff[i]) return 0;
    return 1;
}


int n_expo::operator==(const n_expo& X) const
{
    int i,n=Max(dim,X.dimension());
    for (i=0;i<n;i++) if (operator[](i)!=X(i)) return 0;
    return 1;
}

int n_expo::operator[](int i) const
{
    if (i<0 || i>=dim) return 0;
    else return coeff[i];
}

int& n_expo::operator[](int i)
{
    if (i<0 || i>=dim) Error::error("\nERROR: n_expo::operator[] : NULL-Pointer!");
    return coeff[i];
}


std::ostream& operator <<(std::ostream& os,const n_expo& X)
{
    if (X.isempty()) {
        os<<'0'<<"  ";
        return os;
    }
    for (int i=0; i<X.dimension();++i) os<<X[i]<<"  ";
    return os;
}

//    -----------------------------------------------------------------------------------------------------------------------------------------


//

int n_expo::ordtype=-2;

matrix<int> n_expo::ordmatrix(0);

//

int n_expo::operator<(const n_expo& y) const
{
    const n_expo& x= *this;
    int i,j,m1,m2,n=Max(x.dimension(),y.dimension());

    m1=1;m2=0;
    if (ordtype==0)  { // totaldegree lexik.

        m1=m2=0;
        for (i=0;i<n;i++) {
            m1+=x[i];
            m2+=y[i];
        }

        if (m1<m2) return 1;
        if (m1>m2) return 0;
    }

    if (ordtype==-1 || m1==m2) { //lexik.
        for (i=0;i<n;i++) {
            if (x[i] != y[i]) return (x[i] < y[i]);
        }
        return 0;
    }

    if (ordtype==-2) {
        m1=m2=0;
        for (i=0;i<n;i++) {
            m1+=x[i];
            m2+=y[i];
        }
        if (m1<m2) return 1;
        if (m1>m2) return 0;

        for (i=n-1;i>0;i--) {
            if (x[i]>y[i]) return 1;
            else if (x[i]<y[i]) return 0;
        }
        return 0;
    }

    if (ordtype!=-1) {
        const matrix<int> &M=ordmatrix;
        int m=M.numberofcolumns();

        for (i=0;i<m;i++) {
            m1=m2=0;
            for (j=0;j<m;j++) {
                m1+=x[j]*M(i,j);
                m2+=y[j]*M(i,j);
            }
            if (m1==m2) {
                int ord=ordtype,k;
                ordtype=-1;
                k=operator<(y);
                ordtype=ord;
                return k;
            }
            else return (m1<m2);
        }
        return 0;
    }
    return 0;
}

//


n_expo operator*(const n_expo& x,const n_expo& y)
{
    int i,n=val::Max(x.dimension(),y.dimension());

    n_expo z(n);

    for (i=0;i<n;i++) z(i)=x(i) + y(i);
    return z;
}


const n_expo& operator*=(n_expo& x,const n_expo& y)
{
    int i,m,n=val::Max(m=x.dimension(),y.dimension());

    if (n>m) {
        n_expo z(y);
        for (i=0;i<m;i++) z(i)+=x(i);
        x=std::move(z);
    }
    else for (i=0;i<n;i++) x(i)+=y(i);
    return x;
}


n_expo operator/(const n_expo& x,const n_expo& y)
{
    int i,n=val::Max(x.dimension(),y.dimension());

    n_expo z(n);

    for (i=0;i<n;i++) z(i)=x(i) - y(i);
    return z;
}


const n_expo& operator/=(n_expo& x,const n_expo& y)
{
    int i,m,n=val::Max(m=x.dimension(),y.dimension());

    if (n>m) {
        n_expo z(n);
        for (i=0;i<m;i++) z(i)=x(i)-y(i);
        for (;i<n;i++) z(i) = -y(i);
        x=std::move(z);
    }
    else for (i=0;i<n;i++) x(i)-=y(i);
    return x;
}

int operator|(const n_expo& x,const n_expo& y)
{
    int i,n=Max(x.dimension(),y.dimension());
    for (i=0;i<n;i++) if (y(i)<x(i)) return 0;
    return 1;
}


n_expo lcm(const n_expo& x,const n_expo &y)
{
    int i,n=Max(x.dimension(),y.dimension());
    n_expo z(n);

    for (i=0;i<n;i++) z(i)=Max(x(i),y(i));
    return z;
}


n_expo char_to_nexpo(const char* s,int l)
{
    n_expo LES;

    if (l==0) return LES;
    int first=0,i,k,n;
    std::string zeichen="";

    for (i=0;i<l;i++)
        if (s[i]!='\n' && s[i]!=' ' ) break;
    first=i;
    n=0;
    //Get number of columns n
    for (;i<l;) {
        if (i<l && s[i]=='\n') {break;}
        while (i<l && (s[i]!=' ' && s[i] != '\n' )) i++;
        n++;
        while (i<l && s[i]==' ') i++;
    }

    if (n==0) return LES;

    LES=val::n_expo(n);

    i=0;
    for (k=first;k<l;) {
        if (s[k]!=' ' && s[k]!='\n') {
            while (k<l && s[k]!=' ' && s[k]!='\n') {
                zeichen+=s[k];
                k++;
            }
            LES(i)=val::FromString<int>(zeichen);
            zeichen="";
            i++;
        }
        while (k<l && s[k]==' ') {
               k++;
        }
    }
    return LES;
}


std::istream& operator >>(std::istream& is,n_expo& x)
{
  char *buf,*hbuf;
  int i,l=0,wert=0,k=1;

  buf =new char[1000];


  while ((wert!=10) && is) {
    wert=is.get();
    if (wert==10){
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

  x=char_to_nexpo(buf,l);
  delete[] buf;
  return is;
}

//

template <>
void n_polynom<integer>::changesign(integer &a)
{
    a.changesign();
}


template<>
void n_polynom<rational>::changesign(rational &a)
{
    a.changesign();
}

//
template <>
void n_polynom<integer>::edivby(const integer& a)
{
    for (term *p=head;p!=NULL;p=p->next) p->coeff.EDIVBY(a);
}


template<>
const n_polynom<int>& n_polynom<int>::normalize()
{
    if (head==NULL) return *this;
    int ggT=head->coeff;
    if (ggT==1 || ggT==-1) return *this;
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::gcd(ggT,p->coeff);
        if (ggT==1 || ggT==-1) return *this;
    }
    for (p=head;p!=NULL;p=p->next) p->coeff/=ggT;

    return *this;
}


template<>
const n_polynom<integer>& n_polynom<integer>::normalize()
{
    if (head==NULL) return *this;
    integer eins(1),minuseins(-1);
    if (head->coeff==eins || head->coeff==minuseins) return *this;

    integer ggT(head->coeff);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::ggTspez(ggT,p->coeff);
        if (ggT==eins || ggT==minuseins) return *this;
    }
    for (p=head;p!=NULL;p=p->next) p->coeff.EDIVBY(ggT);

    return *this;
}


template <>
int n_polynom<int>::content() const
{
    if (head==NULL) return 1;
    if (head->coeff==1) return 1;
    if (head->coeff==-1) return -1;

    int ggT(head->coeff);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::gcd(ggT,p->coeff);
        if (ggT==1 || ggT==-1) return ggT;
    }

    return ggT;
}


template <>
integer n_polynom<integer>::content() const
{
    integer eins(1),minuseins(-1);
    if (head==NULL) return eins;
    if (head->coeff==eins) return eins;
    if (head->coeff==minuseins) return minuseins;

    integer ggT(head->coeff);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::ggTspez(ggT,p->coeff);
        if (ggT==eins || ggT==minuseins) return ggT;
    }

    return ggT;
}

template <>
rational n_polynom<rational>::content() const
{
    term *pf;
    integer z(1),n(1),eins(1),minuseins(-1);

    if (head==NULL) return rational(z,n);
    pf=head;
    z=nominator(pf->coeff);
    n = denominator(pf->coeff);
    pf=pf->next;
    for (; pf!=NULL;pf=pf->next)
    {
        if (z!=eins && z!= minuseins) z=ggTspez(z,nominator(pf->coeff));
        n=lcm(n,denominator(pf->coeff));
    }
    return rational(z,n);
}


template <>
int n_polynom<integer>::reduction(const Glist<n_polynom<integer> >& P,int top,int interred)
{
 int reduzbar,reduced=0,interreduced=0;
 unsigned maxnG=0-1,i;
 integer lcm,a1,gcd,div,hgcd,b;
 n_expo X;
 n_polynom<integer> g,h;
 term *p,*r;
 GlistIterator<n_polynom<integer> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {normalize();return 0;}

 if (interred) top=0;
 if (top) interred=0;

 // 1.: top-reductions
 if (!interred) {
    do {

        reduzbar=0;
        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)|head->X) {
                reduzbar=reduced=1;
                break;
            }
            if (head->X<ItP.getelement().head->X){  
                maxnG=i;
                break;
            }
        }
        if (reduzbar) {
            nreduction++;
            g.head=ItP.getelement().head;
            hgcd=ggTspez(head->coeff,g.head->coeff);
            a1=EDIV(g.head->coeff,hgcd);
            X=head->X/g.head->X;
            b=EDIV(head->coeff,hgcd);
            g.head=g.head->next;
            p=head;
            head=head->next;
            delete p;
            operator*=(a1);
            minusmalmonom(g,b,X);
            normalize();
            if (head==NULL) break;
        }
    }
    while (reduzbar);
 }

 if (top || head==NULL) {
     g.head=NULL;
     return reduced;
 }

 // 2. inter-reductions:
 p=head;

 while (p->next!=NULL) {
     reduzbar=0;
     for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
         if (ItP.getelement().head->X|p->next->X) {
             reduzbar=reduced=1;
             if (!interreduced) {
                 r=p->next;p->next=NULL;
                 gcd=content();
                 p->next=r;r=NULL;
                 interreduced=1;
             }
             break;
         }
         if (p->next->X<ItP.getelement().head->X){  
             maxnG=i;
             break;
         }

     }
     if (reduzbar) {
         nreduction++;
         h.head=p->next;p->next=NULL;
         g.head=ItP.getelement().head;
         hgcd=ggTspez(h.head->coeff,g.head->coeff);
         a1=EDIV(g.head->coeff,hgcd);
         X=h.head->X/g.head->X;
         b=EDIV(h.head->coeff,hgcd);
         g.head=g.head->next;
         r=h.head;
         h.head=h.head->next;
         delete r;
         p->next=h.head;
         operator*=(a1);
         gcd*=a1;
         h.minusmalmonom(g,b,X);
         div=gcd;
         for (r=h.head;r!=NULL;r=r->next) {
             if (div==integer(1) || div==integer(-1)) break;
             else (div=ggTspez(div,r->coeff));
         }

         gcd.EDIVBY(div);
         p->next=h.head;
         edivby(div);
     }
     else {
         p=p->next;
         if (p!=NULL && interreduced) gcd=ggTspez(gcd,p->coeff);
     }
 }
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


template <>
int n_polynom<modq>::reduction(const Glist<n_polynom<modq> >& P,int top,int interred)
{
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modq b;
 n_expo X;
 n_polynom<modq> g,h;
 term *p,*r;
 GlistIterator<n_polynom<modq> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {
     normalize();
     return 0;
 }

 if (interred) top=0;
 if (top) interred=0;

 // 1. top-reductions:
 if (!interred) {
    do {
        reduzbar=0;

        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)|head->X) {//(divisible(Y,X,m.X)) {//(X|Y) {
             reduzbar=reduced=1;
             break;
            }

            if (head->X<ItP.getelement().head->X){  
                maxnG=i;
                break;
            }
        }
        if (reduzbar) {
            nreduction++;
            g.head=ItP.getelement().head;
            b=head->coeff/g.head->coeff;
            X=head->X/ItP.getelement().head->X;  
            g.head=g.head->next;
            p=head;
            head=head->next;
            delete p;
            minusmalmonom(g,b,X);
            if (head==NULL) break;
        }
    }
    while (reduzbar);
 }
 if (top || head==NULL) {
     normalize();
     g.head=NULL;
     return reduced;
 }

 // 2. inter-reductions
 p=head;

 while (p->next!=NULL) {
     reduzbar=0;

     for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
         if (ItP.getelement().head->X|p->next->X){
             reduzbar=reduced=1;
             break;
         }
         if (p->next->X < ItP.getelement().head->X){  
             maxnG=i;
             break;
         }
     }
     if (reduzbar) {
         nreduction++;
         h.head=p->next;
         g.head=ItP.getelement().head;
         b=h.head->coeff/g.head->coeff;
         X=p->next->X/g.head->X;
         p->next=NULL;
         g.head=g.head->next;
         r=h.head;
         h.head=h.head->next;
         delete r;
         h.minusmalmonom(g,b,X);
         p->next=h.head;
     }
     else p=p->next;
 }
 normalize();
 h.head=NULL;
 g.head=NULL;
 return reduced;
}

}// end namespace val
