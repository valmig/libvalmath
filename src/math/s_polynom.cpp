
#include <s_polynom.h>
#include <Glist.h>
#include <modq.h>
#include <rational.h>


namespace val
{


//


template <>
std::atomic<int> s_polynom<modq>::mnumber(0);

template <>
std::atomic<int> s_polynom<integer>::mnumber(0);

template <>
std::atomic<int> s_polynom<modq>::nreduction(0);

template <>
std::atomic<int> s_polynom<integer>::nreduction(0);

template <>
const val::modq s_polynom<val::modq>::zero(0);

template <>
s_expo s_polynom<modq>::neutral(0);

//


template <>
template <> const s_polynom<modq>& s_polynom<modq>::convertfrom(const s_polynom<integer> &f)
{
    if (f.iszero()) {
        del();
        return *this;
    }

    integer d(modq::q);
    modq wert;

    term *p;
    s_polynomIterator<integer> q=f;

    // Set head:

    for (;q;q++) {
        if ((wert=modq(int(q.actualcoef()%d)))!=zero) break;
    }
    if (!q) {
        del();
        return *this;
    }
    else if (head!=NULL) {
        head->coeff=wert;
        head->X=q.actualterm();
    }
    else {
        head=new term(wert,q.actualterm());
    }
    p=head;q++;

    for (;q;q++) {
        if ((wert=modq(int(q.actualcoef()%d)))!=zero) {
            if (p->next!=NULL) {
                p->next->coeff=wert;
                p->next->X=q.actualterm();
            }
            else p->next=new term(wert,q.actualterm());
            p=p->next;
        }
    }
    if (p->next!=NULL) {
        s_polynom<modq> g;
        g.head = p->next;
        p->next=NULL;
    }
    return *this;
}


//

template <>
void s_polynom<val::modq>::minusmalmonom(const s_polynom<val::modq> &g,const val::modq &a,const s_expo &Y)
{
 int c;
 val::modq b,a1=-a;
 s_expo Z;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 Z=q->X;
	 Z*=Y;
	 c=expocompare(Z,head->X);
	 if (c==-1) break;
	 else if (c==1) {		// head->X < Z 
		 b=q->coeff;
		 b*=a1;
		 r=new term(b,Z,head);
		 head=r;
		 q=q->next;
		 break;
	 }
	 else {
         b=q->coeff;
         b*=a1;
		 head->coeff+=b; 
		 if (head->coeff==zero) {
			 r=head;
			 head=head->next;
			 delete r;
		 }
		 else { q=q->next; break; }
	 }
 }

 p=head;
 if (p!=NULL) {
	 // Invariant: p->X > q->X :
	 while (q!=NULL && p->next!=NULL) {
		 Z=q->X;
		 Z*=Y;
		 c=expocompare(Z,p->next->X);
		 if (c==-1) p=p->next;
		 else if (c==1) {		// p->next->X < Z 
			 b=q->coeff;
			 b*=a1;
			 r=new term(b,Z,p->next);
			 p->next=r;
			 p=p->next;
			 q=q->next;
		 }
		 else {
             b=q->coeff;
             b*=a1;
			 p->next->coeff+=b;
			 if (p->next->coeff==zero) {
				 r=p->next;
				 p->next=r->next;
				 delete r;
			 }
			 else p=p->next;
			 q=q->next;
		 }
	 }
 }
 if (head==NULL) {
	 for(;q!=NULL;q=q->next) {
		 b=q->coeff;
		 b*=a1;
		 if (b!=zero) {
			 head = new term(b,q->X*Y);
			 p=head;
			 q=q->next;
			 break;
		 }
	 }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff;
     b*=a1;
	 if (b!=zero) {
		 r = new term(b,q->X*Y);
		 p->next=r;
		 p=p->next;
	 }
 }

 return;
}


template <>
void s_polynom<modq>::minusmalmonom(const s_polynom<modq> &g,const modq &a,const s_expo &Y,const integer& d)
{
 int c;
 integer d1;
 modq b,a1=-a;
 s_expo Z;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 Z=q->X*Y;
	 d1=q->deg +d;
	 c=wexpocompare(Z,d1,head->X,head->deg);
	 if (c==-1) break;
	 else if (c==1) {     // head->X < Z
		 b=q->coeff;
		 b*=a1;
		 r=new term(b,Z,std::move(d1),head);
		 this->head=r;
		 q=q->next;
		 break;
	 }
	 else {

         b=q->coeff;
         b*=a1;
		 head->coeff+=b; 
		 if (head->coeff==zero) {
			 r=this->head;
			 head=head->next;
			 delete r;
		 }
		 else { q=q->next; break; }
	 }
 }

 p=head;
 if (p!=NULL) {
	 // Invariant: p->X > q->X 
	 while (q!=NULL && p->next!=NULL) {
		 Z=q->X*Y;
		 d1=q->deg+d;
		 c=wexpocompare(Z,d1,p->next->X,p->next->deg);
		 if (c==-1) p=p->next;
		 else if (c==1) {		// p->next->X < Z
			 b=q->coeff;
			 b*=a1;
			 r=new term(b,Z,std::move(d1),p->next);
			 p->next=r;
			 p=p->next;
			 q=q->next;
		 }
		 else {
             b=q->coeff*a1;
			 p->next->coeff+=b;
			 if (p->next->coeff==zero) {
				 r=p->next;
				 p->next=r->next;
				 delete r;
			 }
			 else p=p->next;
			 q=q->next;
		 }
	 }
 }
 if (this->head==NULL) {
	 for(;q!=NULL;q=q->next) {
		 b=q->coeff*a1;
		 if (b!=zero) {
			 head = new term(b,q->X*Y,q->deg+d);
			 p=head;
			 q=q->next;
			 break;
		 }
	 }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff*a1;
	 if (b!=zero) {
		 r = new term(b,q->X*Y,q->deg+d);
		 p->next=r;
		 p=p->next;
	 }
 }

 return;
}


//

template <>
void s_polynom<integer>::changesign(integer &a)
{
    a.changesign();
}


template<>
void s_polynom<rational>::changesign(rational &a)
{
    a.changesign();
}
//

template <>
void s_polynom<integer>::edivby(const integer& a)
{
    if (a==integer(-1)) {
        for (term *p=head;p!=NULL;p=p->next) p->coeff.changesign();
        return;
    }
    for (term *p=head;p!=NULL;p=p->next) p->coeff.EDIVBY(a);
}


template<>
const s_polynom<int>& s_polynom<int>::normalize()
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
const s_polynom<integer>& s_polynom<integer>::normalize()
{
    if (head==NULL) return *this;
    if (head->coeff.isunit()) return *this;

    integer ggT(head->coeff);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::ggTspez(ggT,p->coeff);
        if (ggT.isunit()) return *this;
    }
    for (p=head;p!=NULL;p=p->next) p->coeff.EDIVBY(ggT);

    return *this;
}


template <>
int s_polynom<int>::content() const
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
integer s_polynom<integer>::content() const
{
    integer eins(1),minuseins(-1);
    if (head==NULL) return eins;
    if (head->coeff.isunit()) return eins;

    integer ggT(head->coeff);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::ggTspez(ggT,p->coeff);
        if (ggT.isunit()) return eins;
    }

    return ggT;
}

template <>
rational s_polynom<rational>::content() const
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

//

template <>
int s_polynom<integer>::reduction(const Glist<s_polynom<integer> >& P,int top,int interred,int degoption)
{
 int reduzbar,reduced=0,interreduced=0;
 unsigned maxnG=0-1,i;
 integer lcm,a1,gcd,div,hgcd,b,eins(1),minuseins(-1);
 s_expo X;
 s_polynom<integer> g,h;
 term *p,*r;
 GlistIterator<s_polynom<integer> > ItP;
 integer d;

 if (head==NULL) return 0;
 if (P.isempty()) {normalize();return 0;}

 if (interred) top=0;
 if (top) interred=0;

 // 1. top-reductions
 if (!interred) {
    do {

        reduzbar=0;
        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)|head->X) {
                reduzbar=reduced=1;
                break;
            }
        }
        if (reduzbar) {
            nreduction++;
            g.head=ItP.getelement().head;
            hgcd=ggTspez(head->coeff,g.head->coeff);
            X=head->X/g.head->X;
            if (degoption) d=head->deg-g.head->deg;
            if (hgcd.isunit()) {
                b=std::move(head->coeff);
                p=head;
                head = head->next;
                delete p;
                if (g.head->coeff!=eins) operator *=(g.head->coeff);
            }
            else {
                a1=EDIV(g.head->coeff,hgcd);
                b=EDIV(head->coeff,hgcd);
                p=head;
                head=head->next;
                delete p;
                if (a1!=eins) operator*=(a1);
            }
            g.head=g.head->next;
            if (degoption) minusmalmonom(g,b,X,d);
            else minusmalmonom(g,b,X);
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

 // 2. inter-reductions
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
	 }
	 if (reduzbar) {
         nreduction++;
		 h.head=p->next;p->next=NULL;
         g.head=ItP.getelement().head;
		 hgcd=ggTspez(h.head->coeff,g.head->coeff);
         X=h.head->X/g.head->X;
         if (degoption) d=h.head->deg - g.head->deg;
         a1=EDIV(g.head->coeff,hgcd);
         b=EDIV(h.head->coeff,hgcd);
         r=h.head;
         h.head=h.head->next;
         delete r;
         p->next=h.head;
         if (a1 != eins) {operator*=(a1); gcd*=a1;}
		 g.head=g.head->next;
		 if (degoption) h.minusmalmonom(g,b,X,d);
		 else h.minusmalmonom(g,b,X);
		 div=gcd;
		 for (r=h.head;r!=NULL;r=r->next) {
			 if (div.isunit()) break;
			 else (div=ggTspez(div,r->coeff));
		 }
		 p->next=h.head;
		 if (!div.isunit()) {gcd.EDIVBY(div);edivby(div);}
	 }
	 else {
		 p=p->next;
		 if (p!=NULL && interreduced && !gcd.isunit()) gcd=ggTspez(gcd,p->coeff);
	 }
 }
 //normalize();
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


template <>
int s_polynom<modq>::reduction(const Glist<s_polynom<modq> >& P,int top,int interred,int degoption)
{
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modq b;
 s_expo X;
 s_polynom<modq> g,h;
 term *p,*r;
 GlistIterator<s_polynom<modq> > ItP;
 integer d;


 if (head==NULL) return 0;
 if (P.isempty()) {
	 normalize();
	 return 0;
 }

 if (interred) top=0;
 if (top) interred=0;

 // 1. top-reductions
 if (!interred) {
    do {
        reduzbar=0;

        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)|head->X) {
			 reduzbar=reduced=1;
			 break;
            }
        }
        if (reduzbar) {
		    nreduction++;
            g.head=ItP.getelement().head;
            b=head->coeff/g.head->coeff;
            X=head->X/ItP.getelement().head->X;  // Y/X
            if (degoption) d=head->deg - ItP.getelement().head->deg;
            g.head=g.head->next;
            p=head;
            head=head->next;
            delete p;
            if (degoption) minusmalmonom(g,b,X,d);
            else minusmalmonom(g,b,X);
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
	 }
	 if (reduzbar) {
		 nreduction++;
		 h.head=p->next;
		 g.head=ItP.getelement().head;
		 b=h.head->coeff/g.head->coeff;
		 X=p->next->X/g.head->X;  // Y/X
		 if (degoption) d=p->next->deg-g.head->deg;
		 p->next=NULL;
		 g.head=g.head->next;
		 r=h.head;
		 h.head=h.head->next;
		 delete r;
		 if (degoption) h.minusmalmonom(g,b,X,d);
		 else h.minusmalmonom(g,b,X);
		 p->next=h.head;
	 }
	 else p=p->next;
 }
 normalize();
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


template<>
void s_polynom<modq>::liftpolynom(const Glist<s_polynom<modq> > &inwG,const Glist<s_polynom<modq> > &G,val::vector<char> &remain)
{
 int reduzbar,i,ismonom=0,red=0;
 integer d;
 s_expo X;
 modq b;
 s_polynom<modq> g,F;
 term *p;
 GlistIterator<s_polynom<modq> > pinwG,pG;

 if (head==NULL) return;
 if (inwG.isempty()) return;

 if (head->next==NULL) ismonom=1;

 // top-reductions are sufficient
 do {
	 reduzbar=0;

	 for (pinwG=inwG,pG=G,i=0;pinwG;pinwG++,pG++,i++) {
		 if (pinwG().head->X|head->X) {
			 reduzbar=1;
			 break;
		 }
	 }
	 if (reduzbar) {
		 if (!red && ismonom && (pinwG().head->next==NULL) && (pinwG().head->X==head->X)) {
			 del();
			 remain[i]=1;
			 return;
		 }
		 red=1;
		 nreduction++;
		 g.head=pinwG().head;
		 b=head->coeff/g.head->coeff;
		 X=head->X/pinwG().head->X;  
		 d=head->deg - pinwG().head->deg;
		 g.head=g.head->next;
		 p=head;
		 head=head->next;
		 delete p;
		 minusmalmonom(g,b,X,d);
		 F.minusmalmonom(pG(),-b,X,d);
		 if (head==NULL) break;
	 }
 }
 while (reduzbar);

 if (head!=NULL) {
	 Error::error("\nBy lifting polynomial, it does not reduced to 0!");
 }
 head=F.head;
 F.head=NULL;
 g.head=NULL;
}



template<>
void s_polynom<integer>::liftpolynom(const Glist<s_polynom<integer> > &inwG,const Glist<s_polynom<integer> > &G,val::vector<char> &remain)
{
 int reduzbar,i,ismonom=0,red=0,notfirst=0;
 integer d;
 integer a1,hgcd,b,eins(1);
 s_expo X;
 s_polynom<integer> g,F;
 term *p;
 GlistIterator<s_polynom<integer> > ptoinwH,ptoH;

 if (head==NULL) return;
 if (inwG.isempty()) return;

 if (head->next==NULL) ismonom=1;
 // top-reductions are sufficient
 do {
	 reduzbar=0;


	 for (ptoinwH=inwG,ptoH=G,i=0;ptoinwH;ptoinwH++,ptoH++,i++) {
		 if ( ptoinwH().head->X | head->X) {
			 reduzbar=1;
			 break;
		 }
	 }
	 if (reduzbar) {
		 if (!red && ismonom && (ptoinwH().head->next==NULL) && (ptoinwH().head->X== head->X)) {
			 del();
			 remain[i]=1;
			 return;
		 }
		 red=1;
		 nreduction++;
		 g.head=ptoinwH().head;
		 hgcd=ggTspez(head->coeff,g.head->coeff);
		 X=head->X/g.head->X;
		 d=head->deg - g.head->deg;

        if (hgcd.isunit()) {
                b=std::move(head->coeff);
                p=head;
                head = head->next;
                delete p;
                if (g.head->coeff!=eins) {
                    operator *=(g.head->coeff);
                    if (notfirst) F*=(g.head->coeff);
                }
        }
        else {
                a1=EDIV(g.head->coeff,hgcd);
                b=EDIV(head->coeff,hgcd);
                p=head;
                head=head->next;
                delete p;
                if (a1!=eins) {
                   operator*=(a1);
                   if (notfirst) F*=a1;
                }
        }

		 g.head=g.head->next;
		 minusmalmonom(g,b,X,d);
		 notfirst=1;
		 b.changesign();
		 F.minusmalmonom(ptoH(),b,X,d);
		 if (head==NULL) break;
	 }
 }
 while (reduzbar);

 if (head!=NULL) {
	 Error::error("\nBy lifting polynomial, it does not reduced to 0!");
 }
 F.normalize();
 head=F.head;
 F.head=NULL;
 g.head=NULL;
}


} // end namespace val
