
#include <s_groebner/w_polynom.h>
#include <modq.h>

namespace val
{


integer degw(const s_expo& X,const val::vector<integer> &w)
{
    int i,n=s_expo::getdim();
    integer d;

    for (i=0;i<n;i++) d+= integer(X[i]) * w(i);
    return d;
}



template <>
int w_polynom<integer>::reduction(const Glist<w_polynom<integer> >& P,int top,int interred)
{
 int reduzbar,reduced=0,interreduced=0;
 unsigned maxnG=0-1,i;
 integer lcm,a1,gcd,div,hgcd,b,eins(1),minuseins(-1);
 s_expo X;
 w_polynom<integer> g,h;
 term *p,*r;
 GlistIterator<w_polynom<integer> > ItP;

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
            if (head->X<ItP.getelement().head->X){  
                maxnG=i;
                break;
            }
        }
        if (reduzbar) {
            nreduction++;
            g.head=ItP.getelement().head;
            hgcd=ggTspez(head->coeff,g.head->coeff);
            X=head->X/g.head->X;
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
            s_polynom<integer>::minusmalmonom(g,b,X);
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
         X=h.head->X/g.head->X;
         a1=EDIV(g.head->coeff,hgcd);
         b=EDIV(h.head->coeff,hgcd);
         r=h.head;
         h.head=h.head->next;
         delete r;
         p->next=h.head;
         if (a1 != eins) {operator*=(a1); gcd*=a1;}
		 g.head=g.head->next;
		 h.s_polynom<integer>::minusmalmonom(g,b,X);
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
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


template <>
int w_polynom<modq>::reduction(const Glist<w_polynom<modq> >& P,int top,int interred)
{
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modq b;
 s_expo X;
 w_polynom<modq> g,h;
 term *p,*r;
 GlistIterator<w_polynom<modq> > ItP;


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
            s_polynom<modq>::minusmalmonom(g,b,X);
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
		 h.s_polynom<modq>::minusmalmonom(g,b,X);
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
void w_polynom<modq>::liftpolynom(const Glist<w_polynom<modq> > &inwG,const Glist<w_polynom<modq> > &G,val::vector<char> &remain)
{
 int reduzbar,i,ismonom=0,red=0;
 integer d;
 s_expo X;
 modq b;
 w_polynom<modq> g,F;
 term *p;
 GlistIterator<w_polynom<modq> > pinwG,pG;

 if (head==NULL) return;
 if (inwG.isempty()) return;

 if (head->next==NULL) ismonom=1;

 // Only top-reductions are necessary
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
			 head = pG().head;
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
	 Error::error("\nPolynomial does not reduce to zero by lifting!");
 }
 head=F.head;
 F.head=NULL;
 g.head=NULL;
}



template<>
void w_polynom<integer>::liftpolynom(const Glist<w_polynom<integer> > &inwG,const Glist<w_polynom<integer> > &G,val::vector<char> &remain)
{
 int reduzbar,i,ismonom=0,red=0,notfirst=0;
 integer d;
 integer a1,hgcd,b,eins(1);
 s_expo X;
 w_polynom<integer> g,F;
 term *p;
 GlistIterator<w_polynom<integer> > ptoinwH,ptoH;

 if (head==NULL) return;
 if (inwG.isempty()) return;

 if (head->next==NULL) ismonom=1;
 // Only top-reductions are necessary
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
			 head = ptoH().head;
			 return;
		 }
		 red=1;
		 nreduction++;
		 g.head=ptoinwH().head;
		 hgcd=ggTspez(head->coeff,g.head->coeff);
		 X=head->X/g.head->X;

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

		 d=head->deg - g.head->deg;
		 g.head=g.head->next;
		 p=head;
		 head=head->next;
		 delete p;
		 minusmalmonom(g,b,X,d);
		 notfirst=1;
		 b.changesign();
		 F.minusmalmonom(ptoH(),b,X,d);
		 if (head==NULL) break;
	 }
 }
 while (reduzbar);

 if (head!=NULL) {
	 Error::error("\nPolynomial does not reduce to zero by lifting!");
 }
 F.normalize();
 head=F.head;
 F.head=NULL;
 g.head=NULL;
}

} //end namespace

