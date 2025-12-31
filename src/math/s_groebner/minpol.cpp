
#include <s_groebner/minpol.h>


namespace minpol
{


template <>
val::integer m_polynom<val::integer>::normalize_cont()
{
    if (head==NULL) return val::integer(1);
    val::integer ggT=head->coeff;
    if (ggT==1 || ggT==-1) return val::integer(1);
    term *p=head->next;

    for (;p!=NULL;p=p->next) {
        ggT=val::ggTspez(ggT,p->coeff);
        if (ggT==1 || ggT==-1) return val::integer(1);
    }
    for (p=head;p!=NULL;p=p->next) p->coeff.EDIVBY(ggT);

    return ggT;
}




template <>
void m_polynom<val::integer>::minusmalconst(const m_polynom<val::integer> &g,const val::integer &a)
{
 using namespace val;
 int c;
 integer b;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Se head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 c=expocompare(q->X,head->X);
	 if (c==-1) break;
	 else if (c==1) {
		 b=q->coeff*a;
		 b.changesign();
		 r=new term(std::move(b),q->X,head);
		 head=r;
		 q=q->next;
		 break;
	 }
	 else {
         b=q->coeff*a;
         b.changesign();
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
	 // Invariant p->X > q->X :
	 while (q!=NULL && p->next!=NULL) {
		 c=expocompare(q->X,p->next->X);
		 if (c==-1) p=p->next;
		 else if (c==1) {
			 b=q->coeff*a;
			 b.changesign();
			 r=new term(std::move(b),q->X,p->next);
			 p->next=r;
			 p=p->next;
			 q=q->next;
		 }
		 else {
             b=q->coeff*a;
             b.changesign();
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

 // Insert evt. rest:
 if (head==NULL) {
	 for(;q!=NULL;q=q->next) {
		 b=q->coeff*a;
		 b.changesign();
		 if (b!=zero) {
			 head = new term(std::move(b),q->X);
			 p=head;
			 q=q->next;
			 break;
		 }
	 }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff*a;
     b.changesign();
	 if (b!=zero) {
		 r = new term(std::move(b),q->X);
		 p->next=r;
		 p=p->next;
	 }
 }

 return;
}


template <>
void m_polynom<val::modq>::minusmalconst(const m_polynom<val::modq> &g,const val::modq &a)
{
 using namespace val;
 int c;
 modq b,a1=-a;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 c=expocompare(q->X,head->X);
	 if (c==-1) break;
	 else if (c==1) {//(head->dat<m1) {
		 b=q->coeff;
		 b*=a1;
		 r=new term(b,q->X,head);
		 head=r;
		 q=q->next;
		 break;
	 }
	 else {
         b=q->coeff;
         b*=a1;
		 head->coeff+=b; //h
		 if (head->coeff==zero) {//dat==NULL) {
			 r=head;
			 head=head->next;
			 delete r;
		 }
		 else { q=q->next; break; }
	 }
 }
 p=head;
 if (p!=NULL) {
	 // Invariant p->X > q->X :
	 while (q!=NULL && p->next!=NULL) {
		 c=expocompare(q->X,p->next->X);
		 if (c==-1) p=p->next;
		 else if (c==1) {
			 b=q->coeff;
			 b*=a1;
			 r=new term(b,q->X,p->next);
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

 // Insert evt. rest :
 if (head==NULL) {
	 for(;q!=NULL;q=q->next) {
		 b=q->coeff;
		 b*=a1;
		 if (b!=zero) {
			 head = new term(b,q->X);
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
		 r = new term(b,q->X);
		 p->next=r;
		 p=p->next;
	 }
 }

 return;
}



template <>
template <>
int m_polynom<val::integer>::linearreduction(const val::Glist<m_polynom<val::integer> >& P,const val::Glist<val::pol<val::rational> > &H,val::pol<val::rational> &minpol)
{
 using namespace val;
 int reduzbar,reduced=0,interreduced=0;
 unsigned maxnG=0-1,i;
 integer lcm,a1,gcd,div,hgcd,b,eins(1),minuseins(-1);
 //expo X,Y;
 s_expo X(0);
 m_polynom<integer> g,h;
 term *p,*r;
 GlistIterator<m_polynom<integer> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {return 0;}

 do {
        reduzbar=0;
        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)==head->X) {
                reduzbar=reduced=1;
                break;
            }
        }
        if (reduzbar) {
            nreduction++;
            g.head=ItP.getelement().head;
            hgcd=ggTspez(head->coeff,g.head->coeff);
            if (hgcd.isunit()) {
                b=std::move(head->coeff);
                p=head;
                head = head->next;
                delete p;
                if (g.head->coeff!=eins) {operator *=(g.head->coeff);minpol*=rational(g.head->coeff);}
            }
            else {
                a1=EDIV(g.head->coeff,hgcd);
                b=EDIV(head->coeff,hgcd);
                p=head;
                head=head->next;
                delete p;
                if (a1!=eins) {operator*=(a1);minpol*=rational(a1);}
            }
            g.head=g.head->next;
            minusmalconst(g,b);
            minpol-= rational(b)*H[i];
            minpol*=rational(integer(1),normalize_cont());
            if (head==NULL) break;
        }
    }
 while (reduzbar);
 if (head==NULL) {
	 g.head=NULL;
	 return reduced;
 }

 // 2. Inter-reductions;
 p=head;

 while (p->next!=NULL) {
	 reduzbar=0;
	 for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
		 if (ItP.getelement().head->X==p->next->X) {
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
         a1=EDIV(g.head->coeff,hgcd);
         b=EDIV(h.head->coeff,hgcd);
         r=h.head;
         h.head=h.head->next;
         delete r;
         p->next=h.head;
         if (a1 != eins) {operator*=(a1); gcd*=a1;minpol*=rational(a1);}
         
		 g.head=g.head->next;
		 h.minusmalconst(g,b);
		 minpol-=rational(b) * H[i];
		 div=gcd;
		 for (r=h.head;r!=NULL;r=r->next) {
			 if (div.isunit()) break;
			 else (div=ggTspez(div,r->coeff));
		 }

		 p->next=h.head;
		 if (!div.isunit()) {gcd.EDIVBY(div);edivby(div); minpol*=rational(integer(1),div); }
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
template <>
int m_polynom<val::modq>::linearreduction(const val::Glist<m_polynom<val::modq> >& P,const val::Glist<val::pol<val::modq> > &H,val::pol<val::modq> &minpol)
{
 using namespace val;
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modq b;
 s_expo X(0);
 m_polynom<modq> g,h;
 term *p,*r;
 GlistIterator<m_polynom<modq> > ItP;


 minpol.del();
 if (head==NULL) return 0;
 if (P.isempty()) {
	 normalize();
	 return 0;
 }

    do {
        reduzbar=0;

        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)==head->X) {
			 reduzbar=reduced=1;
			 break;
            }
        }
        if (reduzbar) {
		    nreduction++;
            g.head=ItP.getelement().head;
            b=head->coeff/g.head->coeff;
            g.head=g.head->next;
            p=head;
            head=head->next;
            delete p;
            minpol+=b*H[i];
            minusmalconst(g,b);
            if (head==NULL) break;
        }
    }
    while (reduzbar);

 if (head==NULL) {
	 g.head=NULL;
	 return reduced;
 }

 // 2.Inter-reductions
 p=head;

 while (p->next!=NULL) {
	 reduzbar=0;

	 for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
		 if (ItP.getelement().head->X==p->next->X){
			 reduzbar=reduced=1;
			 break;
		 }
	 }
	 if (reduzbar) {
		 nreduction++;
		 h.head=p->next;
		 g.head=ItP.getelement().head;
		 b=h.head->coeff/g.head->coeff;
		 p->next=NULL;
		 g.head=g.head->next;
		 r=h.head;
		 h.head=h.head->next;
		 delete r;
		 h.minusmalconst(g,b);
		 minpol+=b*H[i];
		 p->next=h.head;
	 }
	 else p=p->next;
 }
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


template<>
template<>
int m_polynom<val::integer>::normalform(const val::Glist<m_polynom<val::integer> > &P,val::rational &c)
{
 using namespace val;
 int reduzbar,reduced=0,interreduced=0;
 unsigned maxnG=0-1,i;
 integer lcm,a1,gcd,div,hgcd,b,eins(1),minuseins(-1);
 s_expo X;
 m_polynom<integer> g,h;
 term *p,*r;
 GlistIterator<m_polynom<integer> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {return 0;}



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
             if (g.head->coeff!=eins) {operator *=(g.head->coeff);c*=g.head->coeff;}
         }
         else {
             a1=EDIV(g.head->coeff,hgcd);
             b=EDIV(head->coeff,hgcd);
             p=head;
             head=head->next;
             delete p;
             if (a1!=eins) {operator*=(a1);c*=a1;}
         }
         g.head=g.head->next;
         minusmalmonom(g,b,X);
         c/=rational(normalize_cont());
         if (head==NULL) break;
     }
 }
 while (reduzbar);

 if (head==NULL) {
	 g.head=NULL;
	 return reduced;
 }

 // 2.Schritt Interrreduktion;
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
         if (a1 != eins) {operator*=(a1); gcd*=a1;c*=a1;}
		 g.head=g.head->next;
		 h.minusmalmonom(g,b,X);
		 div=gcd;
		 for (r=h.head;r!=NULL;r=r->next) {
			 if (div.isunit()) break;
			 else (div=ggTspez(div,r->coeff));
		 }
		 p->next=h.head;
		 if (!div.isunit()) {gcd.EDIVBY(div);edivby(div);c/=rational(div);}
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

template<>
template<>
int m_polynom<val::modq>::normalform(const val::Glist<m_polynom<val::modq> > &P,val::modq&)
{
 using namespace val;
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modq b;
 s_expo X;
 m_polynom<modq> g,h;
 term *p,*r;
 GlistIterator<m_polynom<modq> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {
	 return 0;
 }

 do {
	 reduzbar=0;
     for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
         if ((ItP.getelement().head->X)|head->X) {//(divisible(Y,X,m.X)) {//(X|Y) {
			 reduzbar=reduced=1;
			 break;
         }
		 if (head->X<ItP.getelement().head->X){  // X<Y
			 maxnG=i;
			 break;
         }

     }
     if (reduzbar) {
		 nreduction++;
         g.head=ItP.getelement().head;
         b=head->coeff/g.head->coeff;
         X=head->X/ItP.getelement().head->X;  // Y/X
         g.head=g.head->next;
         p=head;
         head=head->next;
         delete p;
         minusmalmonom(g,b,X);
         if (head==NULL) break;
     }
 }
 while (reduzbar);

 if (head==NULL) {
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
		 X=p->next->X/g.head->X;  // Y/X
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
 h.head=NULL;
 g.head=NULL;
 return reduced;
}


val::pol<val::rational> minimalpolynom(const val::Glist<m_polynom<val::integer> > &G,int k,int comment)
{
    using namespace val;
    int degXk,i,n=s_expo::getdim(),maxdeg=0;
    Glist<m_polynom<integer> > GredTerm;
    Glist<pol<rational> > H;
    s_expo X(0);
    rational c(1);
    m_polynom<integer> g,h;
    pol<rational> minpol;

    if (k<0 || k>=n) k=n-1;

    // Check if G is zero-dimensional
    val::vector<int> degX(0,n);
    if (!iszerodimensional(G,degX)) {
        common_bb::MyMessage("\nIdeal is not zero-dimensional!");
        maxdeg=5000;
    }

    degXk = degX[k];
    // ---------------------------------

    // Preparation:
    for (i=0;i<degXk;i++) {
        g.del();
        X[k] = i;
        g.insert(integer(1),X);
        GredTerm.inserttoend(g);
        H.inserttoend(pol<rational>(rational(1),i));
    }
    i=degXk;
    g.del();
    X[k]=i;
    g.insert(integer(1),X);

    X[k]=1;

    do {
        if (comment) {common_bb::Clear();common_bb::WriteText("\nActual degree: " + ToString(i) + "\nComputing normal-form!");}
        g.normalform(G,c);
        if (g.iszero()) {
            if (comment) common_bb::WriteText("\n"+ToString(pol<rational>(1,i)));
            GredTerm.dellist();
            return pol<rational>(1,i);
        } // else:
        minpol=pol<rational>(c,i);
        if (comment) common_bb::WriteText("\nComputing linear reductions!");
        h=g;
        h.linearreduction(GredTerm,H,minpol);
        if (!h.iszero()) {
            GredTerm.inserttoend(std::move(h));
            H.inserttoend(std::move(minpol));
        }
        else {
            minpol.normalize();
            if (comment) common_bb::WriteText("\n" + ToString(minpol));
            break;
        }
        i++;
        g*=X;
        if (maxdeg && i>maxdeg) {minpol.del();break;}
    }
    while (!g.iszero());

    GredTerm.dellist();
    return minpol;
}


val::pol<val::modq> minimalpolynom(const val::Glist<m_polynom<val::modq> > &G,int k,int comment)
{
    using namespace val;
    int degXk,i,n=s_expo::getdim(),maxdeg=0;
    Glist<m_polynom<modq> > GredTerm;
    Glist<pol<modq> > H;
    s_expo X(0);
    modq c(1);
    m_polynom<modq> g,h;
    pol<modq> minpol;

    if (k<0 || k>=n ) k=n-1;

    // Check if G is zero-dimensional:
    val::vector<int> degX(0,n);
    if (!iszerodimensional(G,degX)) {
        common_bb::MyMessage("\nIdeal is not zero-dimensional!");
        maxdeg=5000;
    }

    degXk = degX[k];
    // ---------------------------------
    // Preparation:
    for (i=0;i<degXk;i++) {
        g.del();
        X[k] = i;
        g.insert(modq(1),X);
        GredTerm.inserttoend(g);
        H.inserttoend(pol<modq>(modq(1),i));
    }
    i=degXk;
    g.del();
    X[k]=i;
    g.insert(modq(1),X);
    X[k]=1;

    do {
        if (comment) {common_bb::Clear(); common_bb::WriteText("\nActual degree: " + ToString(i) + "\nComputing normal-form!");}
        g.normalform(G,c);
        if (g.iszero()) {
            if (comment) common_bb::WriteText("\n"+ToString(pol<modq>(1,i)));
            GredTerm.dellist();
            return pol<modq>(1,i);
        } // else:
        if (comment) common_bb::WriteText("\nComputing linear reductions!");
        h=g;
        h.linearreduction(GredTerm,H,minpol);
        minpol=pol<modq>(modq(1),i)-minpol;
        if (!h.iszero()) {
            GredTerm.inserttoend(std::move(h));
            H.inserttoend(std::move(minpol));
        }
        else {
            minpol.normalize();
            if (comment) common_bb::WriteText("\n" + ToString(minpol));
            break;
        }
        i++;
        g*=X;
        if (maxdeg && i>maxdeg) {minpol.del(); break;}
    }
    while (!g.iszero());

    GredTerm.dellist();
    return minpol;
}

val::pol<val::rational> minimalpolynom(const val::Glist<val::s_polynom<val::integer> > &G,int k,int comment)
{
    using namespace val;
    Glist<m_polynom<integer> > H;
    convert(G,H);
    return minimalpolynom(H,k,comment);
}


val::pol<val::modq> minimalpolynom(const val::Glist<val::s_polynom<val::modq> > &G,int k,int comment)
{
    using namespace val;
    Glist<m_polynom<modq> > H;
    convert(G,H);
    return minimalpolynom(H,k,comment);
}


} //end namespace
