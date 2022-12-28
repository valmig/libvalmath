
#include <s_groebner/minpol.h>
#include <modint.h>

namespace modint_minpol
{

int get_q(const val::Glist<val::s_polynom<val::modint>>& G)
{
	int q=2;
	for (const auto& f : G) {
		for (const auto& pf : f) {
			q = pf.actualcoef().get_q();
			if (int (pf.actualcoef()) != 1) return q;
		}
	}
	return q;
}

} // end namespace modint_minpol


namespace minpol
{

template <>
void m_polynom<val::modint>::minusmalconst(const m_polynom<val::modint> &g,const val::modint &a)
{
 using namespace val;
 int c;
 modint b,a1=-a;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 c=expocompare(q->X,head->X);
	 if (c==-1) break;
	 else if (c==1) {
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
int m_polynom<val::modint>::linearreduction(const val::Glist<m_polynom<val::modint> >& P,const val::Glist<val::pol<val::modint> > &H,val::pol<val::modint> &minpol)
{
 using namespace val;
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modint b;
 s_expo X(0);
 m_polynom<modint> g,h;
 term *p,*r;
 GlistIterator<m_polynom<modint> > ItP;

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

 // 2. step inter-reductions
 p=head;

 while (p->next!=NULL) {
	 reduzbar=0;

	 for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
		 if (ItP.getelement().head->X==p->next->X){//(divisible(Y,X,m.X)){ //(X|Y) {
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
int m_polynom<val::modint>::normalform(const val::Glist<m_polynom<val::modint> > &P,val::modint&)
{
 using namespace val;
 int reduzbar,reduced=0;
 unsigned maxnG=0-1,i;
 modint b;
 s_expo X;
 m_polynom<modint> g,h;
 term *p,*r;
 GlistIterator<m_polynom<modint> > ItP;

 if (head==NULL) return 0;
 if (P.isempty()) {
	 return 0;
 }

    do {
        reduzbar=0;

        for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
            if ((ItP.getelement().head->X)|head->X) {
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

 if (head==NULL) {
	 g.head=NULL;
	 return reduced;
 }

 // 2. step inter-reductions:
 p=head;

 while (p->next!=NULL) {
	 reduzbar=0;

	 for (ItP.settohead(P),i=0;ItP.actualvalid() && i<maxnG;ItP.moveactual(),i++) {
		 if (ItP.getelement().head->X|p->next->X){
			 reduzbar=reduced=1;
			 break;
		 }
         if (p->next->X < ItP.getelement().head->X){  // X<Y
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
 h.head=NULL;
 g.head=NULL;
 return reduced;
}

val::pol<val::modint> minimalpolynom(const val::Glist<m_polynom<val::modint> > &G,int k,int q,int comment=1)
{
    using namespace val;
    int degXk,i,n=s_expo::getdim(),maxdeg=0;
    Glist<m_polynom<modint> > GredTerm;
    Glist<pol<modint> > H;
    s_expo X(0);
    modint c(1);
    m_polynom<modint> g,h;
    pol<modint> minpol;

    if (k<0 || k>=n ) k=n-1;

    // Check first if G is zero dimensional:
    val::vector<int> degX(0,n);
    if (!iszerodimensional(G,degX)) {
        common_bb::MyMessage("\nIdeal is not zero-dimensional!");
        maxdeg=5000;
    }

    degXk = degX[k];

    // Preparation:
    for (i=0;i<degXk;i++) {
        g.del();
        X[k] = i;
        g.insert(modint(1,q),X);
        GredTerm.inserttoend(g);
        H.inserttoend(pol<modint>(modint(1,q),i));
    }
    i=degXk;
    g.del();
    X[k]=i;
    g.insert(modint(1,q),X);

    X[k]=1;

    do {
        if (comment) {common_bb::Clear(); common_bb::WriteText("\nActual degree: " + ToString(i) + "\nComputing normal-form!");}
        g.normalform(G,c);
        if (g.iszero()) {
            if (comment) common_bb::WriteText("\n"+ToString(pol<modint>(modint(1,q),i)));
            GredTerm.dellist();
            return pol<modint>(modint(1,q),i);
        } // else:
        if (comment) common_bb::WriteText("\nComputing linear-reduction!");
        h=g;
        h.linearreduction(GredTerm,H,minpol);
        minpol=pol<modint>(modint(1,q),i)-minpol;
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

val::pol<val::modint> minimalpolynom(const val::Glist<val::s_polynom<val::modint> > &G,int k,int q = 0,int comment=1)
{
    using namespace val;
    Glist<m_polynom<modint> > H;
    if (!q) q = modint_minpol::get_q(G);
    convert(G,H);
    return minimalpolynom(H,k,q,comment);
}

} // end namespace minpol
