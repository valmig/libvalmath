#ifndef W_POLYNOM_H_INCLUDED
#define W_POLYNOM_H_INCLUDED

#include <s_polynom.h>
#include <Glist.h>

namespace val
{


template <class T>
class w_polynom : public s_polynom<T>
{
public:
    using s_polynom<T>::s_polynom;
    const w_polynom<T>& operator =(const s_polynom<T> &f) {s_polynom<T>::operator =(f);return *this;}
    const w_polynom<T>& operator =(val::s_polynom<T> &&f) {s_polynom<T>::operator =(std::move(f));return *this;}
    void moveto(s_polynom<T>& f) {s_polynom<T> *g=this; f=std::move(*g);}
    //
    static int wexpocompare(const s_expo& X1,const integer& d1,const s_expo& X2,const integer &d2);
    //
    // Prem: Leading monomial belongs to initial-form, f=h=0;
    // f keeps degrees, while the degree of monomials in h is set to 0.
    void makeinitforms(w_polynom<T> &f,w_polynom<T> &h,const val::vector<integer> &w) const;
    //
    void minusmalmonom(const w_polynom<T>& g,const T& a,const s_expo& Y,const integer& d);
    void liftpolynom(const Glist<w_polynom<T> > &inwG,const Glist<w_polynom<T> > &G,val::vector<char> &remain);
    int reduction(const Glist<w_polynom<T> >& P,int top=1,int interred=0);
};



integer degw(const s_expo& X,const val::vector<integer> &w);


//
template <class T>
void convertanddelete(Glist<s_polynom<T> > &G,Glist<w_polynom<T> > &H);

template <class T>
void convertanddelete(Glist<w_polynom<T> > &G,Glist<s_polynom<T> > &H);
//




//
template <class T>
int w_polynom<T>::wexpocompare(const s_expo& X1,const integer& d1,const s_expo& X2,const integer &d2)
{
 if (d1<d2) return -1;
 else if (d2<d1) return 1;
 else return expocompare(X1,X2);
}

//

template <class T>
void w_polynom<T>::makeinitforms(w_polynom<T> &f,w_polynom<T> &h,const val::vector<integer> &w) const
{
 integer d;
 struct w_polynom<T>::term *p,*rf,*qf,*rh,*qh;

 if (this->head==NULL) return;

 f.head=new struct w_polynom<T>::term(this->head->coeff,this->head->X,this->head->degree);
 h.head=new struct w_polynom<T>::term(this->head->coeff,this->head->X);

 qf=f.head;
 qh=h.head;
 d=degw(this->head->X,w);
 for (p=this->head->next;p!=NULL;p=p->next) {
	 if (degw(p->X,w)==d) {
		 rf=new struct w_polynom<T>::term(p->coeff,p->X,p->degree);
		 qf->next=rf;
		 qf=qf->next;
		 rh=new struct w_polynom<T>::term(p->coeff,p->X);
		 qh->next=rh;
		 qh=qh->next;
	 }
 }
}

//

template <class T>
void w_polynom<T>::minusmalmonom(const w_polynom<T> &g,const T &a,const s_expo &Y,const integer& d)
{
 int c;
 integer d1;
 T b,zero=w_polynom<T>::zero;
 s_expo Z;
 struct w_polynom<T>::term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (this->head==NULL) break;
	 Z=q->X*Y;
	 d1=q->deg +d;
	 c=expocompare(Z,this->head->X);
	 if (c==-1) break;
	 else if (c==1) {			// head->X < Z
		 b=q->coeff*a;
		 w_polynom<T>::changesign(b);
		 r=new struct w_polynom<T>::term(std::move(b),Z,std::move(d1),this->head);
		 this->head=r;
		 q=q->next;
		 break;
	 }
	 else {
         b=q->coeff*a;
         w_polynom<T>::changesign(b);
		 this->head->coeff+=b; 
		 if (this->head->coeff==zero) {
			 r=this->head;
			 this->head=this->head->next;
			 delete r;
		 }
		 else { q=q->next; break; }
	 }
 }
 p=this->head;
 if (p!=NULL) {
	 // Invariant: p->X > q->X :
	 while (q!=NULL && p->next!=NULL) {
		 Z=q->X*Y;
		 d1=q->deg+d;
		 c=expocompare(Z,p->next->X);
		 if (c==-1) p=p->next;
		 else if (c==1) {			// p->next->X < Z
			 b=q->coeff*a;
			 w_polynom<T>::changesign(b);
			 r=new  struct w_polynom<T>::term(std::move(b),Z,std::move(d),p->next);
			 p->next=r;
			 p=p->next;
			 q=q->next;
		 }
		 else {
             b=q->coeff*a;
             w_polynom<T>::changesign(b);
			 p->next->coeff+=b;//h
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
		 b=q->coeff*a;
		 w_polynom<T>::changesign(b);
		 if (b!=zero) {
			 this->head = new struct w_polynom<T>::term(std::move(b),q->X*Y,q->deg+d);
			 p=this->head;
			 q=q->next;
			 break;
		 }
	 }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff*a;
     w_polynom<T>::changesign(b);
	 if (b!=zero) {
		 r = new struct w_polynom<T>::term(std::move(b),q->X*Y,q->deg+d);
		 p->next=r;
		 p=p->next;
	 }
 }

 return;
}
//


template <class T>
void convertanddelete(Glist<s_polynom<T> > &G,Glist<w_polynom<T> > &H)
{
    w_polynom<T> f;
    if (!H.isempty()) H.dellist();
    for (G.resetactual();G.actualvalid();) {
        f = std::move(G.actualvalue());
        G.skiphead();
        H.inserttoend(std::move(f));
    }
}

template <class T>
void convertanddelete(Glist<w_polynom<T> > &G,Glist<s_polynom<T> > &H)
{
    s_polynom<T> f;
    if (!H.isempty()) H.dellist();
    for (G.resetactual();G.actualvalid();) {
        G.actualvalue().moveto(f);
        G.skiphead();
        H.inserttoend(std::move(f));
    }
}

//

} // end namespace

#endif // W_POLYNOM_H_INCLUDED
