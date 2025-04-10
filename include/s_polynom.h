#ifndef S_POLYNOM_H_INCLUDED
#define S_POLYNOM_H_INCLUDED

//#include <matrix.h>
#include <s_expo.h>
#include <atomic>
#include <integer.h>
#include <vector.h>




namespace val
{
class modq; class integer; class rational;

template <class T> class s_polynomIterator;
template <class T> class s_polynom;
template <class T> class Glist;
template <class T> class d_array;




//
template <class T> std::ostream& operator <<(std::ostream&,const s_polynom<T>&);
template <class T> std::istream& operator >>(std::istream&,s_polynom<T>&);

//


// class for multivariate polynomials with multivariate exponents (terms) of the class s_expo (static dimension).
// Data structure is a simple linked list, ordered by the order of the terms (s_expo)

template <class T>
class s_polynom
{
protected:
    struct term{
        T coeff;
        s_expo X;
        integer deg;
        term *next;
        //
        term() = delete;
        term(const T &wert,term* nex=NULL) : coeff(wert), X(0), next(nex) {s_polynom<T>::mnumber++;}
        term(T&& wert,term* nex=NULL) : coeff(std::move(wert)),X(0), next(nex) {s_polynom<T>::mnumber++;}
        term(const T& wert,const s_expo& Z,term* nex=NULL) : coeff(wert), X(Z) , next(nex) {s_polynom<T>::mnumber++;}
        term(const T& wert,const s_expo& Z,const integer &d,term* nex=NULL) : coeff(wert), X(Z) ,deg(d), next(nex) {s_polynom<T>::mnumber++;}
        term(T&& wert,const s_expo& Z,integer&& d,term* nex=NULL) : coeff(std::move(wert)),X(Z),deg(std::move(d)),next(nex) {s_polynom<T>::mnumber++;}
        term(T&& wert,const s_expo& Z,term* nex=NULL) : coeff(std::move(wert)), X(Z) , next(nex) {s_polynom<T>::mnumber++;}
        term(const T& wert,s_expo&& Z, term* nex = NULL) : coeff(wert),X(std::move(Z)),next(nex) {s_polynom<T>::mnumber++;}
        term(T&& wert,s_expo&& Z,term* nex=NULL) : coeff(std::move(wert)), X(std::move(Z)) , next(nex) {s_polynom<T>::mnumber++;}
        term(T&& wert,s_expo&& Z,integer&& d,term *nex=NULL) :coeff(std::move(wert)),X(std::move(Z)),deg(std::move(d)),next(nex) {s_polynom<T>::mnumber++;}
        ~term() {s_polynom<T>::mnumber--;}
    };

    term* head = nullptr;
    //
    static std::atomic<int> mnumber;
    static const T zero;
    static s_expo neutral;
    static const integer zerodegree;
    //
    static void changesign(T&);
    static int bigger(term *p1,term *p2,int degoption);
    static void merge(term* &basis1,term* &basis2,int degoption);
    static void sort(term* &basis,int l,int degoption);
    //
    s_polynom<T> add(const s_polynom<T>& p,int plus=1) const;
    void addto(const s_polynom<T>& t,int plus=1);
public:
    // Constructors:
    s_polynom() = default;
    s_polynom(const T& wert) {if (wert==zero) head=NULL; else head= new term(wert); }
    s_polynom(T&& wert) {if (wert==zero) head=NULL; else head= new term(std::move(wert)); }
    s_polynom(const T& wert,const s_expo& Z) {if (wert==zero) head=NULL; else head= new term(wert,Z); }
    s_polynom(T&& wert,const s_expo& Z) {if (wert==zero) head=NULL; else head=new term(std::move(wert),Z); }
    s_polynom(T&& wert,s_expo&& Z) {if (wert==zero) head=NULL; else head= new term(std::move(wert),std::move(Z)); }
    s_polynom(const s_polynom<T>&);
    s_polynom(s_polynom<T>&& f) {head=f.head;f.head=NULL;}
    // Destructors:
    void del();
    ~s_polynom() {del();}
    //
    template <class S> const s_polynom<T>& convertfrom(const s_polynom<S>&);
    const s_polynom<T>& homogenize();
    //
    int iszero() const {return (head==NULL);}
    int ismonomial() const {if (head==NULL) return 1;if (head->next==NULL) return 1; else return 0;}
    //
    const s_polynom<T>& operator =(const s_polynom<T>&);
    const s_polynom<T>& operator =(s_polynom<T>&& f) {if (head!=f.head)
                                                        {s_polynom<T> g;g.head=head;head=f.head;f.head=NULL;}
                                                      return *this;}
    const s_polynom<T>& operator =(const T& wert) {s_polynom<T> f,g(wert);f.head=head;head=g.head;g.head=NULL;return *this;}
    const s_polynom<T>& operator =(T&& wert) {s_polynom<T> f,g(std::move(wert));f.head=head;head=g.head;g.head=NULL;return *this;}
    //
    int operator ==(const s_polynom<T>&) const;
    int operator ==(const T&) const;
    int operator !=(const s_polynom<T>& f) const {return !(*this==f);}
    int operator !=(const T& a) const {return !(*this==a);}
    int operator <(const s_polynom<T>&) const;
    //
    const s_polynom<T>& insert(const T&,const s_expo&);  // ordered
    const s_polynom<T>& insert(T&&,const s_expo&);       // ordered
    const s_polynom<T>& insert(T&&,s_expo&&);            // ordered
    const s_polynom<T>& reord(int degoption=0) {sort(head,length(),degoption);return *this;}
    const s_polynom<T>& inserttohead(const T& coeff,const s_expo &X) {if (coeff!=zero) head=new term(coeff,X,head); return *this;}
    const s_polynom<T>& inserttohead(T&& coeff,const s_expo &X) {if (coeff!=zero) head = new term(std::move(coeff),X,head); return *this;}
    const s_polynom<T>& inserttohead(const T& coeff,const s_expo &X,const integer& deg) {if (coeff!=zero) head=new term(coeff,X,deg,head); return *this;}
    //
    int length() const;
    int totdeg() const;
    const T& LC() const {if (head==NULL) return zero; else return head->coeff;}
    const T& coefficient(const s_expo& X) const;
    const s_expo& LT() const {if (head==NULL) {return neutral=s_expo(0);} else return head->X;}
    const s_polynom<T>& setdegree(const vector<integer> &w);
    const s_polynom<T>& setzerodegree();
    const s_polynom<T>& setLTdegree(const vector<integer> &w);
    const s_polynom<T>& setLTdegree(const integer& d) {if (head==NULL) return *this; head->deg=d; return *this;}
    const s_polynom<T>& setLTzerodegree() {if (head==NULL) return *this; head->deg=integer(0);return *this;}
    const integer& getLTdegree() const {if (head==NULL) return zerodegree; else return head->deg;}
    void makeinitforms(s_polynom<T> &f,s_polynom<T> &h,const val::vector<integer> &w) const;
    //
    static int wexpocompare(const s_expo& X1,const integer& d1,const s_expo& X2,const integer &d2);
    static integer degw(const s_expo& X,const val::vector<integer> &w);
    //
    s_polynom<T> operator +(const s_polynom<T>& p) const {return add(p,1);}
    s_polynom<T> operator -(const s_polynom<T>& p) const {return add(p,0);}
    s_polynom<T> operator-() const;
    const s_polynom<T>& operator +=(const s_polynom<T>& f) {addto(f,1);return *this;}
    const s_polynom<T>& operator -=(const s_polynom<T>& f) {addto(f,0);return *this;}
    DLL_PUBLIC void minusmalmonom(const s_polynom<T>& g,const T& a,const s_expo& Y);                     // *this -= a*Y*g
    DLL_PUBLIC void minusmalmonom(const s_polynom<T>& g,const T& a,const s_expo& Y,const integer& d);
    //
    const s_polynom<T>& operator*=(const T&);
    const s_polynom<T>& operator*=(const s_expo&);
    DLL_PUBLIC void edivby(const T&);
    void malmonom(const T& a,const s_expo& Y);
    DLL_PUBLIC const s_polynom<T>& normalize();
    DLL_PUBLIC T content() const;
    //
    static std::atomic<int> nreduction;
    DLL_PUBLIC int reduction(const Glist<s_polynom<T> >& P,int top=1,int interred=0,int degoption=0);     // Reduction mod <P> for Buchberger-Algorithm
    DLL_PUBLIC int reduction(const d_array<s_polynom<T> >& P,int top=1,int interred=0,int degoption=0);
    DLL_PUBLIC void liftpolynom(const Glist<s_polynom<T> > &inwG,const Glist<s_polynom<T> > &G,val::vector<char> &remain); // Lift-function for groebner-walk-algorithm
    //
    static int getstaticexpodim() {return s_expo::dim;}
    static int getmnumber() {return int(mnumber);}
    //
    friend std::ostream& operator <<<T>(std::ostream&,const s_polynom&);
    friend std::istream& operator >><T>(std::istream&,s_polynom&);
    //
    //
    s_polynomIterator<T> begin() const {return s_polynomIterator<T>(*this);}
    s_polynomIterator<T> end() const {return s_polynomIterator<T>();}
    //
    //
    friend class s_polynomIterator<T>;
};



template <class T>
class s_polynomIterator
{
private:
    struct s_polynom<T>::term *actual;
public:
    s_polynomIterator() {actual=NULL;}
    s_polynomIterator(const s_polynom<T>& f) {actual=f.head;}
    const s_polynomIterator& operator = (const s_polynom<T>& f) {actual=f.head;return *this;}
    operator int() const {return (actual!=NULL);}
    void operator++(int) {if (actual!=NULL) actual=actual->next;}
    void operator++() {if (actual!=NULL) actual=actual->next;}
    const s_polynomIterator<T>& operator* () const {return *this;}
    const T& actualcoef() const {if (actual==NULL) return s_polynom<T>::zero;else return actual->coeff;}
    const s_expo& actualterm() const {if (actual==NULL) {return s_polynom<T>::neutral=s_expo(0);} else return actual->X;}
    const integer& actualdegree() const {if (actual==NULL) return s_polynom<T>::zerodegree; else return actual->deg;}
    int actualvalid() const {return (actual!=NULL);}
    int nextvalid() const {if (actual==NULL) return 0; if (actual->next==NULL) return 0; else return 1;}
    int moveactualtolast();   // zeigt nun auf letztes Monom
};



//
template <>
DLL_PUBLIC std::atomic<int> s_polynom<modq>::mnumber;

template <>
DLL_PUBLIC std::atomic<int> s_polynom<integer>::mnumber;

template <class T>
std::atomic<int> s_polynom<T>::mnumber(0);


template <> DLL_PUBLIC const val::modq s_polynom<val::modq>::zero;

template <class T>
const T s_polynom<T>::zero(val::zero_element<T>());



template <>
DLL_PUBLIC s_expo s_polynom<modq>::neutral;


template <class T>
s_expo s_polynom<T>::neutral(0);


template <class T>
const integer s_polynom<T>::zerodegree(0);


template <>
DLL_PUBLIC std::atomic<int> s_polynom<modq>::nreduction;

template <>
DLL_PUBLIC std::atomic<int> s_polynom<integer>::nreduction;

template <class T>
std::atomic<int> s_polynom<T>::nreduction(0);

//


template <>
DLL_PUBLIC void s_polynom<integer>::changesign(integer&);

template <>
DLL_PUBLIC void s_polynom<rational>::changesign(rational&);

template <class T>
void s_polynom<T>::changesign(T& a)
{
    a=-a;
}

//

template <class T>
int s_polynom<T>::wexpocompare(const s_expo& X1,const integer& d1,const s_expo& X2,const integer &d2)
{
 if (d1<d2) return -1;
 else if (d2<d1) return 1;
 else return expocompare(X1,X2);
}

template <class T>
integer s_polynom<T>::degw(const s_expo& X,const val::vector<integer> &w)
{
    int i,n=s_expo::getdim();
    integer d;

    for (i=0;i<n;i++) d+= integer(X[i]) * w(i);
    return d;
}

//

template <class T>
int s_polynom<T>::bigger(s_polynom<T>::term *p1,s_polynom<T>::term *p2,int degoption)
{
    if (degoption) {
        if (p1->deg > p2->deg) return 1;
        else if (p2->deg > p1->deg) return 0;
        else return (p1->X > p2->X);
    }
    else return (p1->X > p2->X);
}


// Merge two ordered lists. basis1 != NULL, basis2 != NULL and basis1->X != basis2->X
template <class T>
void s_polynom<T>::merge(term* &basis1,term* &basis2,int degoption)
{
 term *p1,*p2,*r;

 p1=basis1;p2=basis2;
 if (bigger(p2,basis1,degoption)) {      // p2->X > basis1->X
     while ((p2->next!=NULL) && bigger(p2->next,basis1,degoption)) p2=p2->next;   // p2->next->X > basis1->X
     basis1=basis2;
     basis2=p2->next;
     p2->next=p1;
 }
 // now p2<p1
 if (basis2==NULL) return;
 p2=basis2;
 while (p1->next!=NULL) {
     if ((basis2!=NULL) && bigger(basis2,p1->next,degoption) ) {                        // basis2->X > p1->next->X
         while ((p2->next!=NULL) && bigger(p2->next,p1->next,degoption) ) p2=p2->next;  // p2->next->X > p1->next->X
         r=p1->next;
         p1->next=basis2;
         basis2=p2->next;
         p2->next=r;
         p1=r;p2=basis2;;
     }
     else p1=p1->next;
 }
 p1->next=basis2;
 basis2=NULL;
 return;
}


template <class T>
void s_polynom<T>::sort(term* &basis,int l,int degoption)
{
 int l1,l2,i;
 term *p1,*p2;

 if (l<=1) return;
 l1=l/2; l2=l-l1;
 p1=basis;
 for (i=1;i<l1;i++) p1=p1->next;
 p2=p1->next;
 p1->next=NULL;
 p1=basis;
 sort(p1,l1,degoption);sort(p2,l2,degoption);
 merge(p1,p2,degoption);
 basis=p1;
 p1=p2=NULL;
 return;
}

//

template <class T>
s_polynom<T>::s_polynom(const s_polynom<T>& f)
{
    head=NULL;
    if (f.head==NULL) return;

    term *q,*t,*r;

    q=f.head;
    t=new term(q->coeff,q->X);
    head=t;
    while (q->next != NULL) {
        q=q->next;
        r=new term(q->coeff,q->X,q->deg);
        t->next=r;
        t=t->next;
    }
}

template <>
template <> DLL_PUBLIC const s_polynom<modq>& s_polynom<modq>::convertfrom(const s_polynom<integer>&);

template <class T>
template <class S> const s_polynom<T>& s_polynom<T>::convertfrom(const s_polynom<S> &f)
{
    if (head!=NULL) del();
    head=NULL;
    if (f.iszero()) return *this;

    T wert;

    s_polynomIterator<S> q;
    term *t,*r;

    q=f;
    while (q) {
        if ((wert=T(q.actualcoef()))!=zero) break;
        q++;
    }
    if (!q) return *this;
    t=new term(wert,q.actualterm());
    head=t;
    while (q.nextvalid()) {
        q++;
        if ((wert=T(q.actualcoef()))!=zero) {
            r=new term(wert,q.actualterm());
            t->next=r;
            t=t->next;
        }
    }
    return *this;
}


template <class T>
const s_polynom<T>& s_polynom<T>::homogenize()
{
    if (head==NULL || s_expo::dim==1) return *this;

    s_expo::dim--;

    int deg=totdeg();

    s_expo::dim++;

    for (term *p=head;p!=NULL;p=p->next) {
        p->X[s_expo::dim-1]=0;
        p->X[s_expo::dim-1] = deg - p->X.totdeg();
    }

    return *this;
}


//

template <class T>
void s_polynom<T>::del()
{
    term *p;

    while (head!=NULL) {
        p=head;
        head=head->next;
        delete p;
    }
}

//
template <class T>
const T& s_polynom<T>::coefficient(const s_expo &X) const
{
    for (term *p=head;p!=nullptr;p=p->next) {
        if (X>p->X) return zero;
        else if (X==p->X) return p->coeff;
    }
    return zero;
}

//
template <class T>
const s_polynom<T>& s_polynom<T>:: operator =(const s_polynom<T>& s)
{

 if (head==s.head) return *this;
 s_polynom<T> help(s);

 val::swap(head,help.head);

 return *this;
}

//

template <class T>
int s_polynom<T>:: operator==(const s_polynom<T>& f) const
{
 if (head==f.head) return *this;
 term *p,*q;

 p=head;
 q=f.head;

 while (p!=NULL) {
       if (q==NULL) return 0;
       if (p->coeff!=q->coeff) return 0;
       if (p->X!=q->X) return 0;
       p=p->next;
       q=q->next;
 }
 if (q==NULL) return 1;
 else return 0;
}


template <class T>
int s_polynom<T>:: operator ==(const T& i) const
{
 if (head==NULL && i==zero) return 1;
 else if (head->next==NULL && head->coeff==i) {
      for (int i=0;i<s_expo::dim;i++) if (head->X.coeff[i]) return 0;
      return 1;
 }
 else return 0;
}



template <class T>
int s_polynom<T>::operator <(const s_polynom<T>& g) const
{
 term *p,*q;

 for (q=g.head,p=head;q!=NULL;q=q->next,p=p->next) {
       if (p==NULL) return 1;
       if (p->deg > q->deg) return 0;
       if (q->deg > p->deg) return 1;
       if (p->X<q->X) return 1;
       else if (p->X==q->X) continue;
       else return 0;
 }
 return 0;
}

//

template <class T>
const s_polynom<T>& s_polynom<T>::insert(const T& wert,const s_expo& g)
{
 int stop=0;
 term *p,*q;

 if (wert==zero) return *this;
 if (head==NULL) {
    p=new term(wert,g);
    head=p;
    return *this;
 }
 else {
    if (g==head->X) {
       head->coeff+=wert;
       if (head->coeff!=zero)                    // coeff == 0
           return *this;                     // => delete term
       else {
           p=head;
           head=head->next;
           delete p;
           return *this;
       }
    }
    if (head->X<g) {
       p=new term(wert,g,head);
       head=p;
       return *this;
    }
    q=head;
    while ((q->next!=NULL) && !stop) {
         if ( g< q->next->X)
            q=q->next;
         else stop=1;
    }
    if (stop)
       if (g==q->next->X) {
            q->next->coeff+=wert;
            if (q->next->coeff!=zero)           // coeff == 0
                return *this;                   // => delete term
            else {
                p=q->next;
                q->next=p->next;
                delete p;
                return *this;
            }
       }
 }
 p= new term(wert,g,q->next);
 q->next=p;
 return *this;
}


template <class T>
const s_polynom<T>& s_polynom<T>::insert(T&& wert,const s_expo& g)
{
 int stop=0;
 term *p,*q;

 if (wert==zero) return *this;
 if (head==NULL) {
    p=new term(std::move(wert),g);
    head=p;
    return *this;
 }
 else {
    if (g==head->X) {
       head->coeff+=wert;
       if (head->coeff!=zero)                   // coeff == 0
           return *this;                    // => delete term
       else {
           p=head;
           head=head->next;
           delete p;
           return *this;
       }
    }
    if (head->X<g) {
       p=new term(std::move(wert),g,head);
       head=p;
       return *this;
    }
    q=head;
    while ((q->next!=NULL) && !stop) {
         if ( g< q->next->X)
            q=q->next;
         else stop=1;
    }
    if (stop)
       if (g==q->next->X) {
            q->next->coeff+=wert;
            if (q->next->coeff!=zero)           // coeff == 0
                return *this;                   // => delete term
            else {
                p=q->next;
                q->next=p->next;
                delete p;
                return *this;
            }
       }
 }
 p= new term(std::move(wert),g,q->next);
 q->next=p;
 return *this;
}


template <class T>
const s_polynom<T>& s_polynom<T>::insert(T&& wert,s_expo&& g)
{
 int stop=0;
 term *p,*q;

 if (wert==zero) return *this;
 if (head==NULL) {
    p=new term(std::move(wert),std::move(g));
    head=p;
    return *this;
 }
 else {
    if (g==head->X) {
       head->coeff+=wert;
       if (head->coeff!=zero)                   // coeff == 0
           return *this;                    // => delete term
       else {
           p=head;
           head=head->next;
           delete p;
           return *this;
       }
    }
    if (head->X<g) {
       p=new term(std::move(wert),std::move(g),head);
       head=p;
       return *this;
    }
    q=head;
    while ((q->next!=NULL) && !stop) {
         if ( g< q->next->X)
            q=q->next;
         else stop=1;
    }
    if (stop)
       if (g==q->next->X) {
            q->next->coeff+=wert;
            if (q->next->coeff!=zero)           // coeff == 0
                return *this;                   // => delete term
            else {
                p=q->next;
                q->next=p->next;
                delete p;
                return *this;
            }
       }
 }
 p= new term(std::move(wert),std::move(g),q->next);
 q->next=p;
 return *this;
}

//

//
template <class T>
s_polynom<T> s_polynom<T>::add(const s_polynom<T>& p,int plus) const
{
 s_polynom<T> q;                                     
 term *r,*s,*t,*u;                                   
 T z;


 r=head;
 s=p.head;
 do {                                    // Set first term in q
    if (r==NULL && s==NULL) return q;    
    if (r==NULL){                        
       if (plus) q.head = new term(s->coeff,s->X);
       else q.head= new term(-s->coeff,s->X);
       s=s->next;
    }
    else if (s==NULL){
       q.head=new term(r->coeff,r->X);
       r=r->next;
    }
    else { 
       if ( (r->X) > (s->X) ){
      q.head = new term(r->coeff,r->X);
      r=r->next;
       }
       else if ( (r->X) == (s->X) ) {
            if (plus) z = r->coeff + s->coeff;
            else z = r->coeff - s->coeff;
            if (z!=zero) q.head = new term(std::move(z),r->X);
            r=r->next;
            s=s->next;
       }
       else {
            if (plus) q.head = new term(s->coeff,s->X);
            else q.head = new term(-s->coeff,s->X);
            s=s->next;
       }
    }
  }
 while (q.head==NULL);

 t=q.head;
 while (r!=NULL && s!=NULL) {
       if ( (r->X) > (s->X) ) {
          u= new term(r->coeff,r->X);
          t->next=u;
          r=r->next;
          t=t->next;
       }
       else {
            if ( (r->X) == (s->X) ) {
                if (plus) z= (r->coeff) + (s->coeff);
                else z= (r->coeff) - (s->coeff);
                if ( z==zero){
                    r=r->next;
                    s=s->next;
                }
                else {
                    u=new term(std::move(z),r->X);
                    t->next=u;
                    t=t->next;
                    r=r->next;
                    s=s->next;
                }
            }
            else {
                if (plus) u=new term(s->coeff,s->X);
                else u=new term(-s->coeff,s->X);
                t->next=u;
                t=t->next;
                s=s->next;
            }
       }
 }
 if (s==NULL) {s=r;plus=1;}
 while (s!=NULL) {
       if (plus) u=new term(s->coeff,s->X);
       else u=new term(-s->coeff,s->X);
       t->next=u;
       t=t->next;
       s=s->next;
 }
 return q;
}


template <class T>
s_polynom<T> s_polynom<T>::operator -() const
{
 term* p,*q,*r;
 s_polynom<T> f;

 p=head;
 if (p==NULL) return f;
 q = new term(-(p->coeff),p->X);
 f.head=q;
 while ((p->next)!=NULL) {
       p=p->next;
       r = new term(-(p->coeff),p->X);
       q->next = r;
       q=q->next;
 }
 return f;
}


template <class T>
void s_polynom<T>::addto(const s_polynom<T>& t,int plus)
{
 term *p,*q,*r;
 s_polynom<T> s;
 int stop=0,fall;
 T z;

 if (head==NULL) {
    if (plus) *this=t;
    else *this=-t;
    return;
 }
 if (t.head==NULL) return;
 s.head=t.head;
 q=s.head;
 if ( (head->X) < (q->X) ){
    if (plus) r=new term(q->coeff,q->X,head);
    else r=new term(-q->coeff,q->X,head);
    head=r;
    q=q->next;
 }
 else if ( (head->X) == (q->X) ) {
    if (plus) z=(head->coeff) + (q->coeff);
    else z=(head->coeff) - (q->coeff);
    if ( z == zero ){
       r=head;
       head=head->next;
       delete r;
       s.head=s.head->next;
       addto(s,plus);
       //*this+=s;
       s.head=NULL;
       return;
    }
    else {
       head->coeff=std::move(z);
       q=q->next;
    }
 }
 p=head;
 while (p!=NULL && q!=NULL && !stop)
       if (p->next==NULL) stop=1;      // q must be appended to p 
       else {
      fall= !( (q->X) < (p->next->X) );
      if ( (q->X) > (p->next->X) ) fall=2;
      switch (fall) {
          case 0 : p=p->next; break;       // q->X < p->next->X

          case 1 :
              if (plus) z = (p->next->coeff) + (q->coeff);
              else z = (p->next->coeff) - (q->coeff);
              if ( z == zero ) {  // == !
                  r=p->next;
                  p->next=r->next;
                  delete r;
              }
              else {
                    p->next->coeff=std::move(z);
                    p=p->next;
              }
              q=q->next;
              break;

          case 2 :
               if (plus) r = new term(q->coeff,q->X,p->next);
               else r = new term(-q->coeff,q->X,p->next);    // < !
               p->next=r;
               p=p->next;
               q=q->next;
               break;
      }
       }
 if (stop)
    while (q!=NULL) {
      if (plus) r=new term(q->coeff,q->X);
      else r=new term(-q->coeff,q->X);
      p->next=r;
      p=p->next;
      q=q->next;
    }
 s.head=NULL;
 return;
}

//

template <class T>
int s_polynom<T>::length() const
{
    term *p;
    int i;

    for (i=0,p=head;p!=NULL;p=p->next,i++);
    return i;
}

template <class T>
int s_polynom<T>::totdeg() const
{
    int d=0;

    if (head==NULL) return -1;
    for (term *p=head;p!=NULL;p=p->next) d=val::Max(d,p->X.totdeg());

    return d;
}


template <class T>
const s_polynom<T>& s_polynom<T>::setdegree(const vector<integer> &w)
{
    integer degree(0);
    int i,n=s_expo::dim;
    term *p;
    for (p=head;p!=NULL;p=p->next) {
        //degree=integer(0);
        for (i=0;i<n;i++) degree+= integer(p->X.coeff[i])*w(i);
        p->deg=std::move(degree); // dann is degree=0;
    }
    return *this;
}


template <class T>
const s_polynom<T>& s_polynom<T>::setzerodegree()
{
    term *p;
    integer zero(0);

    for (p=head;p!=NULL;p=p->next) p->deg=zero;

    return *this;
}


template <class T>
const s_polynom<T>& s_polynom<T>::setLTdegree(const vector<integer> &w)
{
    if (head==NULL) return *this;
    int i,n=s_expo::dim;
    integer degree(0);
    for (i=0;i<n;i++) degree+= integer(head->X.coeff[i])*w(i);
    head->deg=std::move(degree);
    return *this;
}


template <class T>
void s_polynom<T>::makeinitforms(s_polynom<T> &f,s_polynom<T> &h,const val::vector<integer> &w) const
{
 integer d;
 term *p,*rf,*qf,*rh,*qh;

 if (head==NULL) return;

 f.head=new term(head->coeff,head->X,head->deg);
 h.head=new term(head->coeff,head->X);

 qf=f.head;
 qh=h.head;
 d=degw(head->X,w);
 for (p=head->next;p!=NULL;p=p->next) {
     if (degw(p->X,w)==d) {
         rf=new term(p->coeff,p->X,p->deg);
         qf->next=rf;
         qf=qf->next;
         rh=new term(p->coeff,p->X);
         qh->next=rh;
         qh=qh->next;
     }
 }
}


//

template <>
DLL_PUBLIC void s_polynom<modq>::minusmalmonom(const s_polynom<modq>&,const modq&,const s_expo&);


template <class T>
void s_polynom<T>::minusmalmonom(const s_polynom<T> &g,const T &a,const s_expo &Y)
{
 int c;
 T b;
 s_expo Z;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Setze head:
 for (q=g.head;q!=NULL;q=q->next) {
     if (head==NULL) break;
     Z=q->X*Y;
     c=expocompare(Z,head->X);
     if (c==-1) break;
     else if (c==1) {                       // head->X < Z
         b=q->coeff*a;
         changesign(b);
         r=new term(std::move(b),Z,head);
         head=r;
         q=q->next;
         break;
     }
     else {
         b=q->coeff*a;
         changesign(b);
         head->coeff+=b; //h
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
         Z=q->X*Y;
         c=expocompare(Z,p->next->X);
         if (c==-1) p=p->next;
         else if (c==1) {                    // p->next->X < Z
             b=q->coeff*a;
             changesign(b);
             r=new term(std::move(b),Z,p->next);
             p->next=r;
             p=p->next;
             q=q->next;
         }
         else {
             b=q->coeff*a;
             changesign(b);
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
 // Append rest :
 if (head==NULL) {
     for(;q!=NULL;q=q->next) {
         b=q->coeff*a;
         changesign(b);
         if (b!=zero) {
             head = new term(std::move(b),q->X*Y);
             p=head;
             q=q->next;
             break;
         }
     }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff*a;
     changesign(b);
     if (b!=zero) {
         r = new term(std::move(b),q->X*Y);
         p->next=r;
         p=p->next;
     }
 }

 return;
}


template <>
DLL_PUBLIC void s_polynom<modq>::minusmalmonom(const s_polynom<modq> &g,const modq &a,const s_expo &Y,const integer& d);

template <class T>
void s_polynom<T>::minusmalmonom(const s_polynom<T> &g,const T &a,const s_expo &Y,const integer& d)
{
 int c;
 integer d1;
 T b;
 s_expo Z;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head
 for (q=g.head;q!=NULL;q=q->next) {
     if (head==NULL) break;
     Z=q->X*Y;
     d1=q->deg +d;
     c=wexpocompare(Z,d1,head->X,head->deg);
     if (c==-1) break;
     else if (c==1) {           // head->X < Z
         b=q->coeff*a;
         changesign(b);
         r=new term(std::move(b),Z,std::move(d1),head);
         this->head=r;
         q=q->next;
         break;
     }
     else {
         b=q->coeff*a;
         changesign(b);
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
     // Invariant p->X > q->X :
     while (q!=NULL && p->next!=NULL) {
         Z=q->X*Y;
         d1=q->deg+d;
         c=wexpocompare(Z,d1,p->next->X,p->next->deg);
         if (c==-1) p=p->next;
         else if (c==1) {                   // p->next->X < Z
             b=q->coeff*a;
             changesign(b);
             r=new term(std::move(b),Z,std::move(d1),p->next);
             p->next=r;
             p=p->next;
             q=q->next;
         }
         else {
             b=q->coeff*a;
             changesign(b);
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
         b=q->coeff*a;
         changesign(b);
         if (b!=zero) {
             head = new term(std::move(b),q->X*Y,q->deg+d);
             p=head;
             q=q->next;
             break;
         }
     }
 }
 for (;q!=NULL;q=q->next) {
     b=q->coeff*a;
     changesign(b);
     if (b!=zero) {
         r = new term(std::move(b),q->X*Y,q->deg+d);
         p->next=r;
         p=p->next;
     }
 }

 return;
}


//

template <class T>
const s_polynom<T>& s_polynom<T>::operator *=(const T& a)
{
    if (a==zero) del();
    if (head==NULL) return *this;
    term *p,*r;

    head->coeff*=a;
    while (head->coeff==zero) {
        r=head;
        head=head->next;
        delete r;
        if (head==NULL) break;
        head->coeff*=a;
    }
    if (head==NULL) return *this;
    p=head;
    while (p->next!=NULL) {
        p->next->coeff*=a;
        if (p->next->coeff==zero) {
            r=p->next;
            p->next=r->next;
            delete r;
        }
        else p=p->next;
    }
    return *this;
}

template <class T>
const s_polynom<T>& s_polynom<T>::operator *=(const s_expo& X)
{
    if (head==NULL) return *this;

    for (term *p=head;p!=NULL;p=p->next) p->X*=X;
    return *this;
}



template <>
DLL_PUBLIC void s_polynom<integer>::edivby(const integer&);


template <class T>
void s_polynom<T>::edivby(const T& a)
{
    for (term *p=head;p!=NULL;p=p->next) p->coeff/=a;
}


template <class T>
void s_polynom<T>::malmonom(const T& a,const s_expo& Y)
{
    if (head==NULL) return;
    if (a==zero) {del();return;}

    term  *p=head;

    for(;p!=NULL;p=p->next) {
        p->coeff*=a;
        p->X*=Y;
    }
    return;
}



template <>
DLL_PUBLIC const s_polynom<int>& s_polynom<int>::normalize();


template <>
DLL_PUBLIC const s_polynom<integer>& s_polynom<integer>::normalize();


template <class T>
const s_polynom<T>& s_polynom<T>::normalize()
{
    if (head==NULL) return *this;
    term *p;
    T div(head->coeff);
    if (div==T(1)) return *this;
    for (p=head;p!=NULL;p=p->next) p->coeff/=div;
    return *this;
}


template <>
DLL_PUBLIC int s_polynom<int>::content() const;


template <>
DLL_PUBLIC integer s_polynom<integer>::content() const;

template <>
DLL_PUBLIC rational s_polynom<rational>::content() const;


template <class T>
T s_polynom<T>::content() const
{
    return T(1);
}


//

template <class T>
std::ostream& operator <<(std::ostream& os,const s_polynom<T>& f)
{
    struct s_polynom<T>::term *p;
    for (p=f.head;p!=NULL;p=p->next) {
        os<<p->coeff<<std::endl;
        os<<p->X;
        os<<std::endl;
    }
    os<<0<<std::endl;
    return os;
}


template <class T>
std::istream& operator >>(std::istream& is,s_polynom<T>& f)
{
    if (!f.iszero()) f.del();
    T wert,zero=s_polynom<T>::zero;
    s_expo X;

    do {
        wert=zero;
        is>>wert;
        if (wert!=zero) {
            is>>X;
            f.insert(std::move(wert),X);
        }
        else break;
    }
    while (1);
    return is;
}



}  // end namespace val


#endif // S_POLYNOM_H_INCLUDED
