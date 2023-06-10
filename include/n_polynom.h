#ifndef N_POLYNOM_H_INCLUDED
#define N_POLYNOM_H_INCLUDED

#include <matrix.h>
#include <atomic>


namespace val
{

template <class T> class n_polynomIterator;
template <class T> class n_polynom;
template <class T> class Glist;

// class for multidimensional exponents (terms) with variable dimension.
// Order relation defined by a quadratic int - matrix, or lexicographic (ordtype = -1), or DegRevLex ( ordtype = -2 ),
// or DegLex ( ordtype = 0 ).
class DLL_PUBLIC n_expo
{
private:
    int *coeff=nullptr;
    int dim=0;
    static int ordtype;
    static matrix<int> ordmatrix;        // order-matrix
public:
    //
    n_expo() = default;
    explicit n_expo(int d) { if (d<=0) {dim=0;coeff=nullptr;return;} coeff = new int[dim=d];}
    n_expo(int value,int d);
    n_expo(const n_expo& X);
    n_expo(n_expo&& X);
    //
    ~n_expo() {if (coeff!=nullptr) delete[] coeff;}
    //
    const n_expo& operator =(const n_expo& X);
    const n_expo& operator =(n_expo&& X);
    //
    int isempty() const {return (coeff==nullptr);}
    int iszero() const;
    int dimension() const {if (coeff==nullptr) return 0; else return dim;}
    //
    int operator ==(const n_expo& X) const;
    int operator !=(const n_expo& X) const {return !(*this==X);}
    //
    int operator[](int i) const;
    int& operator[](int i);
    int operator()(int i) const {return operator[](i);}
    int& operator()(int i) {return operator[](i);}
    //
    static int getordtype() {return ordtype;}
    static const matrix<int>& getordmatrix() {return n_expo::ordmatrix;}
    static int setordtype(int ord) {ordtype=ord;return ordtype;}
    static const matrix<int>& setordmatrix(const matrix<int>& M) {ordmatrix=M; return ordmatrix;}
    static const matrix<int>& setordmatrix(matrix<int> &&M) {ordmatrix=std::move(M); return ordmatrix;}
    //
    int operator <(const n_expo& y) const;
    int operator >(const n_expo& y) const {return (y < *this);}

    template <class T> friend class n_polynom;
};


DLL_PUBLIC n_expo operator *(const n_expo&,const n_expo&);
DLL_PUBLIC const n_expo& operator*=(n_expo&,const n_expo&);
DLL_PUBLIC n_expo operator /(const n_expo&,const n_expo&);
DLL_PUBLIC const n_expo& operator /=(n_expo&,const n_expo&);
DLL_PUBLIC int operator |(const n_expo&,const n_expo&);
DLL_PUBLIC n_expo lcm(const n_expo&,const n_expo&);
DLL_PUBLIC n_expo char_to_nexpo(const char* c,int l);
DLL_PUBLIC std::istream& operator >>(std::istream&,n_expo&);
DLL_PUBLIC std::ostream& operator <<(std::ostream& os,const n_expo& X);



//
template <class T> std::ostream& operator <<(std::ostream&,const n_polynom<T>&);
template <class T> std::istream& operator >>(std::istream&,n_polynom<T>&);
// Only top-reductions
template <class T> void divrem(const n_polynom<T>& f,const n_polynom<T> &g,n_polynom<T> &q,n_polynom<T> &r);
template <class T> int totaldegree(const n_polynom<T>& f);
//


// class for multivariate polynomials with variable exponents (terms) of the class n_expo
// Data structure is a simple linked list, ordered by the order of the terms (n_expo)
template <class T>
class n_polynom
{
private:
    struct term{
        T coeff;
        n_expo X;
        term *next;
        //
        term() = delete;
        term(const T &wert,term* nex=NULL) : coeff(wert), next(nex) {mnumber++;}
        term(T&& wert,term* nex=NULL) : coeff(std::move(wert)), next(nex) {mnumber++;}
        term(const T& wert,const n_expo& Z,term* nex=NULL) : coeff(wert), X(Z) , next(nex) {mnumber++;}
        term(T&& wert,const n_expo& Z,term* nex=NULL) : coeff(std::move(wert)), X(Z) , next(nex) {mnumber++;}
        term(T&& wert,n_expo&& Z,term* nex=NULL) : coeff(std::move(wert)), X(std::move(Z)) , next(nex) {mnumber++;}
        ~term() {mnumber--;}
    };

    term* head = nullptr;
    //
    static std::atomic<int> mnumber;
    static const T zero;
    static const n_expo neutral;
    static int staticexpodim;
    //
    static int expocompare(const n_expo&,const n_expo&);
    static void changesign(T&);
    static void merge(term* &basis1,term* &basis2);
    static void sort(term* &basis,int l);
    //
    n_polynom<T> add(const n_polynom<T>& p,int plus=1) const;
    void addto(const n_polynom<T>& t,int plus=1);
public:
    // Constructors:
    n_polynom() = default;
    n_polynom(const T& wert) {if (wert==zero) head=NULL; else head= new term(wert); }
    n_polynom(T&& wert) {if (wert==zero) head=NULL; else head= new term(std::move(wert)); }
    n_polynom(const T& wert,const n_expo& Z) {if (wert==zero) head=NULL; else head= new term(wert,Z); }
    n_polynom(T&& wert,const n_expo& Z) {if (wert==zero) head=NULL; else head=new term(std::move(wert),Z); }
    n_polynom(T&& wert,n_expo&& Z) {if (wert==zero) head=NULL; else head= new term(std::move(wert),std::move(Z)); }
    n_polynom(const n_polynom<T>&);
    n_polynom(n_polynom<T>&& f) {head=f.head;f.head=NULL;}
    // Destructors
    void del();
    ~n_polynom() {del();}
    //
    int iszero() const {return (head==NULL);}
    //
    const n_polynom<T>& operator =(const n_polynom<T>&);
    const n_polynom<T>& operator =(n_polynom<T>&& f) { if (head!=f.head) {n_polynom<T> g;g.head=head;head=f.head;f.head=NULL;}
                                                       return *this;}
    const n_polynom<T>& operator =(const T& wert) {n_polynom<T> f,g(wert);f.head=head;head=g.head;g.head=NULL;return *this;}
    const n_polynom<T>& operator =(T&& wert) {n_polynom<T> f,g(std::move(wert));f.head=head;head=g.head;g.head=NULL;return *this;}
    //
    int operator ==(const n_polynom<T>&) const;
    int operator ==(const T&) const;
    int operator !=(const n_polynom<T>& f) const {return !(*this==f);}
    int operator !=(const T& a) const {return !(*this==a);}
    int operator <(const n_polynom<T>&) const;
    //
    const n_polynom<T>& insert(const T&,const n_expo&);
    const n_polynom<T>& insert(T&&,const n_expo&);
    const n_polynom<T>& insert(T&&,n_expo&&);
    const n_polynom<T>& reord() {sort(head,length());return *this;}
    //
    int length() const;
    const T& LC() const {if (head==NULL) return zero; else return head->coeff;}         // leading coefficient
    const T& coefficient(const n_expo& X) const;                                        // coefficient from monomial with term X
    const n_expo& LT() const {if (head==NULL) return neutral; else return head->X;}
    int getdim() const;                                                                 // Max dimension of n_expo in polynomial
    //
    n_polynom<T> operator +(const n_polynom<T>& p) const {return add(p,1);}
    n_polynom<T> operator -(const n_polynom<T>& p) const {return add(p,0);}
    n_polynom<T> operator-() const;
    const n_polynom<T>& operator +=(const n_polynom<T>& f) {addto(f,1);return *this;}
    const n_polynom<T>& operator -=(const n_polynom<T>& f) {addto(f,0);return *this;}
    DLL_PUBLIC void minusmalmonom(const n_polynom<T>& g,const T& a,const n_expo& Y);     // *this -= a*Y*g
    //
    n_polynom<T> operator *(const n_polynom<T> &f) const;
    const n_polynom<T>& operator *=(const n_polynom<T>& f) {*this= *this*f; return *this;}
    const n_polynom<T>& operator*=(const T&);
    const n_polynom<T>& operator/=(const n_expo& X);
    n_polynom<T> operator /(const n_polynom<T> &g) const {n_polynom<T>q,r; divrem(*this,g,q,r); return q;}
    const n_polynom<T>& operator /=(const n_polynom<T> & g) {*this=*this/g; return *this;}
    DLL_PUBLIC void edivby(const T&);                   // for T = integer : Exact division by constant
    void malmonom(const T& a,const n_expo& Y);          // *this *= a*Y
    DLL_PUBLIC const n_polynom<T>& normalize();         // For integer types: Make primitive. For fields: divide by LC().
    DLL_PUBLIC T content() const;
    friend void divrem<T>(const n_polynom &f,const n_polynom &g,n_polynom &q,n_polynom &r);
    //
    static std::atomic<int> nreduction;
    DLL_PUBLIC int reduction(const Glist<n_polynom<T> >& P,int top=1,int interred=0);   // Reductions modulo <P> for Buchberger-Algorithm
    //
    static int staticexpoutput;
    //
    static void setstaticexpodim(int n);
    static int getstaticexpodim() {return staticexpodim;}
    static int getmnumber() {return mnumber;}
    //
    friend std::ostream& operator <<<T>(std::ostream&,const n_polynom&);
    friend std::istream& operator >><T>(std::istream&,n_polynom&);
    //
    n_polynomIterator<T> begin() const {return n_polynomIterator<T>(*this);}
    n_polynomIterator<T> end() const {return n_polynomIterator<T>();}
    //
    friend class n_polynomIterator<T>;
};



template <class T>
class n_polynomIterator
{
private:
    struct n_polynom<T>::term *actual;
public:
    n_polynomIterator() {actual=NULL;}
    n_polynomIterator(const n_polynom<T>& f) {actual=f.head;}
    const n_polynomIterator& operator = (const n_polynom<T>& f) {actual=f.head;return *this;}
    operator int() const {return (actual!=NULL);}
    void operator++(int) {if (actual!=NULL) actual=actual->next;}
    void operator++() {if (actual!=NULL) actual=actual->next;}
    const n_polynomIterator<T>& operator* () const {return *this;}
    const T& actualcoef() const {if (actual==NULL) return n_polynom<T>::zero;else return actual->coeff;}
    const n_expo& actualterm() const {if (actual==NULL) return n_polynom<T>::neutral;else return actual->X;}
    int actualvalid() const {return (actual!=NULL);}
    int moveactualtolast();   // point to last monomial
};




// Definitions n_polynom:

template <class T>
std::atomic<int> n_polynom<T>::mnumber(0);


template <class T>
const T n_polynom<T>::zero(val::zero_element<T>());

template <class T>
const n_expo n_polynom<T>::neutral;


// Constructors:  ----------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
n_polynom<T>::n_polynom(const n_polynom<T>& f)
{
    head=NULL;
    if (f.head==NULL) return;

    term *q,*t,*r;

    q=f.head;
    t=new term(q->coeff,q->X);
    head=t;
    while (q->next != NULL) {
        q=q->next;
        r=new term(q->coeff,q->X);
        t->next=r;
        t=t->next;
    }
}

// Destructors: ----------------------------------------------------------------------------------------------------------------------------------------------
template <class T>
void n_polynom<T>::del()
{
    term *p;

    while (head!=NULL) {
        p=head;
        head=head->next;
        delete p;
        //mnumber--;
    }
}
// ------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
const n_polynom<T>& n_polynom<T>:: operator =(const n_polynom<T>& s)
{

 if (head == s.head) return *this;

 n_polynom<T> help(s);

 val::swap(head,help.head);

 return *this;
}

// --------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
int n_polynom<T>:: operator==(const n_polynom<T>& f) const
{
 if (head==f.head) return 1;
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
int n_polynom<T>:: operator ==(const T& i) const
{
 if (head==NULL && i==zero)
    return 1;
 else if (head->next==NULL && head->X==neutral && head->coeff==i) return 1;
 else return 0;
}



template <class T>
int n_polynom<T>::operator <(const n_polynom<T>& g) const
{
 term *p,*q;

 for (q=g.head,p=head;q!=NULL;q=q->next,p=p->next) {
       if (p==NULL) return 1;
       if (p->X<q->X) return 1;
       else if (p->X==p->X) continue;
       else return 0;
 }
 return 0;
}

// ------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
const n_polynom<T>& n_polynom<T>::insert(const T& wert,const n_expo& g)
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
       if (head->coeff!=zero)      // Koeffz == 0
           return *this;                     // => loesche Element aus den Polynom
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
            if (q->next->coeff!=zero)     		// coeff == 0
                return *this;                   // => delete term from polynomial
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
const n_polynom<T>& n_polynom<T>::insert(T&& wert,const n_expo& g)
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
       if (head->coeff!=zero)      		// coeff == 0
           return *this;                // => delete term from polynomial
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
            if (q->next->coeff!=zero)     		// coeff == 0
                return *this;                   // => delete term from polynomial
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
const n_polynom<T>& n_polynom<T>::insert(T&& wert,n_expo&& g)
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
       if (head->coeff!=zero)      			// coeff == 0
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

// -----------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
int n_polynom<T>::length() const
{
    int l=0;
    for (term *p=head;p!=NULL;p=p->next) l++;
    return l;
}

template <class T>
int n_polynom<T>::getdim() const
{
    int n=0;
    for (term *p=head;p!=nullptr;p=p->next) n = val::Max(n,p->X.dimension());
    return n;
}
//

template <class T>
const T& n_polynom<T>::coefficient(const n_expo& X) const
{
    for (term *p=head;p!=nullptr;p=p->next) {
        if (X>p->X) return zero;
        else if (X==p->X) return p->coeff;
    }
    return zero;
}

// --------------------------------------------------------------------------------------------------------------------
template <class T>
n_polynom<T> n_polynom<T>::add(const n_polynom<T>& p,int plus) const
{
 n_polynom<T> q;                                     // result is ordered n_polynom
 term *r,*s,*t,*u;
 T z;


 r=head;
 s=p.head;
 do {                                    // Set first term in n_polynom q
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
n_polynom<T> n_polynom<T>::operator -() const
{
 term* p,*q,*r;
 n_polynom<T> f;

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
void n_polynom<T>::addto(const n_polynom<T>& t,int plus)
{
 term *p,*q,*r;
 n_polynom<T> s;
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
       if (p->next==NULL) stop=1;      // q is appended  in p
       else {
	  fall= !( (q->X) < (p->next->X) );
	  if ( (q->X) > (p->next->X) ) fall=2;
	  switch (fall) {
	      case 0 : p=p->next; break;     // q->X < p->next->X

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


template <class T>
void n_polynom<T>::minusmalmonom(const n_polynom<T> &g,const T &a,const n_expo &Y)
{
 int c;
 T b;
 n_expo Z;
 term *p,*q,*r;

 if (g.head==NULL || a==zero) return;

 // 1. Set head:
 for (q=g.head;q!=NULL;q=q->next) {
	 if (head==NULL) break;
	 Z=q->X*Y;
	 c=expocompare(Z,head->X);
	 if (c==-1)  break;             // case: Z < head->X
	 else if (c==1) {               // case: head->X < Z
		 b=q->coeff*a;
		 changesign(b);
		 r=new term(std::move(b),std::move(Z),head);
		 head=r;
		 q=q->next;
		 break;
	 }
	 else {
         b=q->coeff*a;
         changesign(b);
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
 if (p!=NULL) {                               // Invariant: p->X > q->X
	 while (q!=NULL && p->next!=NULL) {
		 Z=q->X*Y;
		 c=expocompare(Z,p->next->X);
		 if (c==-1) p=p->next;                // case: Z < p->next->X
		 else if (c==1) {                     // case: p->next->X < Z
			 b=q->coeff*a;
			 changesign(b);
			 r=new term(std::move(b),std::move(Z),p->next);
			 p->next=r;
			 p=p->next;
			 q=q->next;
		 }
		 else {
             b=q->coeff*a;
             changesign(b);
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

 // Append eventually rest.
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

//  -------------------------------------------------------------------------------------------------------------------------------------------


// Temp.:
template <class T>
n_polynom<T> n_polynom<T>::operator *(const n_polynom<T> & f) const
{
    n_polynom<T> g;
    term *p;

    if (head==NULL || f.head==NULL) return g;
    for (p=head;p!=NULL;p=p->next) {
        g.minusmalmonom(f,-(p->coeff),p->X);     // event. define plusmalmonom
    }
    return g;
}


template <class T>
const n_polynom<T>& n_polynom<T>::operator *=(const T& a)
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
const n_polynom<T>& n_polynom<T>::operator/=(const n_expo& X)
{
    if (head==nullptr) return *this;
    for (term *p=head;p!=nullptr;p=p->next) p->X/=X;
    return *this;
}


template <class T>
void n_polynom<T>::edivby(const T& a)
{
    for (term *p=head;p!=NULL;p=p->next) p->coeff/=a;
}


template <class T>
void n_polynom<T>::malmonom(const T& a,const n_expo& Y)
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


template <class T>
const n_polynom<T>& n_polynom<T>::normalize()
{
    if (head==NULL) return *this;
    term *p;
    T div(head->coeff);
    if (div==T(1)) return *this;
    for (p=head;p!=NULL;p=p->next) p->coeff/=div;
    return *this;
}


template <class T>
T n_polynom<T>::content() const
{
    return T(1);
}

// If T is field.
template <class T>
void divrem(const n_polynom<T> &f,const n_polynom<T> &g,n_polynom<T> &q,n_polynom<T> & r)
{
    r=f;q.del();
    if (g.iszero()) {
        Error::error("n_polynom<T>:: division by zero!");
    }
    if (f.iszero()) return;

    n_expo X;
    const n_expo &Y = g.LT();
    T b;
    const T &c=g.LC();
    n_polynom<T> hg;
    struct n_polynom<T>::term *p;

    do {
        if (Y | r.head->X ) {
            hg.head=g.head;
            X=r.head->X/Y;
            b=r.head->coeff/c;
            p=r.head;
            r.head=r.head->next;
            delete p;
            hg.head=hg.head->next;
            r.minusmalmonom(hg,b,X);
            q.insert(std::move(b),std::move(X));
            hg.head=nullptr;
            if (r.head==nullptr) return;
        }
        else break;
    }
    while (1);

}

// static public-member: -------------------------------------------------------------------------------------------------------------------

template <class T>
std::atomic<int> n_polynom<T>::nreduction(0);

template <class T>
int n_polynom<T>::staticexpodim=1;


template <class T>
int n_polynom<T>::staticexpoutput=0;

// static - fctn:

template <class T>
void n_polynom<T>::setstaticexpodim(int n)
{
    if (n<=0) n=1;
    if (n>2000) n=2000;
    staticexpodim=n;
}


template <class T>
int n_polynom<T>::expocompare (const n_expo& X,const n_expo& Y)
{
    if (X==Y) return 0;
    else if (X<Y) return -1;
    else return 1;
}

template <class T>
void n_polynom<T>::changesign(T& a)
{
    a=-a;
}


// Merge two ordered lists. basis1 != NULL, basis2 != NULL and basis1->X != basis2->X
template <class T>
void n_polynom<T>::merge(term* &basis1,term* &basis2)
{
 term *p1,*p2,*r;

 p1=basis1;p2=basis2;
 if (p2->X>basis1->X) {
	 while ((p2->next!=NULL) && (p2->next->X>basis1->X)) p2=p2->next;
	 basis1=basis2;
	 basis2=p2->next;
	 p2->next=p1;
 }
 // jetzt p2<p1
 if (basis2==NULL) return;
 p2=basis2;
 while (p1->next!=NULL) {
	 if ((basis2!=NULL) && (basis2->X>p1->next->X)) {
		 while ((p2->next!=NULL) && (p2->next->X>p1->next->X)) p2=p2->next;
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
void n_polynom<T>::sort(term* &basis,int l)
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
 sort(p1,l1);sort(p2,l2);
 merge(p1,p2);
 basis=p1;
 p1=p2=NULL;
 return;
}

// --------------------------------------------------------------------------------------------------------------------------------------------

//

template <class T>
std::ostream& operator <<(std::ostream& os,const n_polynom<T>& f)
{
    struct n_polynom<T>::term *p;
    int i;
    for (p=f.head;p!=NULL;p=p->next) {
        os<<p->coeff<<std::endl;
        os<<p->X;
        if (n_polynom<T>::staticexpoutput) {
            if (p->X.isempty()) i=1;
            else i=p->X.dimension();
            for (;i<n_polynom<T>::staticexpodim;i++) os<<"0  ";
        }
        os<<std::endl;
    }
    os<<0<<std::endl;
    return os;
}


template <class T>
std::istream& operator >>(std::istream& is,n_polynom<T>& f)
{
    if (!f.iszero()) f.del();
    T wert;
    n_expo X;

    do {
        is>>wert;
        if (wert!=n_polynom<T>::zero) {
            is>>X;
            f.insert(std::move(wert),X);
        }
        else break;
    }
    while (1);
    return is;
}

// ------------------------------------------------------------------------------------------------------------------------------------------

//
template <class T>
void set_unity_element(n_polynom<T> &one)
{
    T one_T;
    set_unity_element(one_T);
    one = n_polynom<T>(one_T);
}

template <class T>
int totaldegree(const n_polynom<T>& f)
{
    if (f.iszero()) return -1;
    n_polynomIterator<T> it=f;
    int tdeg=0,itdeg=0,i,dim;

    dim=it.actualterm().dimension();
    for (i=0;i<dim;++i) itdeg+=it.actualterm()[i];
    it++;
    tdeg=itdeg;
    for (;it;it++) {
        itdeg=0;
        dim=it.actualterm().dimension();
        for (i=0;i<dim;++i) itdeg+=it.actualterm()[i];
        tdeg=val::Max(tdeg,itdeg);
    }

    return tdeg;
}


} // end namesapace val


#endif // N_POLYNOM_H_INCLUDED
