
// class for univariate polynomials. Data structure is a simple linked list.

#ifndef POL_H
#define POL_H


#include <iostream>
//#include <cstdlib>
#include <val_basics.h>
#include <val_utils.h>
#include <error.h>
#include <atomic>
#include <initializer_list>



namespace val
{



template <class T>
struct term
{
      T coef;
      int grad;
      term<T>* next;
      term(const T &a,int g,term<T>* nex=NULL);
      term(T&& a,int g,term<T>* nex=NULL);
      ~term() {mnumber--;}
      static std::atomic<int> mnumber;   //evtl wieder entfernen.
      //static int mnumber;
};


// ------- Friend-functions of class pol<T>----------------------------
template <class T> class pol;
template <class T> class polIterator;
template <class T> int deg(const pol<T>&);
template <class T> const T& leader(const pol<T>&);
template <class T> int length(const pol<T>&);
template <class T> pol<T> operator *(const T&,const pol<T>&);
template <class T> T* horner(const pol<T>&,const pol<T>&);
template <class T> void divrest(const pol<T>&,const pol<T>&,pol<T>&,pol<T>&);
template <class T> void divrem(const pol<T>&,const pol<T>&,pol<T>&,pol<T>&);
template <class T> std::istream& operator >>(std::istream&,pol<T>&);
template <class T> std::ostream& operator <<(std::ostream&,const pol<T>&);
template <class T> void set_unity_element(pol<T> &one);
//--------------------------------------------------------------------------------------


template <class T>
class pol
{
//public:
private:
    term<T> *head=nullptr;
    static const T zero;
public:
// Konstruktoren und Destruktoren
    pol() = default;
    pol(const T&);                 // f = pol(a) => f(x) = a*x^0
    pol(T&& wert);
    pol(const T&,int);               // f = pol(a,n) => f(x) a*x^n
    pol(T&&,int);
    pol(const pol<T>&);
    pol(pol<T> &&);
    pol(std::initializer_list<val::GPair<T,int> > args);
    ~pol();
// 
    void del();                     //  delete all terms in polynomial
    void insert(const T&,int);      //  insert term ordered in polynomial
    void insert(T&&,int);
// Misc
    int iszero() const {return (head==NULL);}
    int degree() const {if (head==NULL) return -1; else return head->grad;}
    friend const T& leader<T>(const pol&);
    const T& leader() const {if (head==NULL) return zero; else return head->coef;}  // Leading coefficient
    const T& LC() const {return leader();}
    friend int length<T>(const pol&);
    int length() const {return val::length(*this);}
    friend int deg<T>(const pol&);
    const T& getlastcoef() const;
    int getlastdeg() const;
    const T& operator[](int) const;
    T eval(const T&) const;                                                      // Evaluation
    T operator() (const T& wert) const {return eval(wert);}
    pol<T> eval (const pol<T>&) const;                                           // Pol.-evaluation ( f(g)) )
    pol<T> operator() (const pol<T>& f) const {return eval(f);}
    template <class Z> Z operator() (const Z& x) const {return eval(x);}
    template <class Z> Z eval(const Z&) const;
    template <class Z> Z derivation(const Z&) const;
    pol<T> derive(int n=1) const;                                                // n-th derivative
    pol<T> integrate() const;                                                    // Integral.
    template <class Z> Z integrate(const Z& a,const Z& b) const;                 // integrates from a to b.
    void trans(int m=1);                                                         // replace p(X) with p(X^m)
    pol<T> multbypower(int n) const;                                                 // *this * x^n
    const pol<T>& getmultbypower(int n);                                         // *this = x^n * *this
    pol<T> divbypower(int n) const;                                              // division with remainder : *this div x^n
    const pol<T>& getdivbypower(int n);                                          // *this = this->divbypower(n);
    T norm(int n) const;                                                         // not defined here
// Assignment
    const pol<T>& operator =(const pol<T>&);
    const pol<T>& operator =(pol<T>&&);
    const pol<T>& operator =(const T&);
// Comparisons
    int operator ==(const pol<T>&) const;
    int operator !=(const pol<T>&) const;
    int operator ==(const T&) const;
    int operator !=(const T&) const;
// Additive
    pol<T> operator +(const pol<T>&) const;
    const pol<T>& operator +=(const pol<T>&);
    pol<T> operator -() const;
    pol<T> operator -(const pol<T>&) const;
    const pol<T>& operator -=(const pol<T>&);
// Multiplicative
    const pol<T>& operator *=(const T&);                 //  Scalar-mult.
    friend pol<T> operator *<T>(const T&,const pol&);
    pol<T> operator *(const pol<T>&) const;
    const pol<T>& operator *=(const pol<T>&);
    friend T* horner<T>(const pol&,const pol&);
    friend void divrest<T>(const pol&,const pol&,pol&,pol&);   // Same as divrem : division with remainder
    pol<T> operator /(const pol<T>&) const;
    pol<T> operator %(const pol<T>&) const;
    const pol<T>& operator /=(const pol<T>&);
    const pol<T>& operator %=(const pol<T>&);
    const pol<T>& normalize();                               // Division by leader coefficient
//
    polIterator<T> begin() const {polIterator<T> it;it=*this;return it;}
    polIterator<T> end() const {return polIterator<T>();}
// Input/Output
    friend std::istream& operator >><T>(std::istream&,pol<T>&);
    friend std::ostream& operator <<<T>(std::ostream&,const pol<T>&);
    friend class polIterator<T>;
};


template <class T>
class polIterator
{
private:
    term<T>* actual;
public:
    polIterator() {actual=NULL;}
    const polIterator<T>& operator =(const pol<T>& f) {actual=f.head;return *this;}
    operator int() const {return (actual!=NULL);}
    void operator++(int) {if (actual!=NULL) actual=actual->next;}
    void operator++() {if (actual!=NULL) actual=actual->next;}
    const T& operator() (void) const {if (actual==NULL) return pol<T>::zero;else return actual->coef;}
    const T& actualcoef() const {if (actual==NULL) return pol<T>::zero;else return actual->coef;}
    const polIterator<T>& operator* () const {return *this;}
    int actualdegree() const {if (actual==NULL) return -1; else return actual->grad;}
    int actualterm() const {return actualdegree();}
    int actualvalid() const {return (actual!=NULL);}
    int moveactualtolast();   // points to last Monom. != 0 if actual != NULL
};




//
template <class T>
pol<T> gcd(const pol<T>& a,const pol<T>& b);


// =======================================================================


// =================== term<T> =========================


template <>
DLL_PUBLIC std::atomic<int> term<double>::mnumber;


template <class T>
std::atomic<int> term<T>::mnumber(0);


// Constructor
template <class T>
term<T>:: term(const T &a,int g,term<T>* nex)
{
 coef=a;
 grad=g;
 next=nex;
 mnumber++;
}


template <class T>
term<T>:: term(T &&a,int g,term<T>* nex) : coef(std::move(a)),grad(g)
{
 next=nex;
 mnumber++;
}

// =====================================================================================================================================




// ==================== pol<T> ==========================================================================================================

// ---------------------------------------------------------------------------------------------------------------------------------------

template <class T>
const T pol<T>::zero(val::zero_element<T>());


template <class T>
pol<T>:: pol (const T& wert)
{
 if (wert!=zero) head=new term<T>(wert,0);
 else head =NULL;
}


template <class T>
pol<T>:: pol (T&& wert)
{
 if (wert!=zero) head=new term<T>(std::move(wert),0);
 else head =NULL;
}

template <class T>
pol<T>:: pol(const T& coef,int n)
{
 if (coef!=zero)
     head = new term<T>(coef,n);
 else head=NULL;
}


template <class T>
pol<T>:: pol(T&& coef,int n)
{
 if (coef!=zero)
     head = new term<T>(std::move(coef),n);
 else head=NULL;
}


template<class T>
pol<T>::pol(const pol<T>& p)
{
 term<T>* q;
 term<T>* t;
 term<T>* r;

 head=NULL;
 if (p.head==NULL) return;
 q=p.head;
 t=new term<T>(q->coef,q->grad);
 head=t;
 while (q->next != NULL) {
       q=q->next;
       r=new term<T>(q->coef,q->grad);
       t->next=r;
       t=t->next;
 }
}

template<class T>
pol<T>:: pol(pol<T>&& p)
{
    head=p.head;
    p.head=NULL;
}


template <class T>
pol<T>::pol(std::initializer_list<val::GPair<T,int> > args)
{
    head=NULL;
    const GPair<T,int> *wert;
    for (wert=args.begin();wert!=args.end();++wert) insert(wert->x,wert->y);
}

template <class T>
pol<T>:: ~pol()
{
 term<T> *p;
 while (head!=NULL) {
       p=head;
       head=head->next;
       delete p;
 }
 head=NULL;
}

// ------------------------------------------------------------------------------------------------------------------------------------------------



// ---------------------- helper functions -----------------------------------------------------------------------------------------------------

template <class T>
void pol<T>::del()
{
 term<T> *p;
 while (head!=NULL) {
       p=head;
       head=head->next;
       delete p;
 }
 head=NULL;
}

template <class T>
void pol<T>::insert(const T& wert,int g)
{
 int stop=0;
 term<T>* p;
 term<T>* q;

 if (wert == zero) return;
 if (head==NULL) {
    p=new term<T>(wert,g);
    head=p; return;}
 else {
    if (g==head->grad){
       head->coef+=wert;
       if (head->coef!=zero)           // coef == 0
           return;                     // => delete term
       else {
           p=head;
           head=head->next;
           delete p;
           return;
       }
    }
    if (g>head->grad) {
       p=new term<T>(wert,g,head);
       head=p;
       return;
    }
    q=head;
    while ((q->next!=NULL) && !stop) {
         if ( q->next->grad > g)
            q=q->next;
         else stop=1;
    }
    if (stop)
       if (g==q->next->grad) {
            q->next->coef+=wert;
            if (q->next->coef!=zero)        // coef == 0
                return;                     // => delete term
            else {
                p=q->next;
                q->next=p->next;
                delete p;
                return;
            }
       }
 }
 p= new term<T>(wert,g,q->next);
 q->next=p;
 return;
}


template <class T>
void pol<T>::insert(T&& wert,int g)
{
 int stop=0;
 term<T>* p;
 term<T>* q;

 if (wert==zero) return;
 if (head==NULL) {
    p=new term<T>(std::move(wert),g);
    head=p; return;}
 else {
    if (g==head->grad){
       head->coef+=wert;
       if (head->coef!=zero)            // coef == 0
            return;                     // => delete term
       else {
            p=head;
            head=head->next;
            delete p;
            return;
       }
    }
    if (g>head->grad) {
       p=new term<T>(std::move(wert),g,head);
       head=p;
       return;
    }
    q=head;
    while ((q->next!=NULL) && !stop) {
      if ( q->next->grad > g)
         q=q->next;
      else stop=1;
    }
    if (stop)
       if (g==q->next->grad) {
          q->next->coef+=wert;
          if (q->next->coef!=zero)              // coef == 0
             return;                        // => delete term
          else {
             p=q->next;
             q->next=p->next;
             delete p;
             return;
          }
       }
 }
 p= new term<T>(std::move(wert),g,q->next);
 q->next=p;
 return;
}


// ---------------------------------------------------------------------


//------------------------- Other functions --------------------------------


template <class T>
int deg(const pol<T>& f)
{
 if (f.head==NULL) return -1;
 else return f.head->grad;
}


template <class T>
const T& leader(const pol<T>& f)
{
 if (f.head==NULL) return pol<T>::zero;
 else return f.head->coef;
}


template <class T>
int length(const pol<T>& f)
{
 int i=0;
 term<T>* p;

 for (p=f.head;p!=NULL;p=p->next) i++;
 return i;
}


template <class T>
const T& pol<T>::getlastcoef() const
{
    if (head==NULL) return zero;
    term<T> *p=head;

    while (p->next!=NULL) p=p->next;
    return p->coef;
}


template <class T>
int pol<T>::getlastdeg() const
{
    if (head==NULL) return -1;
    term<T> *p=head;

    while (p->next!=NULL) p=p->next;
    return p->grad;
}


template <class T>
const T& pol<T>::operator [](int n) const
{
 term<T>* p=head;

 if (n<0) return zero;
 while (p!=NULL) {
     if (n>p->grad) return zero;
     else if (n==p->grad) return p->coef;
     else p=p->next;
 }
 return zero;
}


template <class T>
T pol<T>:: eval(const T& beta) const
{
 int i;
 T wert;
 term<T> *p;

 if (head==NULL) return zero;
 wert=head->coef;
 p=head->next;
 for (i=deg(*this)-1;i>=0;i--) {
     wert*=beta;
     if (p!=NULL)
         if (p->grad==i) {
             wert+=p->coef;
             p=p->next;
         }
 }
 return wert;
}



template <class T>
pol<T> pol<T>:: eval(const pol<T>& beta) const
{
 int i;
 pol<T> wert;
 term<T> *p;

 if (head==NULL) return pol<T>();
 wert.insert(head->coef,0);
 p=head->next;
 for (i=deg(*this)-1;i>=0;i--) {
     wert*=beta;
     if (p!=NULL)
         if (p->grad==i) {
             wert.insert(p->coef,0);
             p=p->next;
         }
 }
 return wert;
}

template <class T>
template <class Z> Z pol<T>::eval(const Z& beta) const
{
 int i;
 Z wert;
 term<T> *p;

 if (head==NULL) return zero_element<Z>();
 wert=Z(head->coef);
 p=head->next;
 for (i=deg(*this)-1;i>=0;i--) {
     wert*=beta;
     if (p!=NULL)
         if (p->grad==i) {
             wert+=Z(p->coef);
             p=p->next;
         }
 }
 return wert;
}


template <class T>
template <class Z> Z pol<T>::derivation(const Z& beta) const
{
 int i;
 Z wert;
 term<T> *p;

 if (head==NULL) return zero_element<Z>();
 wert=Z(T(head->grad)*head->coef);
 p=head->next;
 for (i=head->grad-2;i>=0;i--) {
     wert*=beta;
     if (p!=NULL)
         if (p->grad==i+1) {
             wert+=Z(T(p->grad)*p->coef);
             p=p->next;
         }
 }
 return wert;
}


template <class T>
pol<T> pol<T>::derive(int n) const
{
 int i,j,faktor;
 term<T> *p,*q,*r;
 pol<T> h;
 T wert;

 if (((p=head)==NULL) || (n<0)) return h;
 if (n==0) return *this;

 for (;p!=NULL;p=p->next) {
     faktor=i=p->grad;
     if (i<n) return h;
     for (j=i-1;j>i-n;j--) faktor*=j;
     if ( (wert=(T(faktor)*p->coef))!=zero ) {
        r = new term<T>(std::move(wert),i-n);
        h.head=r;
        p=p->next;
        break;
     }
 }
 q=h.head;
 for (;p!=NULL;p=p->next) {
     faktor=i=p->grad;
     if (i<n) return h;
     for (j=i-1;j>i-n;j--) faktor*=j;
     if ( (wert=(T(faktor)*p->coef))!=zero ) {
        r = new term<T>(std::move(wert),i-n);
        q->next=r;
        q=q->next;
     }
 }
 return h;
}


template <class T>
pol<T> pol<T>::integrate() const
{
    term<T> *p,*r,*q;
    pol<T> h;

    if (head==NULL) return h;
    p=head;
    h.head= new term<T>(p->coef/T(p->grad + 1),p->grad+1);
    p=p->next;
    q=h.head;
    for (;p!=NULL;p=p->next,q=q->next) {
        r= new term<T>(p->coef/T(p->grad + 1),p->grad+1);
        q->next=r;
    }
    return h;
}



template <class T>
template <class Z> Z pol<T>::integrate(const Z& a,const Z& b) const
{
 int i;
 Z wert1,wert2;
 term<T> *p;

 if (head==NULL) return Z(0);
 wert1=Z(head->coef/T(head->grad+1));
 wert2=Z(head->coef/T(head->grad+1));
 p=head->next;
 for (i=head->grad-1;i>=0;i--) {
     wert1*=a;
     wert2*=b;
     if (p!=NULL)
         if (p->grad==i) {
             wert1+=Z(p->coef/T(i+1));
             wert2+=Z(p->coef/T(i+1));
             p=p->next;
         }
 }
 wert1*=a;
 wert2*=b;
 return wert2-wert1;
}



template <class T>
void pol<T>::trans(int m)
{
 term<T>* p;

 if ((m<=1) || (head==NULL)) return;
 for (p=head;p!=NULL;p=p->next) p->grad*=m;
 return;
}


template <class T>
pol<T> pol<T>::multbypower(int n) const
{
 pol<T> f;
 term<T>* q;
 term<T>* t;
 term<T>* r;

 if (head==NULL || n==0) return f;
 else if (n<0) return divbypower(-n);

 q=head;
 t=new term<T>(q->coef,q->grad+n);
 f.head=t;
 while (q->next != NULL) {
       q=q->next;
       r=new term<T>(q->coef,q->grad+n);
       t->next=r;
       t=t->next;
 }
 return f;
}

template <class T>
const pol<T>& pol<T>::getmultbypower(int n)
{
    if (head==NULL || n==0) return *this;
    else if (n<0) return getdivbypower(-n);

    term<T> *p=head;
    for(;p!=NULL;p=p->next) p->grad+=n;
    return *this;

}


template <class T>
pol<T> pol<T>::divbypower(int n) const
{
 pol<T> f;
 term<T>* p;
 term<T>* t;
 term<T>* r;

 if (head==NULL || n==0 || (head->grad-n)<0) return f;
 else if (n<0) return multbypower(-n);

 f.head=t=new term<T>(head->coef,head->grad-n);
 p=head;
 while (p->next!=NULL) {
     p=p->next;
     if (p->grad-n<0) break;
     else {
        r=new term<T>(p->coef,p->grad-n);
        t->next=r;
        t=t->next;
     }
 }
 return f;
}


template <class T>
const pol<T>& pol<T>::getdivbypower(int n)
{
    if (head==NULL || n==0) return *this;
    else if (n<0) return getmultbypower(-n);

    term<T> *p=head;
    if ((p->grad-n)<0) {
        del();
        return *this;
    }
    else p->grad-=n;
    pol<T> f;
    while (p->next!=NULL) {
           if (p->next->grad-n<0) {
                f.head =p->next;
                p->next = NULL;
           }
           else {
                p=p->next;
                p->grad-=n;
           }
    }
    return *this;
}


// ----------------------------------------------------------------------



// -------------------- Assignment ----------------------------

template <class T>
const pol<T>& pol<T>:: operator =(const pol<T>& s)
{
 if (head==s.head) return *this;
 pol<T> help(s);
 val::swap(this->head,help.head);
 return *this;
}


template <class T>
const pol<T>& pol<T>:: operator =(pol<T>&& s)
{
    if (head==s.head) return *this;
    if (head!=NULL) del();
    head=s.head;
    s.head=NULL;
    return *this;
}


template <class T>
const pol<T>& pol<T>:: operator =(const T& wert)
{

 if (head!=NULL) del();
 if (wert==zero) return *this;
 head = new term<T>(wert,0);

 return *this;
}

// -----------------------------------------------------------------------


// ------------------- Comparisons -----------------------------

template <class T>
int pol<T>::operator ==(const pol<T>& g) const
{
 term<T>* p;
 term<T>* q;

 p=head;
 q=g.head;

 while (p!=NULL) {
       if (q==NULL) return 0;
       if ((p->coef!=q->coef) || (p->grad!=q->grad)) return 0;
       p=p->next;
       q=q->next;
 }
 if (q!=NULL) return 0;
 return 1;
}


template <class T>
int pol<T>::operator !=(const pol<T>& g) const
{
 return !(*this==g);
}


template <class T>
int pol<T>::operator ==(const T& wert) const
{
 if (head==NULL) return (wert==zero);
 if (head->grad==0) return ((head->coef)==wert);
 else return 0;
}


template <class T>
int pol<T>::operator !=(const T& wert) const
{
 return !(*this==wert);
}
// ---------------------------------------------------------------------


// -------------------- Additive operators -----------------------------

// Alg.: Merge two sorted lists
template <class T>
pol<T> pol<T>::operator +(const pol<T>& p) const
{
 pol<T> q;                                     
 term<T> *r;                                   
 term<T> *s;                                   
 term<T> *t;                                   
 term<T> *u;
 T z;


 r=head;
 s=p.head;
 do {                                    // Set first term in q
    if (r==NULL && s==NULL) return q;    
    if (r==NULL){                        
       q.insert(s->coef,s->grad);
       s=s->next;
    }
    else if (s==NULL){
       q.insert(r->coef,r->grad);
       r=r->next;
    }
    else { //(r!=NULL && s!=NULL) {
       if ( (r->grad) > (s->grad) ){
      q.insert(r->coef,r->grad);
      r=r->next;
       }
       else if ( (r->grad) == (s->grad) ) {
      q.insert((r->coef)+(s->coef),r->grad);
      r=r->next;
      s=s->next;
       }
       else {
      q.insert(s->coef,s->grad);
      s=s->next;
       }
    }
  }
 while (q.head==NULL);

 t=q.head;
 while (r!=NULL && s!=NULL) {
       if ( (r->grad) > (s->grad) ) {
          u= new term<T>(r->coef,r->grad);
          t->next=u;
          r=r->next;
          t=t->next;
       }
       else {
            if ( (r->grad) == (s->grad) ) {
                if ( (z= (r->coef) + (s->coef) )==zero){
                    r=r->next;
                    s=s->next;
                }
                else {
                    u=new term<T>(std::move(z),r->grad);
                    t->next=u;
                    t=t->next;
                    r=r->next;
                    s=s->next;
                }
            }
            else {
                u=new term<T>(s->coef,s->grad);
                t->next=u;
                t=t->next;
                s=s->next;
            }
       }
 }
 if (s==NULL) s=r;
 while (s!=NULL) {
       u=new term<T>(s->coef,s->grad);
       t->next=u;
       t=t->next;
       s=s->next;
 }
 return q;
}



template <class T>
const pol<T>& pol<T>::operator +=(const pol<T>& t)
{
 term<T> *p;
 term<T> *q;
 term<T> *r;
 pol<T> s;
 int stop=0,fall;

 if (head==NULL) {
    *this=t;
    return *this;
 }
 if (t.head==NULL) return *this;
 s.head=t.head;
 q=s.head;
 if ( (head->grad) < (q->grad) ){
    r=new term<T>(q->coef,q->grad,head);
    head=r;
    q=q->next;
 }
 else if ( (head->grad) == (q->grad) ) {
    if ( ( (head->coef) + (q->coef) ) == zero ){
       r=head;
       head=head->next;
       delete r;
       s.head=s.head->next;
       *this+=s;
       s.head=NULL;
       return *this;
    }
    else {
       head->coef+=q->coef;
       q=q->next;
    }
 }
 p=head;
 while (p!=NULL && q!=NULL && !stop)
       if (p->next==NULL) stop=1;                   // q must be appended to  p 
       else {
      fall= ( (q->grad) >= (p->next->grad) );
      if ( (q->grad) > (p->next->grad) ) fall=2;
      switch (fall) {
          case 0 : p=p->next; break;                    // q->grad < p->next->grad

          case 1 :
              if ( ( (p->next->coef) + (q->coef) ) == zero ) {  // == !
                  r=p->next;
                  p->next=r->next;
                  delete r;
              }
              else {
                  p->next->coef+=q->coef;
                  p=p->next;
              }
              q=q->next;
              break;

          case 2 : r = new term<T>(q->coef,q->grad,p->next);    // < !
               p->next=r;
               p=p->next;
               q=q->next;
               break;
      }
       }
 if (stop)
    while (q!=NULL) {
      r=new term<T>(q->coef,q->grad);
      p->next=r;
      p=p->next;
      q=q->next;
    }
 s.head=NULL;
 return *this;
}


template <class T>
pol<T> pol<T>::operator -() const
{
 term<T>* p;
 term<T>* q;
 term<T>* r;
 pol<T> f;

 p=head;
 if (p==NULL) return f;
 q = new term<T>(-(p->coef),p->grad);
 f.head=q;
 while ((p->next)!=NULL) {
       p=p->next;
       r = new term<T>(-(p->coef),p->grad);
       q->next = r;
       q=q->next;
 }
 return f;
}



template <class T>
pol<T> pol<T>::operator -(const pol<T>& g) const
{
 return *this+(-g);
}



template <class T>
const pol<T>& pol<T>::operator -=(const pol<T>& t)
{
 term<T> *p;
 term<T> *q;
 term<T> *r;
 pol<T> s;
 int stop=0,fall;

 if (head==NULL) {
    *this=-t;
    return *this;
 }
 if (t.head==NULL) return *this;
 s.head=t.head;
 q=s.head;
 if ( (head->grad) < (q->grad) ){
    r=new term<T>(-(q->coef),q->grad,head);
    head=r;
    q=q->next;
 }
 else if ( (head->grad) == (q->grad) ) {
    if ( ( (head->coef) - (q->coef) ) == zero ){
       r=head;
       head=head->next;
       delete r;
       s.head=s.head->next;
       *this-=s;
       s.head=NULL;
       return *this;
    }
    else {
       head->coef-=q->coef;
       q=q->next;
    }
 }
 p=head;
 while (p!=NULL && q!=NULL && !stop)
       if (p->next==NULL) stop=1;                           // q must be appended tp p
       else {
      fall= ( (q->grad) >= (p->next->grad) );
      if ( (q->grad) > (p->next->grad) ) fall=2;
      switch (fall) {
          case 0 : p=p->next; break;                    // q->grad < p->next->grad

          case 1 :
              if ( ( (p->next->coef) - (q->coef) ) == zero ) {  // == !
                  r=p->next;
                  p->next=r->next;
                  delete r;
               }
              else {
                  p->next->coef-=q->coef;
                  p=p->next;
               }
               q=q->next;
               break;

          case 2 : r = new term<T>(-(q->coef),q->grad,p->next);    // < !
               p->next=r;
               p=p->next;
               q=q->next;
               break;
      }
       }
 if (stop)
    while (q!=NULL) {
      r=new term<T>(-(q->coef),q->grad);
      p->next=r;
      p=p->next;
      q=q->next;
    }
 s.head=NULL;
 return *this;
}

// ---------------------------------------------------------------------



// ------------------ Multiplicative operators ------------------------

template <class T>
const pol<T>& pol<T>::operator *=(const T& wert)
{
 term<T>* p;
 term<T>* r;
 

 if (wert == unity_element<T>()) return *this;
 if (wert == zero) {
    del();
    return *this;
 }
 while (head!=NULL) {
    head->coef*=wert;
    if (head->coef==zero) {
       r=head;
       head=head->next;
       delete r;
    }
    else break;
 }
 if (head==NULL) return *this;
 p=head;
 while (p->next!=NULL) {
    p->next->coef*=wert;
    if (p->next->coef==zero) {
       r=p->next;
       p->next=r->next;
       delete r;
    }
    if ((p=p->next)==NULL) return *this;
 }
 return *this;
}


template <class T>
pol<T> operator *(const T& wert,const pol<T>& f)
{
 term<T> *p,*q,*r;
 pol<T> g;

 p=f.head;
 if (p==NULL || wert==pol<T>::zero) return g;
 if (wert==unity_element<T>()) return f;

 while (p!=NULL)
       if ((wert*p->coef)==pol<T>::zero) p=p->next;
       else {
      g.insert(wert*p->coef,p->grad);
      q=g.head;
      p=p->next;
      break;
       }

 for (;p!=NULL;p=p->next)
      if ((wert*p->coef)!=pol<T>::zero) {
      r = new term<T>(wert*p->coef,p->grad);
      q->next=r;
      q=q->next;
       }
 return g;
}



template <class T>
pol<T> pol<T>::operator *(const pol<T>& f) const
{
 term<T> *p,*q,*r,*s;
 pol<T> g,h;
 int d;
 T wert;

 if (f.head==NULL) return g;
 for (p=head;p!=NULL;p=p->next) {
     if (p==NULL) return g;
     d=p->grad;
     wert=p->coef;
     for (q=f.head;q!=NULL;q=q->next) {
         if ( (wert*q->coef) != zero) {
             s= new term<T>(wert*q->coef,d+q->grad);
             h.head=s;
             r=h.head;
             q=q->next;
             break;
         }
     }
     for (;q!=NULL;q=q->next) {
         if ( (wert*q->coef) != zero ) {
             s= new term<T>(wert*q->coef,d+q->grad);
             r->next=s;
             r=r->next;
         }
     }
     g+=h;
     h.del();
 }
 return g;
}


template <class T>
const pol<T>& pol<T>::operator *=(const pol<T>& f)
{
 *this = *this * f;
 return *this;
}


// Multiplicative functions if T is field:


// Hilfsfktn: Horner-Schema
template <class T>
T* horner(const pol<T>& f,const pol<T>& g)
{
 // Here LC(g) == 1
 int n1=(f.head->grad)+1,m=(g.head->grad)+1,i,j,k;
 term<T>* p;
 T* eps;
 T* fv;
 T* gv;

 fv = new T[n1];
 eps = new T[n1];
 gv = new T[m];
 for (i=0;i<n1;i++) fv[i]=eps[i]=pol<T>::zero;
 for (i=0;i<m;i++) gv[i]=pol<T>::zero;

 for (p=f.head;p!=NULL;p=p->next)        // write polynomial f as vector
     fv[p->grad]=p->coef;                

 for (p=g.head;p!=NULL;p=p->next)        // The same for g.
     gv[p->grad]=p->coef;

 for (i=n1-1;i>=0;i--) {
     eps[i]=fv[i];
     k= (m<n1-i)?  m-1 : n1-1-i;
     //k=min(m-1,n1-1-i);
     j= (m-1-i>1)? m-1-i:1;
     for (;j<=k;j++)
     eps[i]-=gv[m-1-j]*eps[i+j];
 }
 delete[] gv;
 delete[] fv;
 return eps;
}



template <class T>
pol<T> pol<T>::operator /(const pol<T>& g) const
{
 int n1,m,i,grad;
 T wert, one = unity_element<T>();
 T* eps;
 term<T> *p,*q;
 pol<T> h;

 if (g.head==NULL) {
     std::string msg="\n ERROR : " + val::gettypename(*this);
     msg+= "::operator / : Division by zero!!\n";
     Error::error(msg.c_str());
 }
 if (head==0) return h;
 n1=head->grad;
 m=g.head->grad;
 if (m>n1) return h;
 wert=g.head->coef;
 eps=horner(*this,((one/wert)*g));
 p=NULL;
 for (grad=0,i=m;i<=n1;i++,grad++)
      if (eps[i]!=zero) {
     q = new term<T>((one/wert)*eps[i],grad,p);
     p=q;
      }
 h.head=p;
 delete[] eps;
 return h;
}


template <class T>
pol<T> pol<T>::operator %(const pol<T>& g) const
{
 int m,i,grad;
 T wert;
 T* eps;
 term<T> *p,*q;
 pol<T> h;

 if (g.head==NULL) {
     std::string msg="\n ERROR : " + val::gettypename(*this);
     msg+= "::operator % : Division by zero!!\n";
     Error::error(msg.c_str());
 }

 if (head==0) return *this;

 if ((m=g.head->grad)>head->grad) return *this;
 wert=g.head->coef;
 eps=horner(*this,((unity_element<T>()/wert)*g));
 p=NULL;
 for (grad=0,i=0;i<m;i++,grad++) {
     if (eps[i]!=zero) {
        q = new term<T>(eps[i],grad,p);
        p=q;
     }
 }
 h.head=p;
 delete[] eps;
 return h;
}


template <class T>
void divrest(const pol<T>& f, const pol<T>& g,pol<T>& q,pol<T>& r)
{
 q.del();r.del();
 if (g.head==nullptr) {
     std::string msg="\n ERROR : " + val::gettypename(g);
     msg+= ": divrest : Division by zero!!\n";
     Error::error(msg.c_str());
 }

 if (f.head==nullptr) return;
 if ((g.head->grad)>f.head->grad) {
    r = f;
    return;
 }

 int m=g.head->grad,i,grad,n1=f.head->grad;
 const T &wert=g.head->coef;
 T one = unity_element<T>();
 T* eps;
 term<T> *p,*t;
 pol<T> h;

 eps=horner(f,((one/wert)*g));
 p=nullptr;
 for (grad=0,i=m;i<=n1;i++,grad++) {
    if (eps[i]!=pol<T>::zero) {
        t = new term<T>((one/wert)*eps[i],grad,p);
        p=t;
    }
 }
 q.head=p;

 p=nullptr;
 for (grad=0,i=0;i<m;i++,grad++) {
     if (eps[i]!=pol<T>::zero) {
        t = new term<T>(eps[i],grad,p);
        p=t;
     }
 }
 r.head=p;
 delete[] eps;
 return;
}


template <class T>
void divrem(const pol<T>& f, const pol<T>& g,pol<T>& q,pol<T>& r)
{
    divrest(f,g,q,r);
}



template <class T>
const pol<T>& pol<T>::operator /=(const pol<T>& g)
{
 *this = *this / g;
 return *this;
}


template <class T>
const pol<T>& pol<T>::operator %=(const pol<T>& g)
{
 *this = *this % g;
 return *this;
}


template <class T>
const pol<T>& pol<T>::normalize()
{
    if (head==NULL) return *this;
    term<T> *p;
    T div(head->coef);
    if (div == unity_element<T>()) return *this;
    for (p=head;p!=NULL;p=p->next) p->coef/=div;
    return *this;
}

// --------------------------------------------------------------------



// ------------------ Input/Output -----------------------


template <class T>
std::istream& operator >>(std::istream& is,pol<T>& f)
{
 T wert,zero(0);
 //char isch;
 int grad;

 if (f.head != NULL) f.del();

 do {
     wert=zero;
     is>>wert;
     if (wert!= zero) {
        grad=0;
        is>>grad;
        f.insert(wert,grad);
     }
 }
 while (wert!=zero);
 return is;
}


template <class T>
std::ostream& operator <<(std::ostream& os,const pol<T>& f)
{
 term<T>* q;
 for (q=f.head;q!=NULL;q=q->next) os<<q->coef<<"  "<<q->grad<<std::endl;
 os<<0;
 return os;
}

//   -------------------------------------------------------------------------------
 template <class T>
 int polIterator<T>::moveactualtolast()
 {
     if (actual==NULL) return 0;
     while (actual->next!=NULL) actual=actual->next;
     return 1;
 }

// =============================================================================

template <class T>
pol<T> gcd(const pol<T>& a,const pol<T>& b)
{
 pol<T> akt,rest,vorrest;

 if (deg(a)<deg(b)) {akt=b;rest=a;}
 else {akt=a;rest=b;}
 while (!rest.iszero())
 {
       vorrest=std::move(rest);
       rest= akt % vorrest;
       akt= std::move(vorrest);
 }
 akt.normalize();
 return akt;
}



//

template <class T>
void set_unity_element(pol<T> &one)
{
    T one_T;
    set_unity_element(one_T);
    one = pol<T>(one_T);
}


// =============================================================================

} // end namespace val

// ---------------------------------------------------------------------


#endif
