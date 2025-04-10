// Simple linked generic list

#ifndef GLIST_H
#define GLIST_H

//#include <cstdlib>
#include <iostream>
#include <initializer_list>
#include <error.h>
#include <val_basics.h>
#include <val_utils.h>

namespace val
{



template <class T> class GlistIterator;
template <class T> class GlistManipulator;



template <class T>
class Glist
{
protected:
    struct node{
        T value;
        node *next;
        node(const T& v, node *nex = nullptr) : value(v), next(nex) {}
        node(T&& v,node *nex = nullptr) : value(std::move(v)), next(nex) {}
     };

     node *head=nullptr;
     node *actual=nullptr;
     node *last=nullptr;
     int anz=0;                                 // number of nodes (from the German word 'Anzahl')

     static void merge(node* &basis1,node* &basis2);
     static void sort(node* &basis,int l);
public:
    Glist() = default;
    Glist(const Glist<T>&);
    Glist(Glist<T>&&);
    Glist(std::initializer_list<T> args);
    ~Glist() {dellist();}
    //
    const Glist<T>& operator =(const Glist<T>&);
    const Glist<T>& operator =(Glist<T>&&);
    //
    void dellist();                             // deletes whole list
    void skiphead();                            // deletes first node
    void skipnextelement();                     // deletes next node to actual
    int delelement(int);
    void inserttoend(const T&);                 // inserts element at the end;
    void inserttohead(const T&);                // inserts element at the beginning;
    void insertnexttoactual(const T&);
    void sinsert(const T&);                     // ordered insertion in increasing order. Operator < in class T must be defined.
    void inserttoend(T&&);
    void inserttohead(T&&);
    void insertnexttoactual(T&&);
    void push(const T& value) {inserttohead(value);}
    void push(T&& value) {inserttohead(std::move(value));}
    void push_back(const T& value) {inserttoend(value);}
    void push_back(T&& value) {inserttoend(std::move(value));}
    void sinsert(T&&);
    const T& getelement() const;                    // returns value of actual
    const T& getnextelement() const;            // returns value of  actual->next.
    T& actualvalue();                               // returns value of actual
    const T& actualvalue() const;
    void resetactual() {actual=head;}           // Sets actual to head
    void resetlast();                           // Sets last to last element of list
    void moveactual();                          // Moves actual one position forward.
    void setnexttohead();                       // Sets successor of actual to the head of the list.
    int actualvalid() const {if (actual==nullptr) return 0; else return 1;}
    int nextvalid() const {if (actual==nullptr) return 0; else if (actual->next==nullptr) return 0; else return 1;}
    int isempty() const {return (head==nullptr);}
    void copyanddelete(Glist<T>& G);             // equivalent to *this = std::move(G)
    void appendtoend(Glist<T>& G);               // equivalent to append(std::move(G))
    void movenextelementtohead(Glist<T>& G);     // node next to actual is set at the beginning of G
    void moveheadtoend(Glist<T>& G);             // head is set at the end of G
    void append(const Glist<T>& G);
    void append(Glist<T>&& G);
    T& operator[](int);
    const T& operator[](int) const;
    const T& operator() (int i) const {return operator[](i);}
    int length() const {return anz;}
    void sort(int l=-1);
    void merge(Glist<T> &H);                     // list is merged with H. Afterwards H is empty.
    GlistIterator<T> begin() const {GlistIterator<T> It;It.actual=head;return It;}
    GlistIterator<T> end() const {GlistIterator<T> It;It.actual=nullptr;return It;}
    GlistManipulator<T> begin() {GlistManipulator<T> It;It.actual=head;return It;}
    GlistManipulator<T> end() {GlistManipulator<T> It;It.actual=nullptr;return It;}
    friend class GlistIterator<T>;
    friend class GlistManipulator<T>;
};



// Read-Only-Iterator to Glist. Dangling pointer is possible.
template <class T>
class GlistIterator
{
private:
    struct Glist<T>::node *actual=nullptr;
public:
    GlistIterator() = default;
    GlistIterator(const GlistIterator<T>&) = default;
    explicit GlistIterator(const Glist<T> &G) {actual = G.head;}
    void settohead(const Glist<T> &G);
    void settolast(const Glist<T> &G) {actual=G.last;}
    void settoactual(const Glist<T> &G) {actual=G.actual;}
    operator int() const {if (actual==nullptr) return 0;else return 1;}
    int actualvalid() const {if (actual==nullptr) return 0; else return 1;}
    int nextvalid() const {if (actual==nullptr) return 0;if (actual->next==nullptr) return 0; else return 1;}
    int ispointtohead(const Glist<T> &G) const {return (actual==G.head);}
    int ispointtoactual(const Glist<T> &G) const {return (actual==G.actual);}
    int ispointtolast(const Glist<T> &G) const {return (actual==G.last);}
    void moveactual();
    const T& getelement() const;
    const T& actualvalue() const;
    const T& operator() (void) const;
    const T& operator* () const {return actual->value;}
    int operator ==(const GlistIterator<T>& It) const {return actual==It.actual;}
    int operator !=(const GlistIterator<T>& It) const {return actual!=It.actual;}
    void operator ++(int) {if (actual==nullptr) return;actual=actual->next;}
    void operator ++() {if (actual==nullptr) return;actual=actual->next;}
    GlistIterator<T>& operator =(const GlistIterator<T>&);
    GlistIterator<T>& operator =(const Glist<T>& G) {actual=G.head;return *this;}
    friend class Glist<T>;
};



//Read-Write-Iterator zu Glist
template <class T>
class GlistManipulator
{
private:
    struct Glist<T>::node *actual=nullptr;
public:
    GlistManipulator() = default; //{actual=NULL;}
    explicit GlistManipulator(Glist<T> &G) {actual = G.head;}
    void settohead(Glist<T> &G) {actual=G.head;}
    void settolast(Glist<T> &G) {actual=G.last;}
    void settoactual(Glist<T> &G) {actual=G.actual;}
    operator int() const {if (actual==nullptr) return 0;else return 1;}
    int actualvalid() const {if (actual==nullptr) return 0; else return 1;}
    int nextvalid() const {if (actual==nullptr) return 0;if (actual->next==nullptr) return 0; else return 1;}
    int ispointtohead(Glist<T> &G) const {return (actual==G.head);}
    int ispointtoactual(Glist<T> &G) const {return (actual==G.actual);}
    int ispointtolast(Glist<T> &G) const {return (actual==G.last);}
    void moveactual();
    const T& getelement() const;
    T& actualvalue() const;
    T& operator() (void) const;
    T& operator* () const {return actual->value;}
    int operator ==(const GlistManipulator<T>& It) const {return actual==It.actual;}
    int operator !=(const GlistManipulator<T>& It) const {return actual!=It.actual;}
    void operator ++(int) {if (actual==nullptr) return;actual=actual->next;}
    void operator ++() {if (actual==nullptr) return;actual=actual->next;}
    GlistManipulator<T>& operator =(const GlistManipulator<T>&);
    GlistManipulator<T>& operator =(Glist<T>& G) {actual=G.head;return *this;}
    friend class Glist<T>;
};




// Definitions
//=======================================================================================================================================================

//--------------------------------------------------------------------------------------------

// Required function for merge-sort. Operator < must be defined in class T.

template <class T>
void Glist<T>::merge(node* &basis1,node* &basis2)
{
 node *p1,*p2,*r;

 p1=basis1;p2=basis2;
 if (p2->value<basis1->value) {
     while ((p2->next!=nullptr) && (p2->next->value<basis1->value)) p2=p2->next;
     basis1=basis2;
     basis2=p2->next;
     p2->next=p1;
 }
 // jetzt p2<p1
 if (basis2==nullptr) return;
 p2=basis2;
 while (p1->next!=nullptr) {
     if ((basis2!=nullptr) && (basis2->value<p1->next->value)) {
         while ((p2->next!=nullptr) && (p2->next->value<p1->next->value)) p2=p2->next;
         r=p1->next;
         p1->next=basis2;
         basis2=p2->next;
         p2->next=r;
         p1=r;p2=basis2;;
     }
     else p1=p1->next;
 }
 p1->next=basis2;
 basis2=nullptr;
 return;
}


template <class T>
void Glist<T>::sort(node* &basis,int l)
{
 int l1,l2,i;
 node *p1,*p2;

 if (l<=1) return;
 l1=l/2; l2=l-l1;
 p1=basis;
 for (i=1;i<l1;i++) p1=p1->next;
 p2=p1->next;
 p1->next=nullptr;
 p1=basis;
 sort(p1,l1);sort(p2,l2);
 merge(p1,p2);
 basis=p1;
 p1=p2=nullptr;
 return;
}

//--------------------------------------------------------------------------------------------------



template <class T>
Glist<T>::Glist(Glist<T>&& H)
{
    head=H.head; last=H.last; actual = H.actual,anz=H.anz;
    H.head=H.actual=H.last=nullptr;H.anz=0;
}

template <class T>
Glist<T>::Glist(const Glist<T>& H)
{
    anz=0;head=actual=last=nullptr;
    node *p=head;
    for (p=H.head;p!=nullptr;p=p->next) inserttoend(p->value);
}

template <class T>
Glist<T>::Glist(std::initializer_list<T> args)
{
    head=actual=last=nullptr;
    const T *wert;
    for (wert=args.begin();wert!=args.end();wert++) inserttoend(*wert);
}


template <class T>
const Glist<T>& Glist<T>::operator =(Glist<T>&& H)
{
    if (head==H.head) return *this;
    dellist();
    head = H.head; actual=H.actual; last = H.last;anz=H.anz;
    H.head=nullptr; H.actual=nullptr;H.last=nullptr;
    H.anz=0;
    return *this;
}


template <class T>
const Glist<T>& Glist<T>::operator =(const Glist<T>& H)
{
    if (head==H.head) return *this;
    Glist<T> help(H);

    val::swap(help.head,head);
    actual=help.actual; last=help.last; anz=help.anz;
    return *this;
}



template <class T>
void Glist<T>:: dellist()
{
 node *p;
 while (head!=nullptr) {
       p=head;
       head=head->next;
       delete p;
 }
 head=nullptr;actual=head;last=head;anz=0;
}


template <class T>
void Glist<T>::skiphead()
{
    if (head==nullptr) return;
    node *p;

    p=head;
    if (actual==head) actual=actual->next;
    head=head->next;
    delete p;
    if (head==nullptr) last=head;
    --anz;
}


template <class T>
void Glist<T>::skipnextelement()
{
    if (actual==nullptr) return;
    if (actual->next==nullptr) return;

    if (actual->next==last) last=actual;

    node *p=actual->next;
    actual->next=p->next;
    delete p;
    --anz;
}


template <class T>
int Glist<T>::delelement(int i)
{
    if (head==nullptr || i<0) return 0;
    if (i==0) {
        skiphead();
        return 0;
    }
    node *p;
    int j;

    for (p=head,j=0;p->next!=nullptr;p=p->next,++j) {
        if (j==i-1) {
            if (actual==p->next) actual=p;
            if (last==p->next) last=p;
            node *r = p->next;
            p->next = r->next;
            delete r;
            --anz;
            return 1;
        }
    }
    return 0;
}


template <class T>
void Glist<T>:: inserttoend(const T &val)
{
 node *p;

 ++anz;
 if (head==nullptr) {
    head= new node(val);
    actual=head;last=head;
    return;
 }
 p= new node(val);
 last->next=p;
 last=last->next;
}


template <class T>
void Glist<T>:: inserttohead(const T &val)
{
 node *p;

 if (head==nullptr) {
     head= new node(val);
     actual=last=head;
     ++anz;
     return;
 }

 p=new node(val,head);
 head=p;
 ++anz;
}


template <class T>
void Glist<T>:: insertnexttoactual(const T &val)
{
  ++anz;
  node *p = new node(val);
  if (head==nullptr) {head=p;actual=last=head;}
  else if (actual==nullptr) {last->next=p;actual=last=last->next;}
  else {
      p->next = actual->next;
      actual->next=p;
      if (last==actual) last = last->next;
  }

}


template <class T>
void Glist<T>::sinsert(const T& s)
{
 node *p,*r;

 ++anz;
 if (head==nullptr) {
    head= new node(s);
    actual=head;last=head;
    return;
 }
 r=new node(s);
 if (s<head->value) {
     r->next=head;
     head=r;
     return;
 }
 if (head->next==nullptr) {
     head->next=r;
     r->next=nullptr;
     last=r;
     return;
 }
 if (last->value<s) {
     last->next=r;
     r->next=nullptr;
     last=last->next;
     return;
 }
 // now: list has at least 2 nodes, s has to be inserted between them.
 for(p=head;(p->next->value)<s;p=p->next);
 r->next=p->next;
 p->next=r;
 return;
}


template <class T>
void Glist<T>:: inserttoend(T &&val)
{
 node *p;
 ++anz;
 if (head==nullptr) {
    head= new node(std::move(val));
    actual=head;last=head;
    return;
 }
 p= new node(std::move(val));
 last->next=p;
 last=last->next;
}


template <class T>
void Glist<T>:: inserttohead(T &&val)
{
 node *p;

 if (head==nullptr) {
     head= new node(std::move(val));
     actual=last=head;
     ++anz;
     return;
 }

 ++anz;
 p=new node(std::move(val),head);
 head=p;
}


template <class T>
void Glist<T>:: insertnexttoactual(T &&val)
{
  ++anz;
  node *p = new node(std::move(val));
  if (head==nullptr) {head=p;actual=last=head;}
  else if (actual==nullptr) {last->next=p;actual=last=last->next;}
  else {
      p->next = actual->next;
      actual->next=p;
      if (last==actual) last = last->next;
  }

}



template <class T>
void Glist<T>::sinsert(T&& s)
{
 node *p,*r;

 ++anz;
 if (head==nullptr) {
    head= new node(std::move(s));
    actual=head;last=head;
    return;
 }
 r=new node(std::move(s));
 if (r->value<head->value) {
     r->next=head;
     head=r;
     return;
 }
 if (head->next==nullptr) {
     head->next=r;
     r->next=nullptr;
     last=r;
     return;
 }
 if (last->value<r->value) {
     last->next=r;
     r->next=nullptr;
     last=last->next;
     return;
 }
 // now: list has at least 2 nodes, s has to be inserted between them.
 for(p=head;(p->next->value)<r->value;p=p->next);
 r->next=p->next;
 p->next=r;
 return;
}



template <class T>
const T& Glist<T>:: getelement() const
{
    if (actual==nullptr) {
        /*
        std::string msg= "\n" + val::gettypename(*this);
        msg +="::getelement(): NULL-POINTER!!";
        Error::error(msg.c_str());
        */
        Error::error("::getelement(): NULL-POINTER!!",*this);
    }
    return actual->value;
}

template <class T>
const T& Glist<T>:: getnextelement() const
{
    if (actual==nullptr || actual->next==nullptr) {
        /*
        std::string msg = "\n" + val::gettypename(*this);
        msg += "::getnextelement(): NULL-POINTER!!";
        Error::error(msg.c_str());
        */
        Error::error("::getnextelement(): NULL-POINTER!!",*this);
    }
    return actual->next->value;
}

template <class T>
T& Glist<T>::actualvalue()
{
    if (actual==nullptr) {
        /*
        std::string msg = "\n" + val::gettypename(*this);
        msg += "::actualvalue(): NULL-POINTER!!";
        Error::error(msg.c_str());
        */
        Error::error("::actualvalue(): NULL-POINTER!!",*this);
    }
    return actual->value;
}


template <class T>
const T& Glist<T>::actualvalue() const
{
    if (actual==nullptr) {
        /*
        std::string msg = "\n" + val::gettypename(*this);
        msg += "::actualvalue(): NULL-POINTER!!";
        Error::error(msg.c_str());
        */
        Error::error("::actualvalue(): NULL-POINTER!!",*this);
    }
    return actual->value;
}



template <class T>
void Glist<T>:: resetlast()
{

 if (head==nullptr) {
     last=head;
     return;
 }
 for (last=head;last->next!=nullptr;last=last->next);
}


template <class T>
void Glist<T>:: moveactual()
{
    if (actual==nullptr) {
        /*
        std::string msg = "\n" + val::gettypename(*this);
        msg += "::moveactual(): NULL-POINTER!!";
        Error::error(msg.c_str());
        */
        Error::error("::moveactual(): NULL-POINTER!!",*this);
    }
    actual=actual->next;
}


template <class T>
void Glist<T>:: setnexttohead()
{
 node *q;

 if (actual==nullptr || actual->next==nullptr) {
    /*
    std::string msg = "\n" + val::gettypename(*this);
    msg += "::setnexttohead(): NULL-POINTER!!";
    Error::error(msg.c_str());
    */
    Error::error("::setnexttohead(): NULL-POINTER!!",*this);
 }

 q=actual->next;

 if (q==last) last=actual;
 actual->next=q->next;
 q->next=head;
 head=q;
}




// Danach: *this=H, H leer.
template <class T>
void Glist<T>::copyanddelete(Glist<T> &H)
{
 dellist();
 head=H.head;
 actual=H.actual;
 last=H.last;
 anz=H.anz;
 H.head=H.actual=H.last=nullptr;H.anz=0;
 return;
}

// H wird am Ende der Liste Eingefuegt, danach H leer!
template <class T>
void Glist<T>::appendtoend(Glist<T> &H)
{
    if (H.head==nullptr) return;
    if (head==nullptr) {
        head=H.head;actual=head;
    }
    else last->next=H.head;
    last=H.last;
    anz+=H.anz;
    H.head=H.last=H.actual=nullptr;H.anz=0;
    return;
}

template <class T>
void Glist<T>::append(const Glist<T> &G)
{
    //std::cout<<"\nhier!";
    if (G.head==nullptr) return;
    for (node *it=G.head;it!=nullptr;it=it->next) {
        //std::cout<<"\nhier!";
        push_back(it->value);
    }
}

template <class T>
void Glist<T>::append(Glist<T>&& G)
{
    if (G.head==nullptr) return;
    if (head==nullptr) {
        actual=head=G.head;
        last=G.last;
        anz=G.anz;
        G.head=G.actual=G.last=nullptr;
        G.anz=0;
    }
    else {
        last->next=G.head;
        last=G.last;
        anz+=G.anz;
    }
    G.head=G.actual=G.last=nullptr;
    G.anz=0;
}




template <class T>
void Glist<T>::movenextelementtohead(Glist<T> &G) // Glied nach actual wird an den  Kopf von G gesetzt;
{
    if (isempty()) return;
    if (actual==nullptr) return;
    if (actual->next==nullptr) return;

    node *p=actual->next;

    if (last==actual->next) last = actual;

    actual->next=p->next;
    p->next=G.head;
    G.head=p;
    ++(G.anz);
    --anz;
}


template <class T>
void Glist<T>::moveheadtoend(Glist<T> &G)
{
    if (head==nullptr) return;
    if (actual==head) actual=actual->next;
    if (last==head) last=last->next;

    node *p=head;
    head=head->next;
    --anz;
    p->next=nullptr;

    if (G.head==nullptr) G.head=G.actual=G.last=p;
    else {G.last->next=p;G.last=G.last->next;}
    ++(G.anz);
}


template <class T>
T& Glist<T>:: operator[](int i)
{
 int j=0;
 node *p;

 if (i<0 || i>=anz)  {

     std::string msg = "\n" + val::gettypename(*this);
     msg+=": ERROR: Index out of range!\nIndex: "+ToString(i) + ". Size: " + ToString(anz)+"\n";
     Error::error(msg.c_str());
 }

 for (p=head;(p!=nullptr)&&(j<i);p=p->next,j++);

 return (p->value);
}


template <class T>
const T& Glist<T>:: operator[](int i) const
{
 int j=0;
 node *p;


 if (i<0 || i>=anz)  {
     std::string msg = "\n" + val::gettypename(*this);
     msg+=": ERROR: Index out of range!\nIndex: "+ToString(i) + ". Size: " + ToString(anz)+"\n";
     Error::error(msg.c_str());
 }

 for (p=head;(p!=nullptr)&&(j<i);p=p->next,j++);

 return (p->value);
}



template <class T>
void Glist<T>::sort(int l)
{
 if (l==-1) l=length();
 if (l<=1) return;
 sort(head,l);

 // l>1
 for (last=head;last->next!=NULL;last=last->next);

 return;
}


template <class T>
void Glist<T>::merge(Glist<T> &H)
{
    merge(head,H.head);
    H.head=H.actual=H.last=nullptr;
    resetactual();
    resetlast();
}


// ----------------------------------------------------------------------------------------

template <class T>
void GlistIterator<T>::settohead(const Glist<T>& G)
{
    actual=G.head;
}


template <class T>
void GlistIterator<T>::moveactual()
{

    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::moveactual(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistIterator<T>::moveactual(): NULL-POINTER!!");
        Error::error("::moveactual(): NULL-POINTER!!",*this);
    }
    actual=actual->next;
}

template <class T>
const T& GlistIterator<T>:: getelement() const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::getelement(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistIterator<T>::getlement(): NULL-POINTER!!");
        Error::error("::getlement(): NULL-POINTER!!",*this);

    }
    return actual->value;
}



template <class T>
const T& GlistIterator<T>:: actualvalue() const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::actualvalue(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistIterator<T>::actualvalue(): NULL-POINTER!!");
        Error::error("::actualvalue(): NULL-POINTER!!",*this);
    }

    return actual->value;
}


template <class T>
const T& GlistIterator<T>::operator() (void) const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::operator(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistIterator<T>::operator(): NULL-POINTER!!");
        Error::error("::operator(): NULL-POINTER!!",*this);

    }

    return actual->value;
}



template<class T>
GlistIterator<T>& GlistIterator<T>::operator =(const GlistIterator<T> &It)
{
    actual = It.actual;
    return *this;
}



// ------------------------------------------------------------------------------

template <class T>
void GlistManipulator<T>::moveactual()
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::moveactual(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistManipulator<T>::moveactual(): NULL-POINTER!!");
        Error::error("::moveactual(): NULL-POINTER!!",*this);
    }

    actual=actual->next;
}

template <class T>
const T& GlistManipulator<T>:: getelement() const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::getelement(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistManipulator<T>::getlement(): NULL-POINTER!!");
        Error::error("getlement(): NULL-POINTER!!",*this);
    }
    return actual->value;
}


template <class T>
T& GlistManipulator<T>:: actualvalue() const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::actualvalue(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistManipulator<T>::actualvalue(): NULL-POINTER!!");
        Error::error("::actualvalue(): NULL-POINTER!!",*this);
    }
    return actual->value;
}


template <class T>
T& GlistManipulator<T>::operator() (void) const
{
    if (actual==nullptr) {
        //std::string msg = "\n" + val::gettypename(*this);
        //msg += "::operator(): NULL-POINTER!!";
        //Error::error(msg.c_str());
        //Error::error("\nval::GlistManipulator<T>::getlement(): NULL-POINTER!!");
        Error::error("::getlement(): NULL-POINTER!!",*this);
    }
    return actual->value;
}


template<class T>
GlistManipulator<T>& GlistManipulator<T>::operator =(const GlistManipulator<T> &It)
{
    actual = It.actual;
    return *this;
}



} // end namespace val


#endif //GLIST_H
