#ifndef SET_H_INCLUDED
#define SET_H_INCLUDED

#include <iostream>
#include <initializer_list>
#include <error.h>
#include <val_basics.h>
#include <val_utils.h>

namespace val
{

template <class T> class SetIterator;

template <class T>
class Set
{
protected:
    struct node {
        T value;
        node *prev,*next;
        node (const T& coef,node* pre=nullptr,node* nex=nullptr) : value(coef) , prev(pre) , next(nex) {}
        node (T&& coef,node* pre=nullptr,node* nex=nullptr) : value(std::move(coef)) , prev(pre) , next(nex) {}
    };
    node *head=nullptr;
    node *last=nullptr;
    int anz=0;
public:
    Set() = default;
    Set(const Set<T> &G);
    Set(Set<T>&& G);
    Set(std::initializer_list<T> args);
    ~Set() {del();}
    //
    int isempty() const {return head==nullptr;}
    int length() const {return anz;}
    //
    const T& operator[] (int n) const;
    //
    const Set<T>& operator=(const Set<T>& G);
    const Set<T>& operator=(Set<T>&& G);
    void del();
    void sinsert(const T& value);
    void sinsert(T&& value);
    void insert(const T& value) {sinsert(value);}
    void insert(T&& value) {sinsert(std::move(value));}
    //
    SetIterator<T> begin() const {return SetIterator<T>(*this);}
    SetIterator<T> end() const {return SetIterator<T>();}
    SetIterator<T> endofset() const {SetIterator<T> It;It.actual=last;return It;}

    friend class SetIterator<T>;
};


template <class T>
class SetIterator
{
private:
    struct Set<T>::node *actual = nullptr;
public:
    SetIterator() = default;
    SetIterator(const Set<T>& G) {actual=G.head;}
    SetIterator(const SetIterator<T> &It) = default;
    operator int() const {if (actual==nullptr) return 0; else return 1;}
    const T& operator() (void) const {return actual->value;}
    const T& operator* () const {return actual->value;}
    void operator ++(int) {if (actual) actual=actual->next;}
    void operator ++() {if (actual) actual=actual->next;}
    void operator --(int) {if (actual) actual=actual->prev;}
    void operator --() {if (actual) actual=actual->prev;}
    int operator ==(const SetIterator<T> &It) const {return actual==It.actual;}
    int operator !=(const SetIterator<T> &It) const {return actual!=It.actual;}
    friend class Set<T>;
};




template <class T>
Set<T>::Set(const Set<T> &G)
{
    head=last=nullptr; anz=0;
    if (G.head==nullptr) return;
    // setze head;
    node *p,*q,*r;
    head = new node(G.head->value); last=head;// last->next = nullptr;last->prev = nullptr;
    anz++;
    for (p=head,q=G.head->next;q !=nullptr;p=p->next,q=q->next) {
        r = new node(q->value,p,p->next);
        p->next=r;last = r;
        anz++;
    }
}


template <class T>
Set<T>::Set(Set<T> &&G)
{
    head=G.head;last=G.last;anz=G.anz;
    G.head = G.last = nullptr;
    G.anz=0;
}

template <class T>
Set<T>::Set(std::initializer_list<T> args)
{
    const T *wert;
    head=last=nullptr; anz=0;
    for (wert=args.begin();wert!=args.end();wert++) sinsert(*wert);
}



template <class T>
void Set<T>::del()
{
    node *p;
    while (head!=nullptr) {
        p=head;
        head=head->next;
        delete p;
    }
    head=last=nullptr;
    anz=0;
}


template <class T>
const T& Set<T>::operator[] (int n) const
{
 int i=0;
 node *p;

 if (n<0 || n>=anz)  {
	 std::string msg=val::gettypename(*this);
	 msg+=": ERROR: Index out of range!\nIndex: "+ToString(n) + ". Size: " + ToString(anz)+"\n";
	 Error::error(msg.c_str());
 }

 for (p=head;(p!=NULL)&&(i<n);p=p->next,i++);

 //if (p==NULL) Error::error("\nGlist::operator[]: NULL-POINTER!!");
 return (p->value);

}


template <class T>
const Set<T>& Set<T>::operator=(const Set<T> &G)
{
    if (head==G.head) return *this;
    Set<T> help(G);

    val::swap(help.head,head);
    last=help.last; anz=help.anz;
    return *this;

    /*
    node *p,*q,*r;
    del();
    if (G.head==nullptr) return *this;
    head = new node(G.head->value);
    last = head;
    anz++;
    for (p=G.head->next,q=head;p!=nullptr;p=p->next,q=q->next,++anz) {
        //std::cout<<"n anz = "<<anz;
        r = new node(p->value,q);
        last = r;
        q->next = r;
        //q = q->next;
    }
    return *this;
    */
}

template <class T>
const Set<T>& Set<T>::operator=(Set<T> &&G)
{
    if (head==G.head) return *this;
    del();
    head = G.head; last = G.last; anz = G.anz;
    G.head = G.last = nullptr; G.anz=0;
    return *this;
}



template <class T>
void Set<T>::sinsert(const T& coef)
{
    if (head==nullptr) {
        head = new node(coef);
        last = head;
        anz++;
        return;
    }

    if (coef < head->value) {
        node *r = new node(std::move(coef),nullptr,head);
        head->prev = r;
        head = r;
        anz++;
        return;
    }
    if (coef == head->value) return;

    if (coef > last->value) {
        node *r = new node(std::move(coef),last);
        last->next = r;
        last = last->next;
        anz++;
        return;
    }
    if (coef == last->value) return;

    node *p,*r;

    for (p=head;p->next!=nullptr;p=p->next) {
        if (p->next->value == coef) return;
        if (coef < p->next->value) {
            r = new node(std::move(coef),p,p->next);
            p->next = r;
            return;
        }
    }
    return;
}


template <class T>
void Set<T>::sinsert(T&& coef)
{
    if (head==nullptr) {
        head = new node(std::move(coef));
        last = head;
        anz++;
        return;
    }

    if (coef < head->value) {
        node *r = new node(std::move(coef),nullptr,head);
        head->prev = r;
        head = r;
        anz++;
        return;
    }
    if (coef == head->value) return;

    if (coef > last->value) {
        node *r = new node(std::move(coef),last);
        //last= last->next;
        last->next = r;
        last = last->next;
        anz++;
        return;
    }
    if (coef == last->value) return;

    node *p;

    for (p=head;p->next!=nullptr;p=p->next) {
        if (p->next->value == coef) return;
        if (coef < p->next->value) {
            node *r = new node(std::move(coef),p,p->next);
            p->next = r;
            return;
        }
    }
    return;
}


} //end namespace val

#endif // SET_H_INCLUDED
