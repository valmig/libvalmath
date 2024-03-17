#ifndef D_ARRAY_H_INCLUDED
#define D_ARRAY_H_INCLUDED


#include <val_basics.h>
#include <error.h>
#include <initializer_list>
#include <val_utils.h>

namespace val
{

template <class T> class d_array;



template <class T>
class d_array
{
private:
    T *dat=nullptr;
    int len=0;
    int cap=0;
    static int plus_cap;
    static int relation_order(const T& a,const T& b) {if (increasing_order) return a<b; else return a>b;}
public:
    static int increasing_order;
    d_array() = default;
    explicit d_array(int n);
    d_array(const T &val,int n);
    d_array(const d_array<T>&);
    d_array(d_array<T>&&);
    d_array(std::initializer_list<T> args);
    const d_array<T>& operator =(const d_array<T>&);
    const d_array<T>& operator =(d_array<T>&&);
    ~d_array() {if (dat!=nullptr) {delete[] dat;cap=len=0;}}
    void del() {if (dat!=nullptr) {delete[] dat;cap=len=0;dat=nullptr;}}
    int length() const {return len;}
    int capacity() const {return cap;}
    int isempty() const {return (len==0);}
    void reserve(int n);
    void resize(int n);
    static void set_plus_cap(int n);
    T& operator [](int i);
    const T& operator [](int i) const;
    void inserttoend(const T& value);
    void inserttoend(T&& value);
    void push_back(const T& value) {inserttoend(value);}
    void push_back(T&& value) {inserttoend(std::move(value));}
    void sinsert(const T& value);                             // ordered insertion resp. relation_order. Operators < , > must be defined in class T.
    void sinsert(T&& value);
    void append(const d_array<T>& v);
    void append(d_array<T>&& v);
    d_array<T> operator+ (const d_array<T> &u) const;
    const T* begin() const {return dat;}
    T* begin() {return dat;}
    const T* end() const {if (dat==nullptr) return dat; else return &dat[len];}
    T* end()  {if (dat==nullptr) return dat; else return &dat[len];}
    void sort();
    void sort(int left,int rigth);
};



template <class T>
int d_array<T>::plus_cap=0;

template <class T>
int d_array<T>::increasing_order=1;

template <class T>
d_array<T>::d_array(int n)
{
    if (n<=0) return;
    dat=new T[n+plus_cap];
    len=n;
    cap=n+plus_cap;
}

template <class T>
d_array<T>::d_array(const T& val,int n)
{
    if (n<=0) return;
    dat=new T[n+plus_cap];
    len=n;
    cap=n+plus_cap;
    for (int i=0;i<n;++i) dat[i]=val;
}

template <class T>
d_array<T>::d_array(const d_array<T> &A)
{
    if (A.dat==nullptr) return;
    dat=new T[cap=(A.len + plus_cap)];
    len=A.len;
    for (int i=0;i<len;++i) dat[i]=A.dat[i];
}


template <class T>
d_array<T>::d_array(d_array<T> &&A)
{
    dat=A.dat;
    len=A.len;
    cap=A.cap;
    A.dat=nullptr;A.len=A.cap=0;
}


template <class T>
d_array<T>::d_array(std::initializer_list<T> args)
{
    len=args.size();
    if (len<=0) return;
    cap=len+plus_cap;
    dat=new T[cap];
    int i=0;
    for (const T& v : args) {
        dat[i] = v;
        ++i;
    }
}


template <class T>
const d_array<T>& d_array<T>::operator =(const d_array<T>& A)
{
    if (dat==A.dat) return *this;
    if (A.len<=cap) {
        for (int i=0;i<A.len;++i) dat[i]=A.dat[i];
        len=A.len;
        return *this;
    }
    else {
        //if (dat!=NULL) delete[] dat;
        cap=A.len + plus_cap;
        len=A.len;
        T* help= new T[cap];
        for (int i=0;i<len;++i) help[i]=A.dat[i];
        if (dat!=nullptr) delete[] dat;
        dat=help;
        help=nullptr;
        return *this;
    }
}


template <class T>
const d_array<T>& d_array<T>::operator =(d_array<T>&& A)
{
    if (dat==A.dat) return *this;
    if (dat!=nullptr) delete[] dat;
    dat=A.dat; cap=A.cap; len=A.len;
    A.dat=nullptr; A.len=A.cap=0;
    return *this;
}


template <class T>
void d_array<T>::reserve(int n)
{
    if (dat!=nullptr || n<=0) return;
    dat= new T[cap=n];
    len=0;
}


template <class T>
void d_array<T>::resize(int n)
{
    if (n<0) return;
    if (n<=cap) len=n;
    else {
        cap=n+plus_cap;
        T *h,*hdat = new T[cap];
        for (int i=0;i<len;++i) hdat[i] = std::move(dat[i]);
        len = n;
        h=dat;
        dat=hdat;
        hdat = nullptr;
        delete[] h;
    }
    return;
}

template <class T>
void d_array<T>::set_plus_cap(int n)
{
    if (n<0) return;
    plus_cap=n;
}


template <class T>
T& d_array<T>::operator [](int i)
{
    if (i<0 || i >= cap) {  //Error::error("\nERROR: d_array<T>::operator[] : NULL-Pointer!");
		std::string msg=val::gettypename(*this);
		msg+=": ERROR: Index out of range!\nIndex: "+ToString(i) + ". Capacity: " + ToString(cap)+"\n";
		Error::error(msg.c_str());
	}
    else if (i>=len) len=i+1;
    return dat[i];
}


template <class T>
const T& d_array<T>::operator [](int i) const
{
    if (i<0 || i >= len) { //Error::error("\nERROR: d_array<T>::operator[] : NULL-Pointer!");
		std::string msg=val::gettypename(*this);
		msg+=": ERROR: Index out of range!\nIndex: "+ToString(i) + ". Size: " + ToString(len)+ ". Capacity: "+ ToString(cap)+ "\n";
		Error::error(msg.c_str());
	}
    return dat[i];
}


template <class T>
void d_array<T>::inserttoend(const T& value)
{
    if (cap>len) {
        dat[len] = value;
        ++len;
    }
    else {
        cap+=plus_cap+1;
        T *hdat = new T[cap];
        for (int i=0;i<len;++i) hdat[i] = std::move(dat[i]);
        hdat[len] = value;
        ++len;
        delete[] dat;
        dat=hdat;
        hdat = nullptr;
    }
}


template <class T>
void d_array<T>::inserttoend(T&& value)
{
    if (cap>len) {
        dat[len] = std::move(value);
        ++len;
    }
    else {
        cap+=plus_cap+1;
        T *hdat = new T[cap];
        for (int i=0;i<len;++i) hdat[i] = std::move(dat[i]);
        hdat[len] = std::move(value);
        ++len;
        delete[] dat;
        dat=hdat;
        hdat = nullptr;
    }
}


template <class T>
void d_array<T>::sinsert(const T& value)
{
    if (cap>len) {
        if (len==0) {dat[0] = value; ++len; return;}
        int i,stelle=len;
        for (i=0;i<len;++i) {
            if (relation_order(value,dat[i])) {
                stelle=i;
                break;
            }
        }
        for (i=len;i>stelle;--i) dat[i]=std::move(dat[i-1]);
        dat[stelle] = value;
        ++len;
    }
    else {
        cap+=plus_cap+1;
        if (len == 0) {dat = new T[cap];dat[0]=value;++len;return;}
        T *hdat = new T[cap];
        int i,stelle=len;
        for (i=0;i<len;++i) {
            if (relation_order(value,dat[i])) {
                stelle=i;
                break;
            }
        }
        for (i=0;i<stelle;++i) hdat[i]=std::move(dat[i]);
        hdat[stelle] = value;
        for (i=stelle;i<len;++i) hdat[i+1] = std::move(dat[i]);
        ++len;
        delete[] dat;
        dat = hdat;
        hdat=nullptr;
    }
}


template <class T>
void d_array<T>::sinsert(T&& value)
{
    if (cap>len) {
        if (len==0) {dat[0] = std::move(value); ++len; return;}
        int i,stelle=len;
        for (i=0;i<len;++i) {
            if (relation_order(value,dat[i])) {
                stelle=i;
                break;
            }
        }
        for (i=len;i>stelle;--i) dat[i]=std::move(dat[i-1]);
        dat[stelle] = std::move(value);
        ++len;
    }
    else {
        cap+=plus_cap+1;
        if (len == 0) {dat = new T[cap];dat[0]=std::move(value);++len;return;}
        T *hdat = new T[cap];
        int i,stelle=len;
        for (i=0;i<len;++i) {
            if (relation_order(value,dat[i])) {
                stelle=i;
                break;
            }
        }
        for (i=0;i<stelle;++i) hdat[i]=std::move(dat[i]);
        hdat[stelle] = std::move(value);
        for (i=stelle;i<len;++i) hdat[i+1] = std::move(dat[i]);
        ++len;
        delete[] dat;
        dat = hdat;
        hdat=nullptr;
    }
}


template <class T>
void d_array<T>::append(const d_array<T>& v)
{
    int l=len;
    resize(len + v.len);
    for (int i=l;i<len;++i) dat[i] = v.dat[i-l];
}

template <class T>
void d_array<T>::append(d_array<T>&& v)
{
    int l=len;
    resize(len + v.len);
    for (int i=l;i<len;++i) dat[i] = std::move(v.dat[i-l]);
    delete[] v.dat;
    v.len=v.cap=0;
    v.dat=nullptr;
}


template <class T>
d_array<T> d_array<T>::operator+ (const d_array<T> &u) const
{
    d_array<T> v;
    int i,j;
    v.cap = u.cap + cap;
    v.len = u.len + len;
    if (!v.cap) return v;

    v.dat = new T[v.cap];

    for (i = 0; i < cap; ++i) v.dat[i] = dat[i];
    for (j = 0; j < u.cap; ++i, ++j) v.dat[i] = u.dat[j];

    return v;
}


template <class T>
void d_array<T>::sort(int left,int right)
{
    if (right>=len) right=len-1;
    if (left<0) left=0;
    if (left>=right) return;

    int i,last;

    val::swap(dat[left],dat[(left+right)/2]);
    last=left;
    for (i=left+1;i<=right;i++) {
        if (relation_order(dat[i],dat[left])) val::swap(dat[++last],dat[i]);
    }
    val::swap(dat[left],dat[last]);
    sort(left,last-1);
    sort(last+1,right);
}

template <class T>
void d_array<T>::sort()
{
    if (len==0) return;
    sort(0,len-1);
}


}
#endif // D_ARRAY_H_INCLUDED
