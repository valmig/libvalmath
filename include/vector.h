#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED

#include <iostream>
//#include <cstdlib>
#include <error.h>
#include <val_basics.h>
#include <val_utils.h>
#include <initializer_list>


// val::vector<T> is a class for representing elements of a module of finite rank (e.g. finite dimensional vector spaces).
// So T should be ring. 
// val::vector<T> should not be used as a container class (use val::d_array<T> instead). 


namespace val
{



// Friend-functions of vector -----------------------------------------
template <class T> class vector;
template <class T> class matrix;
template <class T> vector<T> operator*(const T&,const vector<T>&);
template <class T> std::istream& operator >>(std::istream&,vector<T>&);
template <class T> std::ostream& operator <<(std::ostream&,const vector<T>&);

#ifdef MATRIX_H_INCLUDED
template <class T> const vector<vector<T> >& movefrommatrixtovector(matrix<T>& A,vector<vector<T> >& V,int r); // First r rows
#else
template <class T> const vector<vector<T> >& movefrommatrixtovector(matrix<T>& A,vector<vector<T> >& V,int r=0); // First r rows
#endif // MATRIX_H_INCLUDED

template <class T> const matrix<T>& movefromvectortomatrix(vector<vector<T> >&,matrix<T>&);
// ---------------------------------------------------------------------------------------


template <class T>
class vector
{
protected:
    T* coeff=nullptr;
    int dim=0;
    static const T zero;
    static vector<T> char_to_vector(const char* c,int l);
public:
    vector() = default;
    explicit vector(int);        // vector will not be initialized
    vector (const T& x,int n);   // *this = (x,...,x)
    vector(const vector<T>&);
    vector(vector<T>&&);
    vector(std::initializer_list<T> args);
    ~vector();
    void sort();
    void sort(int left,int rigth);
    const vector<T>& make_zero();
    int dimension() const {return dim;}
    const vector<T>& operator =(const vector<T>&);
    const vector<T>& operator =(vector<T>&&);
    int iszero() const;
    int isempty() const {return (coeff==NULL);}
    int operator ==(const vector<T>&) const;
    int operator !=(const vector<T>& v) const {return !(*this==v);}
    vector<T> operator +(const vector<T>&) const;
    vector<T> operator -(const vector<T>&) const;
    const vector<T>& operator +=(const vector<T>&);
    const vector<T>& operator -=(const vector<T>&);
    vector<T> operator-() const;
    const T& operator [](int i) const {return operator()(i);}
    T& operator [](int i) {return operator()(i);}
    const T& operator ()(int) const;
    T& operator ()(int);
    T* begin() {return coeff;}
    const T* end() const {if (coeff==nullptr) return coeff; else return &coeff[dim];}
    T* end()  {if (coeff==nullptr) return coeff; else return &coeff[dim];}
    friend vector val::operator *<>(const T&,const vector&);
    const vector<T>& operator *=(const T&);
	friend std::istream& operator >><T>(std::istream&,vector<T>&);
	friend std::ostream& operator <<<T>(std::ostream&,const vector<T>&);
	friend const vector<vector<T> >& movefrommatrixtovector<T>(matrix<T>& A,vector<vector<T> >& V,int r);
	friend const matrix<T>& movefromvectortomatrix<T>(vector<vector<T> >&,matrix<T>&);
};



template<class T>
const T vector<T>::zero(val::zero_element<T>());

template <class T>
vector<T>::vector(int n)
{
    if (n<=0) {
        coeff=NULL;
        dim=0;
        return;
    }
    coeff= new T[dim=n];
}


template <class T>
vector<T>::vector(const T& x,int n)
{
    if (n<=0) {
        coeff=NULL;
        dim=0;
        return;
    }
    int i;
    coeff= new T[dim=n];
    for (i=0;i<n;i++) coeff[i]=x;
}


template <class T>
vector<T>::vector(const vector<T>& v)
{
    if (v.coeff==NULL) {
        coeff=NULL;
        dim=0;
        return;
    }
    int i;
    dim=v.dim;
    coeff = new T[dim];
    for (i=0;i<dim;i++) coeff[i]=v.coeff[i];
}

template <class T>
vector<T>::vector(vector<T>&& v)
{
    coeff=v.coeff;
    dim=v.dim;
    v.coeff=NULL;
    v.dim=0;
}

template <class T>
vector<T>::vector(std::initializer_list<T> args)
{
    dim=args.size();
    if (dim<=0) {dim=0;coeff=NULL;return;}
    coeff=new T[dim];
    const T *wert;
    int i;
    for (i=0,wert=args.begin();wert!=args.end();++wert,i++) coeff[i]= *wert;
}


template <class T>
vector<T>:: ~vector()
{
    if (coeff!=NULL) delete[] coeff;
}


template <class T>
const vector<T>& vector<T>::make_zero()
{
    if (coeff==NULL) return *this;
    int i;
    for (i=0;i<dim;i++) coeff[i]=zero;
    return *this;
}


template <class T>
void vector<T>::sort(int left,int right)
{
    if (right>=dim) right=dim-1;
    if (left<0) left=0;
    if (left>=right) return;

    int i,last;

    val::swap(coeff[left],coeff[(left+right)/2]);
    last=left;
    for (i=left+1;i<=right;i++) {
        if (coeff[i]<coeff[left]) val::swap(coeff[++last],coeff[i]);
    }
    val::swap(coeff[left],coeff[last]);
    sort(left,last-1);
    sort(last+1,right);
}

template <class T>
void vector<T>::sort()
{
    if (coeff==NULL) return;
    sort(0,dim-1);
}



template <class T>
const vector<T>& vector<T>::operator =(const vector<T>& v)
{
    if (coeff == v.coeff) return *this;
    if (v.coeff==NULL) {
        delete[] coeff;
        coeff=NULL;
        dim=0;
        return *this;
    }
    if (dim!=v.dim) {
        T *hcoeff = new T[dim=v.dim];
        for (int i=0;i<dim;i++) hcoeff[i] =v.coeff[i];
        if (coeff!=NULL) delete[] coeff;
        coeff=hcoeff;
        hcoeff=nullptr;
        return *this;
    }
    for (int i=0;i<dim;i++) coeff[i] =v.coeff[i];
    return *this;
}


template <class T>
const vector<T>& vector<T>::operator =(vector<T>&& v)
{
    if (coeff==v.coeff) return *this;
    if (coeff!=NULL) delete[] coeff;
    coeff=v.coeff;
    dim=v.dim;
    v.coeff=NULL;
    v.dim=0;
    return *this;
}


template <class T>
int vector<T>::iszero() const
{
    if (coeff==NULL) return 1;
    for (int i=0;i<dim;i++)
        if (coeff[i]!=zero) return 0;
    return 1;
}


template <class T>
int vector<T>::operator ==(const vector<T>& v) const
{
    if (coeff==NULL || v.coeff==NULL) return (coeff==v.coeff);
    if (dim!=v.dim) return 0;
    for (int i=0;i<dim;i++)
        if (coeff[i]!=v.coeff[i]) return 0;
    return 1;
}


template <class T>
vector<T> vector<T>::operator +(const vector<T>& v) const
{
    if (v.coeff==NULL) return *this;
    if (coeff==NULL) return vector<T>(v);

    int i,r;
    if (dim!=v.dim) {
        Error::warning("\nWARNING: vector<T>::operator+ : vectors of different dimensions!");
    }
    if (dim<=v.dim) r=dim;
    else r = v.dim;
    vector<T> w(r);

    for (i=0;i<r;i++) w.coeff[i]=coeff[i]+v.coeff[i];
    return w;
}

template <class T>
vector<T> vector<T>::operator -() const
{
    if (coeff==NULL) return *this;
    vector<T> w(dim);

    for (int i=0;i<dim;i++) w.coeff[i]=-coeff[i];
    return w;
}

template <class T>
vector<T> vector<T>::operator -(const vector<T>& v) const
{
    if (v.coeff==NULL) return *this;
    if (coeff==NULL) return -v;

    int i,r;
    if (dim!=v.dim) {
        Error::warning("\nWARNING: vector<T>::operator- : vectors of different dimensions!");
    }
    if (dim<=v.dim) r=dim;
    else r = v.dim;
    vector<T> w(r);

    for (i=0;i<r;i++) w.coeff[i]=coeff[i]-v.coeff[i];
    return w;
}


template <class T>
const vector<T>& vector<T>::operator +=(const vector<T>& v)
{
    if (v.coeff==NULL) return *this;
    if (coeff==NULL) {
        *this=v;
        return *this;
    }

    int i,r;
    if (dim!=v.dim) {
        Error::warning("\nWARNING: vector<T>::operator+= : vectors of different dimensions!");
    }
    if (dim<=v.dim) r=dim;
    else r=v.dim;
    for (i=0;i<r;i++) coeff[i]+=v.coeff[i];
    return *this;
}


template <class T>
const vector<T>& vector<T>::operator -=(const vector<T>& v)
{
    if (v.coeff==NULL) return *this;
    if (coeff==NULL) {
        *this=-v;
        return *this;
    }

    int i,r;
    if (dim!=v.dim) {
        Error::warning("\nWARNING: vector<T>::operator-= : vectors of different dimensions!");
    }
    if (dim<=v.dim) r=dim;
    else r=v.dim;
    for (i=0;i<r;i++) coeff[i]-=v.coeff[i];
    return *this;
}



template <class T>
const T& vector<T>::operator() (int i) const
{
    if (coeff==NULL || i<0 || i>=dim) return zero;
    else return coeff[i];
}


template <class T>
T& vector<T>::operator() (int i)
{
    if (i<0 || i >= dim) {
		std::string msg=val::gettypename(*this);
		msg+=": ERROR: Index out of range!\nIndex: "+ToString(i) + ". Size: " + ToString(dim)+"\n";
		Error::error(msg.c_str());
	}
    return coeff[i];
}


template<class T>
const vector<T>& vector<T>::operator *=(const T& wert)
{
    if (coeff==NULL) return *this;
    for (int i=0;i<dim;i++) coeff[i]*=wert;
    return *this;
}

template <class T>
vector<T> operator*(const T& wert,const vector<T>& v)
{
    vector<T> w(v);

    for (int i=0;i<w.dim;i++) w.coeff[i]=wert*v.coeff[i];
    return w;
}


template <class T>
vector<T> vector<T>::char_to_vector(const char* s,int l)
{
    vector<T> LES;

    if (l==0) return LES;
    int first=0,i,k,n;
    std::string zeichen="";

    for (i=0;i<l;i++)
        if (s[i]!='\n' && s[i]!=' ' ) break;
    first=i;
    n=0;
    //Bestimme Spaltenzahl n
    for (;i<l;) {
        if (i<l && s[i]=='\n') {break;}
        while (i<l && (s[i]!=' ' && s[i] != '\n' )) i++;
        n++;
        while (i<l && s[i]==' ') i++;
    }

    if (n==0) return LES;

    LES=vector<T>(n);

    i=0;
    for (k=first;k<l;) {
        if (s[k]!=' ' && s[k]!='\n') {
            while (k<l && s[k]!=' ' && s[k]!='\n') {
                zeichen+=s[k];
                k++;
            }
            //LES(i,j)=val::string_to_rational(zeichen);
            LES(i)=val::FromString<T>(zeichen);
            zeichen="";
            i++;
        }
        while (k<l && s[k]==' ') {
               k++;
        }
    }
    return LES;
}


template <class T>
std::istream& operator >>(std::istream& is,vector<T>& v)
{
  char *buf,*hbuf;
  int i,l=0,wert=0,k=1;

  buf =new char[1000];


  while ((wert!=10) && is) {
	wert=is.get();
	if (wert==10){
		if (l!=0) buf[l]='\0';
		else wert=0;
	}
	else buf[l++]=wert;
	if (l>=(k*1000)) {
        k++;
        hbuf= new char[k*1000];
        for (i=0;i<((k-1)*1000);i++) hbuf[i]=buf[i];
        delete[] buf;
        buf=hbuf;
        hbuf=NULL;
	}

	if (!is) l--;
  }

  v=vector<T>::char_to_vector(buf,l);
  delete[] buf;
  return is;
}


template <class T>
std::ostream& operator <<(std::ostream& os,const vector<T>& v)
{
    if (v.coeff==NULL) os<<vector<T>::zero<<"  ";
    for (int i=0;i<v.dim;i++) os<<v.coeff[i]<<"  ";
    return os;
}




} //end namespace val



#endif // VECTOR_H_INCLUDED
