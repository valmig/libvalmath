#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <initializer_list>
#include <val_basics.h>
#include <val_utils.h>
#include <error.h>




namespace val
{


// friend-functions of class matrix -----------------------------------------
template <class T> class matrix;
template <class T> class vector;
template <class T> matrix<T> operator*(const T&,const matrix<T>&);
template <class T> std::istream& operator >>(std::istream&,matrix<T>&);
template <class T> std::ostream& operator <<(std::ostream&,const matrix<T>&);
#ifdef VECTOR_H_INCLUDED
template <class T> const vector<vector<T> >& movefrommatrixtovector(matrix<T>& A,vector<vector<T> >& V,int r);
#else
template <class T> const vector<vector<T> >& movefrommatrixtovector(matrix<T>& A,vector<vector<T> >& V,int r=0);
#endif // VECTOR_H_INCLUDED

template <class T> const matrix<T>& movefromvectortomatrix(vector<vector<T> >&,matrix<T>&);
// ------------------------------------------------------------


template <class T>       // T Ring, ==0,==1,==-1,=0,=1,=-1, def.
class matrix
{
private:
	  T** coeff = nullptr;
	  int rows = 0;
	  int columns = 0;
	  void deletematrix();
	  static const T zero;
	  static matrix<T> char_to_matrix(const std::string& s); // converts std::string to matrix<T>
public:
	  matrix()  = default;
	  matrix(int m ,int n);           // Matrix with m rows and n columns. Will not be initialized.
	  explicit matrix(int n);             // quadratic matrix, not initialized.
	  matrix(const T& y,int m,int n);  // m x n Matrix, elements in it are  = y.
	  matrix(const matrix<T>&);
	  matrix(matrix<T>&&);
	  explicit matrix(std::initializer_list<std::initializer_list<T> > args);
	  ~matrix() {deletematrix();}
	  const matrix<T>& make_zero();     //   all elements  = T(0) 
      const matrix<T>& make_identity();  //  T(1) in the diagonal, all the others = T(0)
      int numberofrows() const {return rows;}
      int numberofcolumns() const {return columns;}
	  const matrix<T>& operator =(const matrix<T>&);
	  const matrix<T>& operator =(matrix<T>&&);
	  int iszero() const;
	  int isempty() const {return (coeff==NULL);}
	  int operator ==(const matrix<T>&) const;
	  int operator !=(const matrix<T>&) const;
	  matrix<T> operator +(const matrix<T>&) const;
	  matrix<T> operator -(const matrix<T>&) const;
	  const matrix<T>& operator +=(const matrix<T>&);
	  const matrix<T>& operator -=(const matrix<T>&);
	  matrix<T> operator *(const matrix<T>&) const;
	  const matrix<T>& operator *=(const matrix<T>&);
      const matrix<T>& operator *=(const T&);
	  const T* operator[](int) const;              // use operator() to avoid NULL-pointer assignment
	  T* operator[](int);                          
	  const T& operator()(int,int) const;
	  T& operator()(int,int);
	  const matrix<T>& swaprows(int,int); 
	  const matrix<T>& swapcolumns(int,int); 
	  const matrix<T>& append(const matrix<T>&);
	  const matrix<T>& append(matrix<T>&&);
	  friend matrix val::operator* <>(const T&,const matrix&);
	  friend int rnumber(const matrix<T>& A) {return A.rows;}
	  friend int cnumber(const matrix<T>& A) {return A.columns;}
	  friend const vector<vector<T> >& movefrommatrixtovector<T>(matrix<T>& A,vector<vector<T> >& V,int r);
	  friend const matrix<T>& movefromvectortomatrix<T>(vector<vector<T> >&,matrix<T>&);
      friend std::istream& operator >><T>(std::istream&,matrix<T>&);
      friend std::ostream& operator <<<T>(std::ostream&,const matrix<T>&);
};



template<class T>
const T matrix<T>::zero(val::zero_element<T>());


template <class T>
void matrix<T>::deletematrix()
{
    if (coeff==NULL) return;
    int i;
    for (i=0;i<rows;i++) delete[] coeff[i];
    delete[] coeff;
    coeff=NULL;
    rows=columns=0;
}


template <class T>
matrix<T>::matrix(int m,int n)
{
    if (m<=0 || n<=0) {
        coeff=NULL;
        rows=columns=0;
        return;
    }
    int i;
    coeff= new T*[rows=m];
    for (i=0;i<m;i++) coeff[i]=new T[columns=n];
}


template <class T>
matrix<T>::matrix(const T& y,int m,int n)
{
    if (m<=0 || n<=0) {
        coeff=NULL;
        rows=columns=0;
        return;
    }
    int i,j;
    coeff= new T*[rows=m];
    for (i=0;i<m;i++) {
        coeff[i]=new T[columns=n];
        for (j=0;j<n;j++) coeff[i][j]=y;
    }
}


template <class T>           
matrix<T>:: matrix(int n)
{
    if (n<=0) {
        coeff=NULL;
        rows=columns=0;
        return;
    }
    int i;
    coeff= new T*[rows=columns=n];
    for (i=0;i<n;i++) coeff[i]=new T[n];
}


template <class T>
matrix<T>:: matrix(const matrix<T>& A)
{
 if (A.coeff==NULL) {
     coeff=NULL;
     rows=columns=0;
     return;
 }

 int i,j;
 coeff= new T*[rows=A.rows];
 for (i=0;i<rows;i++) coeff[i]=new T[columns=A.columns];
 for (i=0;i<rows;i++)
     for (j=0;j<columns;j++) coeff[i][j]=A.coeff[i][j];
}


template <class T>
matrix<T>::matrix(matrix<T>&& A)
{
    coeff = A.coeff;
    rows=A.rows;columns=A.columns;
    A.coeff=NULL;
    A.rows=A.columns=0;
}


template <class T>
matrix<T>:: matrix(std::initializer_list<std::initializer_list<T> > args)
{
    coeff=NULL;
    rows=columns=0;
    const std::initializer_list<T> *V;
    const T *wert;

    rows = args.size();
    if (rows<=0) {rows=0;return;}
    // columns:
    for (V=args.begin();V!=args.end();V++) {
        columns = val::Max(columns,int(V->size()));
    }
    if (columns<=0) {columns=rows=0;return;}

    matrix<T> B(zero,rows,columns);

    *this = std::move(B);

    int i,j;


    for (i=0,V=args.begin();V!=args.end();V++,i++) {
        for (j=0,wert=V->begin();wert!=V->end();wert++,j++) coeff[i][j]=*wert;
    }
}


template <class T>
const matrix<T>& matrix<T>::make_zero()
{
    if (coeff==NULL) return *this;
    int i,j;

    for (i=0;i<rows;i++)
        for (j=0;j<columns;j++) coeff[i][j]=zero;
    return *this;
}


template <class T>
const matrix<T>& matrix<T>::make_identity()
{
    if (coeff==NULL) {
         Error::warning("\nWARNING: matrix<T>::make_identity: matrix==NULL, operation will be ignored!");
         return *this;
    }
    int i,j;
    for (i=0;i<rows;i++)
        for (j=0;j<columns;j++) {
            if (i==j) coeff[i][j]=T(1);
            else coeff[i][j]=zero;
        }
    return *this;
}



template <class T>
const matrix<T>& matrix<T>:: operator =(const matrix<T>& A)
{
 if (coeff==A.coeff) return *this;

 if (rows!=A.rows || columns!=A.columns) {
    matrix<T> help(A);
    val::swap(help.coeff,coeff);
    val::swap(help.columns,columns);
    val::swap(help.rows,rows);
    return *this;
 }
 int i,j;

 for (i=0;i<rows;i++)
     for(j=0;j<columns;j++) coeff[i][j]=A.coeff[i][j];
 return (*this);
}


template <class T>
const matrix<T>& matrix<T>:: operator =(matrix<T>&& A)
{
    if (coeff == A.coeff) return *this;
    if (coeff!=NULL) {
        deletematrix();
    }
    coeff=A.coeff;
    rows=A.rows;columns=A.columns;
    A.coeff=NULL;
    A.rows=A.columns=0;
    return *this;
}


template <class T>
int matrix<T>::iszero() const
{
    if (coeff==NULL) return 1;
    int i,j;
    for (i=0;i<rows;i++) {
        for (j=0;j<columns;j++) {
            if (coeff[i][j] != zero) return 0;
        }
    }
    return 1;
}


template <class T>
int matrix<T>:: operator ==(const matrix<T>& A) const
{
 if (coeff==NULL) return (A.coeff==NULL);

 int i,j,w=1;
 if (rows!=A.rows || columns!=A.columns) return 0;
 for (i=0;i<rows;i++)
     for(j=0;(j<columns && w);j++) w*=(coeff[i][j]==A.coeff[i][j]);
 return w;
}



template <class T>
int matrix<T>:: operator !=(const matrix<T>& A) const
{
 return !(*this == A);
}



template <class T>
matrix<T> matrix<T>:: operator +(const matrix<T>& A) const
{
 if (coeff==NULL) return matrix<T>(A);
 if (A.coeff==NULL) return matrix<T>(*this);

 int i,j,r=rows,c=columns;
 if (rows!=A.rows || columns!=A.columns) {
    Error::warning("\nWARNING: matrix<T>::operator +: matrices of different size!");
    if (rows>A.rows) r=A.rows;
    if (columns>A.columns) c=A.columns;
 }
 matrix<T> B(r,c);

 for (i=0;i<r;i++)
     for(j=0;j<c;j++) B.coeff[i][j]=coeff[i][j]+A.coeff[i][j];
 return B;
}


template <class T>
matrix<T> matrix<T>:: operator -(const matrix<T>& A) const
{
 if (coeff==NULL) {
     if (A.coeff==NULL) return matrix<T>(A);
     matrix<T> B(A);
     int i,j;
     for (i=0;i<B.rows;i++)
        for (j=0;j<B.columns;j++) B.coeff[i][j]*=T(-1);
     return B;
 }
 if (A.coeff==NULL) return matrix<T>(*this);

 int i,j,r=rows,c=columns;
 if (rows!=A.rows || columns!=A.columns) {
    Error::warning("\nWARNING: matrix<T>::operator -: matrices of different size!");
    if (rows>A.rows) r=A.rows;
    if (columns>A.columns) c=A.columns;
 }
 matrix<T> B(r,c);

 for (i=0;i<r;i++)
     for(j=0;j<c;j++) B.coeff[i][j]=coeff[i][j]-A.coeff[i][j];
 return B;
}


template <class T>
const matrix<T>& matrix<T>:: operator +=(const matrix<T>& A)
{
 if (A.coeff==NULL) return *this;
 if (coeff==NULL) {
    *this=A;
    return *this;
 }

 int i,j,r=rows,c=columns;
 if (rows!=A.rows || columns!=A.columns) {
    Error::warning("\nWARNING: matrix<T>::operator +=: matrices of different size!");
    if (rows>A.rows) r=A.rows;
    if (columns>A.columns) c=A.columns;
 }
 for (i=0;i<r;i++)
     for(j=0;j<c;j++) coeff[i][j]+=A.coeff[i][j];
 return *this;
}


template <class T>
const matrix<T>& matrix<T>:: operator -=(const matrix<T>& A)
{
 if (A.coeff==NULL) return *this;
 if (coeff==NULL) {
    *this=A;
    for (int i =0;i<rows;i++)
        for (int j=0;j<columns;j++) coeff[i][j]*=T(-1);
    return *this;
 }

 int i,j,r=rows,c=columns;
 if (rows!=A.rows || columns!=A.columns) {
    Error::warning("\nWARNING: matrix<T>::operator -=: matrices of different size!");
    if (rows>A.rows) r=A.rows;
    if (columns>A.columns) c=A.columns;
 }
 for (i=0;i<r;i++)
     for(j=0;j<c;j++) coeff[i][j]-=A.coeff[i][j];
 return *this;
}


template <class T>
matrix<T> matrix<T>:: operator *(const matrix<T>& A) const
{
 if (coeff==NULL || A.coeff==NULL) {
    matrix<T> B;
    return B;
 }
 int i,j,k,r=rows,c=A.columns,h=columns;
 if (columns!=A.rows) {
    Error::warning("\nWARNING: matrix<T>::operator *: matrices of wrong size!");
    if (columns>A.rows) h=A.rows;
 }
 matrix<T> B(r,c);

 for (i=0;i<r;i++)
     for (j=0;j<c;j++) {
        B.coeff[i][j]=zero;
        for (k=0;k<h;k++) B.coeff[i][j]+=coeff[i][k]*A.coeff[k][j];
     }
 return B;
}


template <class T>
const matrix<T>& matrix<T>:: operator *=(const matrix<T>& A)
{
 *this= (*this) * A;
 return *this;
}


template <class T>
const matrix<T>& matrix<T>:: operator *=(const T &x)
{
 if (coeff==NULL) return *this;
 int i,j;
 for (i=0;i<rows;i++)
     for (j=0;j<columns;j++) coeff[i][j]*=x;
 return *this;
}


template <class T>
const T* matrix<T>:: operator[](int i) const
{
 if (coeff==NULL) Error::error("\nError: matrix<T>::operator[]: NULL-POINTER!");
 if (i<0 || i>=rows) return coeff[0];
 return coeff[i];
}


template <class T>
T* matrix<T>:: operator[](int i)
{
 if (coeff==NULL) Error::error("\nError: matrix<T>::operator[]: NULL-POINTER!");
 if (i<0 || i>=rows) return coeff[0];
 return coeff[i];
}


template <class T>
const T& matrix<T>:: operator()(int i,int j) const
{
    if (i<0 || i>=rows || j<0 || j>= columns)
        return zero;
    return coeff[i][j];
}

template <class T>
T& matrix<T>:: operator()(int i,int j)
{
    if (i<0 || i>=rows || j<0 || j>= columns) {
		std::string msg="\n" + val::gettypename(*this);
		msg+=": ERROR: Index out of range!\nIndex i: "+ToString(i) + ". Rows: " + ToString(rows);
		msg+="\nIndex j: "+ToString(j) + ". Columns: " + ToString(columns)+"\n";
		Error::error(msg.c_str());
	}
    return coeff[i][j];
}


template <class T> 
const matrix<T>& matrix<T>:: swaprows(int i,int j)
{
 T* hilf;

 if (coeff==NULL) return *this;
 if (i<0 || i>=rows || j<0 || j>=rows) {
     Error::warning("\nWARNING: matrix<T>::swaprows: indices out of range!\nOperation will be ignored!");
     return *this;
 }

 hilf=coeff[i];
 coeff[i]=coeff[j];
 coeff[j]=hilf;
 return *this;
}


template <class T>
const matrix<T>& matrix<T>::swapcolumns(int i,int j)
{
    if (i<0 || i>=columns || j<0 || j>=columns)  {
     Error::warning("\nWARNING: matrix<T>::swapcolumns: indices out of range!\nOperation will be ignored!");
     return *this;
 }
    for (int k=0;k<rows;k++) val::swap(coeff[k][i],coeff[k][j]);
    return *this;
}

template <class T>
const matrix<T>& matrix<T>::append(const matrix<T>& A)
{
    if (A.coeff==nullptr) return *this;
    if (coeff==nullptr) {
        *this=A;
        return *this;
    }
    if (columns!=A.columns) {
        std::string msg= "\n" + gettypename(*this) + "::append(A) :ERROR: Different sizes of columns\n";
        msg+="\n columns : " + ToString(columns);
        msg+="\n A.columns: " + ToString(A.columns);
        Error::error(msg.c_str());
    }

    int m=rows + A.rows, n = columns,i,j;
    matrix<T> H(m,n);

    for (i=0;i<rows;++i) {
        for (j=0;j<n;++j) H.coeff[i][j] = coeff[i][j];
    }
    for (i=0;i<A.rows;++i) {
        for (j=0;j<n;++j) H.coeff[i+rows][j] = A.coeff[i][j];
    }
    val::swap(coeff,H.coeff);
    val::swap(columns,H.columns);
    val::swap(rows,H.rows);
    return *this;
}


template <class T>
const matrix<T>& matrix<T>::append(matrix<T>&& A)
{
    if (A.coeff==nullptr) return *this;
    if (coeff==nullptr) {
        *this=std::move(A);
        return *this;
    }
    if (coeff==A.coeff) return append(A);
    if (columns!=A.columns) {
        std::string msg= "\n" + gettypename(*this) + "::append(A) :ERROR: Different sizes of columns\n";
        msg+="\n columns : " + ToString(columns);
        msg+="\n A.columns: " + ToString(A.columns);
        Error::error(msg.c_str());
    }

    int m=rows + A.rows,i;
    T **h= new T*[m];

    for (i=0;i<rows;++i) {
        h[i]=coeff[i];
        coeff[i] = nullptr;
    }
    delete[] coeff;
    for (i=0;i<A.rows;++i) {
        h[i+rows] = A.coeff[i];
        A.coeff[i] = nullptr;
    }
    delete[] A.coeff;
    A.coeff=nullptr;
    A.rows=A.columns=0;
    coeff=h;
    rows=m;
    return *this;
}


template <class T>
matrix<T> operator *(const T &x,const matrix<T>& A)
{
    matrix<T> B(A);
    B*=x;
    return B;
}


template <class T>
matrix<T> matrix<T>::char_to_matrix(const std::string &s)
{
    int l=s.size(),first=0,i,j,k,nzeichen=0,m,n;
    std::string zeichen="";
    matrix<T> A;

    for (i=0;i<l;i++)
        if (s[i]!='\n' && s[i]!=' ' ) break;
    first=i;
    m=n=0;
    //Get number of columns n
    for (;i<l;) {
        if (i<l && s[i]=='\n') {break;}
        while (i<l && (s[i]!=' ' && s[i] != '\n' )) i++;
        n++;
        while (i<l && s[i]==' ') i++;
    }
    //Get number or rows m
    if (n==0) {m=0;return A;}
    for (i=first;i<l;) {
        if (s[i]!=' ' && s[i]!='\n') {
            while (i<l && s[i]!=' ' && s[i]!='\n') i++;
            nzeichen++;
        }
        while (i<l && (s[i]==' ' || s[i] =='\n')) i++;
    }
    m=nzeichen/n;
    if (nzeichen%n!=0) m++;

    A=val::matrix<T>(m,n);
    A.make_zero();

    i=j=0;
    for (k=first;k<l;) {
        if (s[k]!=' ' && s[k]!='\n') {
            while (k<l && s[k]!=' ' && s[k]!='\n') {
                zeichen+=s[k];
                k++;
            }
            A(i,j)=val::FromString<T>(zeichen);
            zeichen="";
            j++;
            if (j>=n) {j=0;i++;}
        }
        while (k<l && (s[k]==' ' || s[k] =='\n')) {
               k++;
        }
    }
    return A;
}


template <class T>
std::istream& operator >>(std::istream& is,val::matrix<T>& A)
{
  std::string line="",buf="";

  while (is) {
    line="";
    std::getline(is,line);
    if (line=="") break;
    buf+="\n"+line;
  }

  A=val::matrix<T>::char_to_matrix(buf);
  return is;
}


template <class T>
std::ostream& operator <<(std::ostream &os,const val::matrix<T> &A)
{
    int m=A.numberofrows(),n=A.numberofcolumns(),i,j;
    for (i=0;i<m;i++) {
        for (j=0;j<n;j++) os<<A(i,j)<<"  ";
        if (i != (m-1)) os<<std::endl;
    }
    return os;
}


} // end namespace val




#endif // MATRIX_H_INCLUDED
